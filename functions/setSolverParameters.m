function [ solverSetup , fileHandle ] = setSolverParameters( solver , p , grid , vel , data )
%SETSOLVERPARAMETERS Sets all parameters from user settings
%
%
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////



%% Set up solver
% (1) Choose functions for desired approximation accuracy
%   Same accuracy is used by all components of motion.
switch(solver.accuracy)
 case 'low'
  solverSetup.derivFunc = @upwindFirstFirst;
  solverSetup.integratorFunc = @odeCFL1;
 case 'medium'
  solverSetup.derivFunc = @upwindFirstENO2;
  solverSetup.integratorFunc = @odeCFL2;
 case 'high'
  solverSetup.derivFunc = @upwindFirstENO3;
  solverSetup.integratorFunc = @odeCFL3;
 case 'veryHigh'
  solverSetup.derivFunc = @upwindFirstWENO5;
  solverSetup.integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

% (2) Set up time approximation scheme (maximum CFL number defined by CFL_stabel*factorCFL)
%   stats defines if steps and times are reported
solverSetup.integratorOptions = odeCFLset('factorCFL', solver.factorCFL , 'stats', solver.stats );

% (3) Set up convective motion (only basics here; the rest see below!)
convectiveFunc = @termConvection;
convectiveData.dim = solver.dim;
convectiveData.grid = grid;
convectiveData.derivFunc = solverSetup.derivFunc;
% convectiveData.integratorFunc = solver.integratorFunc;
convectiveData.p = p;

% (4) Set up basic motion in the normal direction.
normalFunc = @termNormal;
normalData.grid = grid;
normalData.speed = -p.s_l_u;
normalData.derivFunc = solverSetup.derivFunc;

% (5) Set up curvature motion if desired
if solver.curvature
  curvatureFunc = @termCurvature;
  curvatureData.grid = grid;
  curvatureData.b = p.marksteinLength * p.s_l_u;
  if strcmp(p.flameType,'slit')
    % Compute curvature slit geometry
    curvatureData.curvatureFunc = @curvatureSecond;
  elseif strcmp(p.flameType,'conical')
    % Compute curvature conical geometry
    curvatureData.curvatureFunc = @curvatureSecond_conical;
  end
end



%% Set up post time step functions
% (1) Add post time step function which reinitializes the G-field if desired
if solver.reinit.do 
  % Do reinitialize level-set
  solverSetup.integratorOptions = odeCFLset(solverSetup.integratorOptions, 'postTimestep', {@reinitializeLS,@writeInOut2File});
  
else
  % Do not reinitialize level-set (not recommended!)
  solverSetup.integratorOptions = odeCFLset(solverSetup.integratorOptions, 'postTimestep', {@writeInOut2File});
  
end

% (2) Options and data for all post time step function
solverSetup.schemeData.solver = solver;
solverSetup.schemeData.counter = 1; % counts how many steps have been performed
solverSetup.schemeData.transientDef = vel.transientDef;

% (3) Open files for output
fileNameA = [p.run_output_folder,filesep,'Surface.out'];
fileNameV = [p.run_output_folder,filesep,'V_ref.out'];

if strcmp(solver.initial,'default') || strcmp(solver.initial,'useInit')
  % Open file for output surface
  solverSetup.schemeData.writeInOut.fileA = fopen(fileNameA,'w');
  fprintf(solverSetup.schemeData.writeInOut.fileA,'time [s] , A(t)');
  
  % Open file for output reference velocity
  solverSetup.schemeData.writeInOut.fileV = fopen(fileNameV,'w');
  fprintf(solverSetup.schemeData.writeInOut.fileV,'time [s] , v_ref(t)');
  
elseif strcmp(solver.initial,'resume')
  % Open and check file for output surface
  solverSetup.schemeData.writeInOut.fileA = openAndCheckFile( solver, fileNameA  );
  
  % Open and check file for output reference velocity
  solverSetup.schemeData.writeInOut.fileV = openAndCheckFile( solver, fileNameV  );
  
end

% Returns file handles separately for convenience
fileHandle = {solverSetup.schemeData.writeInOut.fileA , solverSetup.schemeData.writeInOut.fileV };



%% Set up visualization/ output
% Set up output intervalls and plotting
if solver.outSteps <= 0
  % Infinity plot period: plotting only at beginning and end (tmax)
  solverSetup.tOut = inf;
else
  % Period at which intermediate plots should be produced.
  solverSetup.tOut = (solver.tMax - solver.t0) / (solver.outSteps );
end
  
% Stepper takes only one step: Plot at every single time step (SLOW!).
if (solver.singleStep)
  solverSetup.integratorOptions = odeCFLset(solverSetup.integratorOptions, 'singleStep', 'on');
end



%% Set up convectiveData

% (1) Prepare inlet signal if FOAM list is chosen
if vel.useFOAM_list 
  % First time step of given time series is at desired t0 (solver) + t_transient0 (velocity)
  ListFOAM2GFLAME( vel.fileFOAM , vel.path2File ,'resample',vel.dx,...
    'NoshiftT0',vel.t_transient0 ,'solverT0',solver.t0+vel.t_transient0);
end


% (2) sweep settings/ preparations
if ( strcmpi(vel.transientDef,'sweep') || strcmpi(vel.transientDef,'broad') ) 
  
  % Check if TFDtools with siCreateSineSweep() are added to path
  if ~ (exist('siCreateSineSweep.m','file')==2) && strcmpi(vel.transientDef,'sweep')
    error('Please add TFDtools to your Matlab path!')
  end
  
  % Generate sine sweep and write uInlet and sweep_data.mat file to output folder
  generateSpecialInputSignalGFLAME( vel , p );
  
end


% (3) TaX settings
if strcmpi(vel.transientDef,'TaX')
  % Check if TaX is add to path
  if ~ (exist('tax','class')==8)
    error('Please add TaX to your Matlab path!')
  end
  
  % Check if TaX model exists
  if exist(vel.path2sys,'file')~=2
    error(['Specified TaX model could not be found in <',vel.path2sys,'>!'])
  end
  
  % load TaX model
  tmp = load(vel.path2sys);
  % First field is intepreted as TaX model
  tmp_n = fieldnames(tmp);
  convectiveData.sys_acoustic = tmp.(tmp_n{1});
  if ~isa(convectiveData.sys_acoustic,'tax')
    error(['Invalid TaX model specified in <',vel.path2sys,'>!'])
  end
  
  % Mean flame surface area
  convectiveData.A_mean = vel.A_mean;
  
  % Acoustic CFL number
  convectiveData.CFL_ac = vel.CFL_ac;
  
  % Compute acoustic time step from CFD number
  convectiveData.Ts_ac = CFLtoTs(convectiveData.sys_acoustic, convectiveData.CFL_ac);
  
  % Initialize uref, A_surf and time vector (length of number of elements in x1 direction)
  nFac = 30; % Ensure that time series is long enough! Factor introduced by Markus Brandl 2016 (BA)
  if strcmp(solver.initial,'resume')
    % If simulation was resumed, then read in past time steps
    t0 = solver.t0 - ( solver.x1Lim(2) - solver.x1Lim(1) ) / vel.meanFlowSpeed;
    t1 = solver.t0;
    [ t_vec , uref_vec , A_surf_vec ] = getLastNTimeStepsFromFile( fileNameA , fileNameV ,...
      t0 , t1 , nFac*grid.N(1) );
    
  else
    % If not resumed, then initialise all fields
    t_vec = ones(nFac*grid.N(1),1) * solver.t0;
    uref_vec = ones(nFac*grid.N(1),1) * vel.meanFlowSpeed;
    A_surf_vec = ones(nFac*grid.N(1),1) * convectiveData.A_mean;
    
  end
  % Write to scheme data
  convectiveData.uref   = uref_vec;
  convectiveData.A_surf = A_surf_vec;
  convectiveData.t_vec  = t_vec;

else
  % In not Tax mode only store last value for uref, A_surf and t_vec
  convectiveData.uref   = 0;
  convectiveData.A_surf = 0;
  convectiveData.t_vec  = 0;
  
end


% (4) Open external velocity data (if desired!)
if strcmp(vel.velModel , 'CFD_FOAM')
  % Load velocity data
   load([vel.outputDir,filesep,'U.mat'],'dataFOAM');
   convectiveData.externalVelData = dataFOAM;
end


% (5) Create velocity object by VORTEX (if desired)
if strcmp(vel.velModel , 'VORTEX')
  % Add velocityFiled object to schemeData
  convectiveData.velocityObject = VelocityField();
  % Define starting time
  convectiveData.velocityObject.timeSim = solver.t0;
  % Define bulk velocity
  convectiveData.velocityObject.potentialFlowField.bulkFlow = {vel.meanFlowSpeed,0};
end


% (6) Add data from velocity model settings to convectiveData (except logfile)
fieldNames = fieldnames(vel);
for ii = 1:length(fieldNames)
  if ~strcmp(fieldNames{ii},'logfile')
    convectiveData.(fieldNames{ii}) = vel.(fieldNames{ii});
  end
end


% (7) Initialize velocity field.
[convectiveData.velocityField, schemeData] = feval( convectiveData.velocity , solver.t0 , data , convectiveData ); %Axel (schemeData)

if strcmp(convectiveData.velModel,'FirstPrincipleBased')     %Axel
  convectiveData.FPB = schemeData.FPB;
end
clear schemeData;


%% Combine components of motion.
if solver.curvature
  % If there is a nonzero curvature contribution to speed.
  solverSetup.schemeFunc = @termSum;
  solverSetup.schemeData.innerFunc = { normalFunc; curvatureFunc; convectiveFunc };
  solverSetup.schemeData.innerData = { normalData; curvatureData; convectiveData };
else
  % Otherwise ignore curvature.
  solverSetup.schemeFunc = @termSum;
  solverSetup.schemeData.innerFunc = { normalFunc; convectiveFunc };
  solverSetup.schemeData.innerData = { normalData; convectiveData };
end


end

