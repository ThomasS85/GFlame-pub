function [ GFcase  ] = setUpCase( Fun4GenSettings , Fun4SolSettings , Fun4VelSettings , varargin )
%SETUPCASE Sets all parameters for a GFLAME case
%
% 17.04.15: added Feature to pass parameters to setUpCase via varargin and
%           setUserPaameters()
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Parse varargin for settings
% Should initial G-field be transformed to a signed distance function?
ind = find(strcmpi(varargin,'noInit'),1);
if ~isempty(ind)
  % Do not initialize it to signed distance
  makeSignedDist = 0;
else
  % Do initialize it to signed distance (default)
  makeSignedDist = 1;
end


%% Get user parameters
% Get parameters set by users

% (a) Get physical and geometric parameters
% Get general user parameters 
p = eval([Fun4GenSettings,'()']);
if nargin > 3
  p = setUserParameters( p , varargin );
end

% Calculate derived or fixed parameters
p = setDerivedParamters(p);

% Get folder of main case (for subfunctions)
p.caseRootDir = pwd;

% (b) Get solver parameters
solver = eval([Fun4SolSettings,'(p)']);
if nargin > 3
  solver = setUserParameters( solver , varargin );
end
% Respect flame lift-off
solver.x1Lim = solver.x1Lim + p.liftOff;

% Check if discretisation is too small
if max(p.R_i / solver.dxi, p.H_flame / solver.dxi) > 1e3
  error('Discretisation too small! Please use different spatial discretisation s.dxi.)')
end

% parse solver settings and set values for BC and limits in x1/2 direction
if strcmpi(p.geom,'Vinv')
  
  % Inverse V flames
  if strcmp(p.domainType,'sym')
    if strcmp(solver.holdFlame,'radial')
      % Limits in x2 direction
      solver.x2Lim = [ 0 p.R_flame ];
      % BCs
      solver.boundaryCon{1} = @addGhostExtrapolate; % x1-direction
      solver.boundaryCon{2} = @addGhostVectorDirichletSym; % x2-direction
    elseif strcmp(solver.holdFlame,'axial')
      % Limits in x2 direction
      solver.x2Lim = [ 0 p.R_flame+p.R_flame*0.1 ];
      % BCs
      solver.boundaryCon{1} = @addGhostVectorDirichlet_Extrapolate; % x1-direction
      if solver.curvature
        % Use Symmetry BC at center (lower noise for surface evaluation!)
        solver.boundaryCon{2} = @addGhostExtrapolateSym;  % x2-directionelse
      else
        % Do not use symmetry at center since no curvature!
        solver.boundaryCon{2} = @addGhostExtrapolate; % x2-direction
      end
      
    else
      error('Unknown flame holder chosen!')
    end
    
  elseif strcmp(p.domainType,'full')
    if strcmp(solver.holdFlame,'radial')
      % Limits in x2 direction
      solver.x2Lim = [ -p.R_flame p.R_flame ];
      % BCs
      solver.boundaryCon{1} = @addGhostExtrapolate; % x1-direction
      solver.boundaryCon{2} = @addGhostVectorDirichlet; % x2-direction
    elseif strcmp(solver.holdFlame,'axial')
      % Limits in x2 direction
      solver.x2Lim = [ -(p.R_flame+p.R_flame*0.1) (p.R_flame+p.R_flame*0.1) ];
      % BCs
      solver.boundaryCon{1} = @addGhostVectorDirichlet_Extrapolate; % x1-direction
      solver.boundaryCon{2} = @addGhostExtrapolate; % x2-direction
    else
      error('Unknown flame holder chosen!')
    end
    
  else
    error('No valid domain type chosen (settings file)!')
  end
  
  
elseif strcmpi(p.geom,'V') || strcmpi(p.geom,'M')
  % V flame
  if strcmp(p.domainType,'sym')
    % Limits in x2 direction
    solver.x2Lim = [ 0 p.R_a ];
    % BCs
    solver.boundaryCon{1} = @addGhostVectorDirichlet_Extrapolate; % x1-direction
    if solver.curvature
      % Use Symmetry BC at center (lower noise for surface evaluation!)
      solver.boundaryCon{2} = @addGhostExtrapolateSym;  % x2-directionelse
    else
      % Do not use symmetry at center since no curvature!
      solver.boundaryCon{2} = @addGhostExtrapolate; % x2-direction
    end
    
  elseif strcmp(p.domainType,'full')
    % Limits in x2 direction
    solver.x2Lim = [ -p.R_a p.R_a ];
    % BCs
    solver.boundaryCon{1} = @addGhostVectorDirichlet_Extrapolate; % x1-direction
    solver.boundaryCon{2} = @addGhostExtrapolate; % x2-direction
    
  else
    error('No valid domain type chosen (settings file)!')
  end
  
else
  error('Unknown flame geometry!')
end


%% Set up case
% So far all user settings have been read in and were (pre-) parsed, now
% based on these settings all parameters are derived
% Initialize
caseName = '';
outputDir = ''; 
nameFolds2 = [];

if strcmp(solver.initial,'default') || strcmp(solver.initial,'useInit')
  % Start a new run
  
  % Generate the grid
  [ myGrid ] = generateGrid( solver );
  
  % Set initial conditions
  [ data ] = setInitialConditions( myGrid, p , solver , makeSignedDist );
  
  % Set boundary conditions
  [ myGrid] = setBoundaryConditions( myGrid, p , solver , data );
  
  
  % Create case folder if there is none
  if ~(exist(p.caseFolder,'dir') == 7)
    mkdir(p.caseFolder)
  end
  
  % Define folder for output
  d = dir(p.caseFolder);
  isub = [d(:).isdir]; %# returns logical vector indicating all directories
  nameFolds = {d(isub).name}'; % write out all directories to cell aray of cells
  nameFolds = regexp(nameFolds,'^\d+$','match'); % delete all entries but those with numeric data
  nameFolds(cellfun('isempty',nameFolds)) = []; % delete empty entries
  nameFolds = cellfun(@cell2mat,nameFolds,'UniformOutput',false); % convert to array of strings
  nameFolds = cellfun(@str2num,nameFolds); % Convert to vector of numerics
  
  % Deside how output folder should be called
  if isempty(nameFolds)
    % If no case exits in that folder yet
    p.caseNumber = 1;
  else
    p.caseNumber = max(nameFolds)+1;
  end
  p.run_output_folder = [p.caseFolder,'/',num2str( p.caseNumber )];
  % check purge write and if folder is empty
  if length(nameFolds) >= solver.purgeWrite && solver.purgeWrite ~=0
    % if purge Write is active, delete smallest folder
    rmdir( [p.caseFolder,'/',num2str( min(nameFolds) )] ,'s' );
  elseif isempty(nameFolds)
    p.run_output_folder = [p.caseFolder,'/',num2str( 1 )];
  end
  % Create folder for output
%   mkdir(p.run_output_folder);
  
elseif strcmp(solver.initial,'resume')
  % Resumes to a previous run
  % define path 2 case to which we shall resume
  p.run_output_folder = [p.caseFolder,'/',num2str( solver.resume2case )];
  
  if exist(p.run_output_folder,'dir')~=7
    error(['Could not load case ',num2str( solver.resume2case),' since it does not exist.'])
  end
  
  % Find latest timestep
  d = dir(p.run_output_folder);
  isub = ~[d(:).isdir]; %# returns logical vector indicating all not directories
  nameFolds = {d(isub).name}'; % write out all files to cell array of cells
%   nameFolds = regexp(nameFolds,'^\d+\.\d+','match'); 
  nameFolds = regexp(nameFolds,'^\d+\.{0,1}\d*\.{0,1}','match'); % delete all entries but those with numeric data
%   nameFolds = regexp(nameFolds,'^\d+','match'); % delete all entries but those with numeric data
  nameFolds(cellfun('isempty',nameFolds)) = []; % delete empty entries
  nameFolds = cellfun(@cell2mat,nameFolds,'UniformOutput',false); % convert to array of strings
  nameFolds = cellfun(@(s) s(1:end-1) , nameFolds ,'UniformOutput', false ); % remove the last dot
  nameFolds = cellfun(@str2num,nameFolds); % Convert to vector of numerics  
  latestT = max(nameFolds);
  
  % Now resume to this time step
  if exist([p.run_output_folder,'/',num2str(latestT,15),'.mat'],'file') == 2
    % load data (g, initial conditions and boundary conditions)
    load([p.run_output_folder,'/',num2str(latestT,15),'.mat'],'GFcase','data')
    
    % Only use certain information from settings files to prevent errors
    
    % data is loaded automatically!
    
    % Load general settings and save actual
    p_tmp = p;
    p = GFcase.p;
    % Get actual folder of main case (for subfunctions) 
    p.caseRootDir = pwd;
    % Keep some data
    p.caseNumber = solver.resume2case;
    p.run_output_folder = p_tmp.run_output_folder;
    
    % Load grid data
    myGrid= GFcase.grid;
    
    % Load solver settings 
    % Continue simulation with desired time span
    tspan = solver.tMax - solver.t0;
    % Define t0 and tMax so that simulation continues with the desired time
    % and time span (tmax-t0)
%     solver.t0 = GFcase.solver.tMax;
    solver.t0 = latestT; % Use latest time step to continue simulation!
    solver.tMax = solver.t0 + tspan;
    % Overwrite dimension
    solver.dim = GFcase.solver.dim ;
    % Overwrite grid data with loaded data
    solver.dxi = GFcase.solver.dxi;
    solver.x1Lim = GFcase.solver.x1Lim;
    
    % Velocity settings are NOT loaded!
     
  else
    % Show warning since there is no GFLAME computation to continue. 
    warning(['No file <',[num2str(latestT),'.mat'],'>',' with initial conditions found in output-case folder! Ignore this warning if you only want to load CFD data.'])
    
    % If no data could be loaded initialise grid and data like in default
    % mode:
    % Generate the grid
    [ myGrid] = generateGrid( solver );  
    % Set initial conditions
    [ data ] = setInitialConditions( myGrid, p , solver );   
    % Set boundary conditions
    [ myGrid] = setBoundaryConditions( myGrid, p , solver , data );
  end
  
  % Look for loaded CFD data
  d = dir(p.run_output_folder);
  isub = [d(:).isdir]; %# returns logical vector indicating all not directories
  nameFolds2 = {d(isub).name}'; % write out all folders to cell array of cells
  nameFolds2 = regexp(nameFolds2,'((FOAM)|(VTK))_.*','match'); % delete all entries but those with numeric data
  nameFolds2(cellfun('isempty',nameFolds2)) = []; % delete empty entries
  nameFolds2 = cellfun(@cell2mat,nameFolds2,'UniformOutput',false); % convert to array of strings
  if ~isempty(nameFolds2)
    if length(nameFolds2)>1
      error('Only one folder with CFD data is allowed per case folder. Remove additional folders!')
    end
    caseName = nameFolds2{1};
    outputDir = [p.run_output_folder,filesep,caseName];
  end
  
  
else
  error('Non valid type of initial condition!')
  
end


%% Velocity field
% Get user parameters for velocity field
vel = eval([Fun4VelSettings,'(p)']);
if nargin > 3
  vel = setUserParameters( vel , varargin );
end

% Check if chosen velocity model is applicable
if strcmp(vel.velModel ,'convectiveConfined') || strcmp(vel.velModel ,'convectiveIncompConfined')
  if p.Cr < sqrt( 1 - (p.E-1)/p.E * cos(p.alpha) )
    % Set velocity model to unconfined if Cr is too small
    if strcmp(vel.velModel ,'convectiveConfined')
      vel.velModel = 'convective';
    else
      vel.velModel = 'convectiveIncomp';
    end
    warning('Confined velocity model not applicable since confinement ratio too small! Switched to conventional convective velocity model!')
  end
end

% Set parameters velocity field
vel = setVelocityParameters( vel , solver );

% Overwrite case name if case was resumed and CFD data were found
if ~isempty(nameFolds2)
  vel.caseName = caseName;
  vel.outputDir = outputDir;
  vel.path2Case = 'unknown';
end

%% etc
% If sweep is chosen in velocity settings, then compute final time
if strcmpi(vel.transientDef,'sweep')
  [~ , t_vec , f_orig ] = generateSpecialInputSignalGFLAME( vel , p , 'NoOutput' );
  solver.tMax = t_vec(end); 
  vel.f_min = f_orig(1);
  vel.f_max = f_orig(end);
  
elseif strcmpi(vel.transientDef,'broad')
  [~ , t_vec , f_orig ] = generateSpecialInputSignalGFLAME( vel , p , 'NoOutput' );
  solver.tMax = t_vec(end);
  vel.f_maxBB = f_orig;
  
end


%% Define output of logfile
% Define log file output
[ p,solver,vel ] = defineLogfileOut( p,solver,vel );
% Set indicator if logfile has already been initialized to zero
p.logFileInit = 0;

% % Get filename for logfile
% p.logFileName = [p.run_output_folder,'/log.out'];
% ii = 1;
% while exist(p.logFileName,'file')
%   p.logFileName = [p.run_output_folder,'/log',num2str(ii),'.out'];
%   ii = ii + 1;
% end


%% Output (wrap data)
GFcase.data = data;
GFcase.p = p;
GFcase.solver = solver;
GFcase.grid= myGrid;
GFcase.vel = vel;




end

