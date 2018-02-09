clear 
close all
clc

% Run GFLAME 0.1 

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Add GFLAME to path
addpath(genpath('/home/steinbacher/Documents/MATLAB/GFLAME-pub'));

% For specific functions you might also need Tax/sss and TFDTools added to your path!


%% define functions for settings
Fun4GenSettings = 'settings_Example';
Fun4SolSettings = 'settings_Solver_Example';
Fun4VelSettings = 'settings_VelocityModel_Example';

%% Should solver be started?
doSolve = 1;

%% Set up case
[ GFcase  ] = setUpCase( Fun4GenSettings , Fun4SolSettings , Fun4VelSettings );


%% Integrate level set
if doSolve
  % Integrate Level Set
  disp('Integrating Level set...')
  [ data , t_end , GFcase ] = integrateLevelSet( GFcase );
  disp('Finished!')

  % Post process data
  % Posprocess sweep signal
  if strcmp(GFcase.vel.transientDef,'sweep')
    postProcessGFLAMESweep( GFcase );
  end
  % Posprocess broad band data
  if strcmp(GFcase.vel.transientDef,'broad')
    identifyGFLAMEOutputData( GFcase );
  end
  
  % Finish execution (necessary to close logfile!)
  logging( GFcase , 4  )
end



%% plotSurfaceVelocityData( GFcase )
% deriveSurfaceSignal( GFcase , 'plot' , 'smoothing' , 40 , 'plotRef');
% plotSurfaceVelocityData( GFcase )

%% Plot solution
% visualizeGrid(GFcase.grid);
% plotGFLAMEOutputData( GFcase );

%% Extract iso line
% [ C ] = extractIsoLine( GFcase ,GFcase.p.caseNumber , 0.04 , 0,'changeOriginX1',31e-3);


%% Identify flame
% For these functions TFDtools are required!
%  close all

% (1) Slit Flame
% FIR filter
%  identifyGFLAMEOutputData( GFcase ,'t0',0,'case',1,'timeDelay',0,'maxFreq',1200,'compOrders',[15,10],...
%     'sweepCase',2) 

% ARX filter
% identifyGFLAMEOutputData( GFcase ,'t0',0,'case',1,'timeDelay',0,'maxFreq',1200,'idModel','ARXregul',4,...
%   'compOrders',[4,3],'sweepCase',2) 
