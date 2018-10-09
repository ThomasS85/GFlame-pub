function [ s ] = settings_Solver_Example( p )
%SETTINGS_SOLVER Solver Settings for G-equation

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% General
% Dimensionality of the problem (currently only 2 supported!!)
s.dim = 2;  
% include diffusion due to curvature? -> Set Markstein length if so!
s.curvature = 0;
% Load initial data or start with straight default (line)? 
%     'resume':   - will load a file 'init.mat' in directory specified by
%                   resume2case (see below!)  
%                 - This init.mat file should contain data and GFcase
%                 - Only solver settings can be changed (except grid related
%                   ones)!
%     'default':  - Starts simulation with a straight line as initial
%                   condition
%     'useInit':  - Uses a file 'init.mat' stored at the case root
%                   directory which contains a field 'data' as IC
s.initial = 'default';
% Number of the folde which contains the case which shall be resumed
s.resume2case = 1;
% Maximum number of output folder. If maximum is reached the first folder will 
% be overwritten! Set to '0' to deactivate this behaviour!
s.purgeWrite = 0;


%% Integration parameters% Compute flame surface area or heat release? 'area' or 'heatRelease'
s.outputY = 'area';
s.t0 = 0;             % start time
s.tMax = 0.01;        % end time
% How close (relative) do we need to get to tMax to be considered finished?
s.small = 100 * eps;
% Time approximation scheme: maximum CFL number; defined by CFL_stabel*factorCFL
%   Use smaller value if curvature is switched off (otherwise flame surface time series might be noisy)!
s.factorCFL = 0.5;
% Should steps and times be reported to command window?
s.stats = 'on';

% Order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2. (Best
%                                 compromise between time and accuracy)
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
s.accuracy = 'medium';


%% Spatial Discretisation
% solver.dxi = 2e-4;
% s.dxi = 4e-4;   % Cuquel
% s.dxi = 5e-5;   % Kornilov
s.dxi = 2e-4; %Axel


%% Settings reinitialization (G-solver)
% Reinitialisation prevents the gradient of the G-gield close to the zero iso-line do become smaller which
% would reduce accuracy. The Reinitialisation needed for the surface integration is not affected by this
% settings here (see settings below!).
s.reinit.do = 1;               % Do reinitialisation after timestep?
s.reinit.intervalStep = 20;    % Do reinitialization every n steps (best but slowest: 1)
s.reinit.accuracy = 'medium';  % Same options as above (in tests medium was better than veryHigh and faster)
% Time at which to halt the reinitialization iteration (default = max(grid.max - grid.min)).
%  If tMax < 0, it is interpreted as the number of CFL limited reinitialization timesteps 
%  to take: number of steps = -round(tMax).
s.reinit.tMax = s.dxi*10;   % about 10 cells are reinitialized
% If the average update of nodes drops below errorMax * max(grid.dx), then assume that 
%   reinitialization has converged and return early.
s.reinit.errorMax = 1e-3;


%% Settings surface integration 
% What method to use?  
%     'deltaFun' = use delta function  or  'contour' = integrate along contour line (recommended)
s.SurfInt_method = 'contour';
% order of the delta function: 1 or 2 
s.orderDelta = 2; % order 1 is recommended for low resolution setups since it is more reliable here
% Evaluate surface and reference velocity every n steps:
s.intervalWriteOI = 5;
% If curvature is 1: Define standard deviation of Gaussian filter for surface intergation (0.8 recommended). Set
% to -1 to switch it off.
s.sigmaG = 0.8;


%% Settings reinitialization (Surface Integration using a delta function)
% The used delta function requires the G-field to be a signed distance function close to the zero iso-line.
% Thus, reinitialization is neccessary! It is recommended to use high accuracy here to remove noise!
s.deltaFunReinit.accuracy = 'veryHigh';  
s.deltaFunReinit.tMax = s.dxi*4;   
% If the average update of nodes drops below errorMax * max(grid.dx), then assume that 
%   reinitialization has converged and return early.
s.deltaFunReinit.errorMax = 1e-5;


%% Output settings
% How many intermediate outputs to produce? '0': No output, '1': Output at
% beginning and end. Output can be plot or data which is written to disc
s.outSteps =  20; 
% Output at each timestep (overrides outSteps)? '1': yes '0': no
s.singleStep = 0;     

% Should be plotted (a gif is generated!)?
s.doPlot = 1;
s.exportPlot = 'gif'; % ''-> No export; 'gif'-> A gif is generated; 'png' -> png files are saved
% Delete previous plot before showing next?
s.deleteLastPlot = 1;
% What kind of display?
s.displayType = 'contour';
% s.displayType = 'surf';
% Define position size of figure: [x1 x2 widthX1 widthX2] in pixel
s.figurePositionProps = [0,100,1000,1000];

% Define what should be plotted
s.plotVelPcolor = 1;               % Do Velocity pcolor plot?
s.pcolorLim = [0.5,1.5];           % Limits for pcolor plot
s.pcolorNTicks = 5;                % Number of ticks for pcolor plot
s.plotVelQuiver = 1;               % Do Velocity quiver plot?
s.quiverOptions = [ 5 , 2 , 0.3 ]; % [dPlotX1,dPlotX2,scale]
s.plotGrid = 0;                    % Do plot of grid?
s.doPlotBurner = 2;                % Should burner be plotted? Set to 2 if also confinement should be included
s.lablesOff = 1;                   % Plots no lables 

% Should data be written to disc (same intervall as plotting?)
s.doWriteData = 0;


%% Domain settings and boundary conditions
% Minima/ Maxima of geometry in x1 direction (for grid generation, plotting,etc) WITHOUT flame lift-off!
% s.x1Lim = [ 0 1.3*p.H_flame+p.liftOff ];
% s.x1Lim = [ 0 , 40e-3 ];  % Cuquel
% s.x1Lim = [ 0 , 8e-3 ];  % Kornilov
% s.x1Lim = [ -30e-3 , 30e-3 ]; % Axel duct
s.x1Lim = [ 0 , 30e-3 ]; % Axel backwardFacingStep

% Boundary Condition in x1 direction:
%   radial    - Holds the flame at the radial walls. Domain width is
%                p.R_flame
%   axial     - Holds the flame at inlet at radius p.R_i. Domain width
%                is p.R_a!
s.holdFlame = 'axial';

end

