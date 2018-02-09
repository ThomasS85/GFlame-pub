function h = plotGFLAMEOutputData( GFcase , varargin )
%PLOTGFLAMEOUTPUTDATA plots output of GFLAME or creats gif
%
% Press
%     ->  : move on time step further
%     <-  : go back on step
%     ESC : exit
%     l   : jump to last time step
%     f   : Jump to first time step
%
%% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (tsteinbacher@tfd.mw.tum.de).    //
% // Created, 02.03.2015 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

change2caseRootDir( GFcase )

%% Parse varargin
% Create gif file?
ind = find(strcmpi(varargin,'makeGIF'),1);
if ~isempty(ind)
  % User input
  makeGIF = 1;
else
  % Default: No
  makeGIF = 0;
end

% Should case be change?
ind = find(strcmp(varargin,'case'),1);
if ~isempty(ind)
  % Change run output folder
  GFcase.solver.resume2case = varargin{ind+1};
  GFcase.p.caseNumber = varargin{ind+1};
else
  % Load case specified by case number
  GFcase.solver.resume2case = GFcase.p.caseNumber;
end

% Change dir to case root
change2caseRootDir( GFcase )

%% Get info about data
% define where to load data from
caseNumber = GFcase.p.caseNumber;
run_output_folder = ['.',filesep,num2str( caseNumber )];

% Find available timesteps
d = dir(run_output_folder);
isub = ~[d(:).isdir]; %# returns logical vector indicating all not directories
nameFolds = {d(isub).name}'; % write out all files to cell aray of cells
% nameFolds = regexp(nameFolds,'^\d+\.\d+','match'); % delete all entries but those with numeric data
nameFolds = regexp(nameFolds,'^\d+\.{0,1}\d*\.{0,1}','match'); % delete all entries but those with numeric data
nameFolds(cellfun('isempty',nameFolds)) = []; % delete empty entries
nameFolds = cellfun(@cell2mat,nameFolds,'UniformOutput',false); % convert to array of strings
nameFolds = cellfun(@(s) s(1:end-1) , nameFolds ,'UniformOutput', false ); % remove the last dot
% nameFolds = cellfun(@str2num,nameFolds); % Convert to vector of numerics


%% Plot data: Loop over time
% initialise plotting
f = figure('visible','on','color','w');
load([run_output_folder,filesep,nameFolds{1},'.mat'])
% unwrap data
grid = GFcase.grid;
solver = GFcase.solver;
solverSetup = GFcase.solverSetup;
p = GFcase.p;
% Set figure properties by solver settings
set(gcf,'Position',solver.figurePositionProps)
% Create visualization.
if isfield(solverSetup,'schemeData')
  visualizeVelocityField( solver , solverSetup , grid )
end
h = visualizeLevelSet(grid, GFcase.data, solver.displayType, 0, 'Init');
set(h,'Color','r','LineWidth',4)
% set(h,'linewidth',3)

if ~makeGIF
  waitforbuttonpress
else
  % Get filename for gif file
  filename = [run_output_folder,'/XiFlame_',p.flameType,'.gif'];
  
end

% loop over time
tt = 1;
loopTime = 1;
while loopTime
  figure(f);
  % Load data
  load([run_output_folder,filesep,nameFolds{tt},'.mat'])
  % unwrap data
  grid = GFcase.grid;
  solver = GFcase.solver;
  solverSetup = GFcase.solverSetup;
  % delete old visualisation
  cla( get(f,'CurrentAxes') )
  % Create new visualization.
  if isfield(solverSetup,'schemeData')
    visualizeVelocityField( solver , solverSetup , grid )
  end
  h = visualizeLevelSet(grid, data, solver.displayType, 0, [ 't = ', nameFolds{tt} ]);
  set(h,'Color','r','LineWidth',4)
%   set(h,'linewidth',3)
  
  if ~makeGIF
    % Interactive plotting
    waitforbuttonpress
    userInput = double( get(f,'CurrentCharacter') );
    
    if userInput == 27
      % Exit function if ESC is pressed
      loopTime = 0;
    elseif userInput == 29
      % Arrow key right: move one step further
      tt = min( tt + 1 , length(nameFolds) );
    elseif userInput == 28
      % Arrow key left: move one step back
      tt = max( tt - 1 , 1 );
    elseif userInput == 108
      % l: jump to last time step
      tt = length(nameFolds);
    elseif userInput == 102
      % f: jump to first time step
      tt = 1;
    end
    
  else
    % Interactive plotting: A gif is generated
    % Generate Gif: initial state
    drawnow
    frame = getframe(f);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if tt == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    tt = tt + 1;
    if tt > length(nameFolds)
      loopTime = 0;
    end
    
  end
  
end

end

