function [ data , tNow , GFcase ] = integrateLevelSet( GFcase , varargin )
%INTEGRATELEVELSET Integrates G-equation temporally

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Change dir to case root
change2caseRootDir( GFcase )

% Set user  defined parameters from vararin if specified
if nargin > 1
  % Each parameter has an unique name. Therefor varargin can just be
  % handles and matched with solver, p and vel
  % ATTENTION: varargin is not parsed for consistency or validity!
  [ GFcase.p ] = setUserParameters( GFcase.p , varargin );
  [ GFcase.solver ] = setUserParameters( GFcase.solver , varargin );
  [ GFcase.vel ] = setUserParameters( GFcase.vel , varargin );
  % Set indicator if initial logfile has been written to 0 so that changes
  % are documented!
  GFcase.p.logFileInit = 0;
end

% Prepare for output (create folders and initialize log file)
[ GFcase ] = prepare4output( GFcase , 'integrateLevelSet');

% Unwrap parameters
data = GFcase.data;
p = GFcase.p;
solver = GFcase.solver;
myGrid = GFcase.grid;
vel = GFcase.vel;

% If no local CFd data exists and there is a folder with CFD data in the
% root directory then copy data from root to local folder. Otherwise, if a
% local folder exists take this data
if strcmpi( vel.velModel , 'CFD_FOAM')
  if ~exist([p.run_output_folder,filesep,'FOAM_',vel.caseName],'dir') ...
        && exist([p.caseRootDir,filesep,'FOAM_',vel.caseName],'dir')
    % copy directory for velocity data from root to case dir
    sourceLocation = [p.caseRootDir,filesep,'FOAM_',vel.caseName];
    vel.outputDir = [p.run_output_folder,filesep,'FOAM_',vel.caseName];
    mkdir(vel.outputDir);
    copyfile([sourceLocation,filesep,'*'],vel.outputDir);
  elseif exist([p.run_output_folder,filesep,'FOAM_',vel.caseName],'dir')
    % take case velocity data if it exists
    vel.outputDir = [p.run_output_folder,filesep,'FOAM_',vel.caseName];
  else
    error('CFD velocity data not found!');
  end
  
end

% Update log file
logging( GFcase , 1 , 'Starting integration of G-equation...' )

% Set up solver (also opens files for output)
[ solverSetup , fileHandle ] = setSolverParameters( solver , p , myGrid , vel , data );


% Initialize Display
if solver.doPlot
  f = figure('visible','on','color','w');
  set(gcf,'Position',solver.figurePositionProps)
  
  % Pcolor or quiver plot?
  if solver.plotVelPcolor || solver.plotVelQuiver
    visualizeVelocityField( solver , solverSetup , myGrid , data );
  end
  
  % Plot level set (always)
  h = visualizeLevelSet(myGrid, data, solver.displayType, 0);
  set(h,'Color','r','LineWidth',4)
  
  % Plot grid?
  if solver.plotGrid
    visualizeGrid(myGrid);
  end
    
  % Export figure
  if strcmpi(solver.exportPlot,'gif')
    % Get filename for gif file
    filename = [p.run_output_folder,filesep,'XiFlame_',p.flameType,'.gif'];
    ii = 1;
    while exist(filename,'file')
      filename = [p.run_output_folder,filesep,'XiFlame_',p.flameType,num2str(ii),'.gif'];
      ii = ii + 1;
    end
    % Generate Gif: initial state
    drawnow
    frame = getframe(f);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    
  elseif strcmpi(solver.exportPlot,'png')
    % Get starting number for files
    filename = [p.run_output_folder,filesep,'out_',sprintf('%.4i',0)];
    indShift = 1;
    while exist(filename,'file')
      filename = [p.run_output_folder,filesep,'out_',sprintf('%.4i',indShift)];
      indShift = indShift + 1;
    end
    % Save snapshots
    print(f,[p.run_output_folder,filesep,'out_',sprintf('%.4i',indShift)],'-dpng','-r100','-opengl');
  end
  
end


% Loop until tMax (subject to a little roundoff).
tNow = solver.t0;
startTime = cputime;
myCount = 1;
while (solver.tMax - tNow > solver.small * solver.tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(solver.tMax, tNow + solverSetup.tOut) ];
  
  % Take a timestep.
  [ t , y , solverSetup.schemeData ] = feval(solverSetup.integratorFunc, solverSetup.schemeFunc, tSpan, y0,...
    solverSetup.integratorOptions, solverSetup.schemeData);
  tNow = t(end);
  
  % Get back the correctly shaped data array
  data = reshape(y, myGrid.shape);
  
  
  % Plotting
  if solver.doPlot
    % Get correct figure, and remember its current view.
    figure(f);
       
    % Delete last visualization if necessary.
    if(solver.deleteLastPlot)
      delete(h);
    end   
    
    % Pcolor or quiver plot?
    if solver.plotVelPcolor || solver.plotVelQuiver
      cla( get(f,'CurrentAxes') )
      visualizeVelocityField( solver , solverSetup , myGrid , data );
    end
     
    % Create new visualization.
    h = visualizeLevelSet(myGrid, data, solver.displayType, 0);
    set(h,'Color','r','LineWidth',4)
    
    % Plot grid?
    if solver.plotGrid
      visualizeGrid(myGrid);
    end
    
    % Export figure
    if strcmpi(solver.exportPlot,'gif')
      % Generate Gif
      drawnow
      frame = getframe(f);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,filename,'gif','WriteMode','append');
    elseif strcmpi(solver.exportPlot,'png')
      % Save snapshots
      print(f,[p.run_output_folder,filesep,'out_',sprintf('%.4i',indShift+myCount)],'-dpng','-r100','-opengl');
    end
    
  end
  
  % Writing of data (write always last time step out)
  if solver.doWriteData || ~(solver.tMax - tNow > solver.small * solver.tMax)
    % Wrap data and write 2 disc (only solver and vel could have been changed)
    GFcase.solver = solver;
    GFcase.vel = vel;
    GFcase.solverSetup = solverSetup;
    save([p.run_output_folder,'/',num2str(tNow,8),'.mat'],'data','GFcase')    
  end
  
  if ~solver.doWriteData && ~solver.doPlot
    warning('Neither ploting nor writing of data switched on. Consider to reduce <outSteps> to zero!')
  end
  
  myCount = myCount + 1;
end

endTime = cputime;
runtime = sprintf('Total execution time %g seconds\n', endTime - startTime);
logging( GFcase , 1 , 'Integration finished' )
logging( GFcase , 1 ,runtime )

% Close files
for ii = 1:length(fileHandle)
  fclose(fileHandle{ii});
end

% Final outputs
if strcmpi(solver.exportPlot,'png')
  disp('Convert series of png to gif animation:')
  disp('convert -layers optimize-plus -colors 128 -resize 70% -delay 20 -loop 0 out_*.png animation.gif')
  disp('Compress gif:')
  disp('gifsicle-static --lossy=100 -o animation_compressed.gif animation.gif')
end

% Wrap parameters (only solver and vel could have been changed)
GFcase.solver = solver;
GFcase.vel = vel;
GFcase.solverSetup = solverSetup;

end

