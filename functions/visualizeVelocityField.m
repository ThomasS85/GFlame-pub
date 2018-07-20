function h = visualizeVelocityField( solver , solverSetup , myGrid , data , varargin )
%VISUALIZEVELOCITYFIELD Does a pcolor and quiver plot of the velocity field
%   
% Inputs: 
%         solver      - GFLAME solver object
%         solverSetup - GFLAME solver setup struct containing
%             - schemeData - struct with innerFunc and innerData of termConvection
%           
%         grid        - GFLAME grid struct
%
% Only for 2D data!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.de).        //
% // Created, 11.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% Parse varargin
% Limits for velocity data
ind = find(strcmpi(varargin,'lim'),1);
if ~isempty(ind)
  % Use user specified limits
  limitsVel = varargin{ind+1};
else
  % use default limits from solver
  limitsVel = solver.pcolorLim;
end

% What velocity component to plot (pcolr)?
ind = find(strcmpi(varargin,'vComp'),1);
if ~isempty(ind)
  % Use user specified limits
  vComp = varargin{ind+1};
else
  % Plot absolute values
  vComp = 'abs';
end

% pcolor plot? 
ind = find(strcmpi(varargin,'doPcolor'),1);
if ~isempty(ind)
  % User input
  doPcolor = varargin{ind+1};
else
  % Default
  doPcolor = solver.plotVelPcolor;
end

% Quiver plot? 
ind = find(strcmpi(varargin,'doQuiver'),1);
if ~isempty(ind)
  % User input
  doQuiver = varargin{ind+1};
else
  % Default from solver
  doQuiver = solver.plotVelQuiver;
end

% Should burner be plotted?
ind = find(strcmpi(varargin,'doPlotBurner'),1);
if ~isempty(ind)
  % User input
  doPlotBurner = varargin{ind+1};
else
  % Default from solver
  doPlotBurner = solver.doPlotBurner;
end

% How many ticks should be used?
ind = find(strcmpi(varargin,'pcolorNTicks'),1);
if ~isempty(ind)
  % User input
  pcolorNTicks = varargin{ind+1};
else
  % Default from solver
  pcolorNTicks = solver.pcolorNTicks;
end

% Options for quiver plot: vector with [ dPlotX1 , dPlotX2 , scale ]
ind = find(strcmpi(varargin,'quiverOpt'),1);
if ~isempty(ind)
  % No arguments
  quiverOptions = varargin{ind+1};
  
else
  % options for quiver plot from solver
  quiverOptions = solver.quiverOptions;  
end

% What time is it (for title string)?
ind = find(strcmpi(varargin,'tNow'),1);
if ~isempty(ind)
  % Use user specified time
  tNow = varargin{ind+1};
else
  % Do not plot time
  tNow = 'none';
end


%% Extract data from solver object
% Index of convective schemeData
tmp = cellfun(@(x) strcmp(func2str(x),'termConvection') , solverSetup.schemeData.innerFunc );
indSD = find(tmp,1,'first');

if isempty(indSD)
  error('No convective term specified or name of convective term unknown!')
end

% Velocity model
% velMod = solverSetup.schemeData.innerData{indSD}.velModel;

% velocity field
velocity = solverSetup.schemeData.innerData{indSD}.velocityField;



%% Prepare velocity data
% Expand each field
for ii = 1:length(velocity)
  if length(velocity{ii})==1
    % Uniform velocity field -> Expand to whole field
    velocity{ii} = ones(myGrid.shape) * velocity{ii};
  end
end

% calculate absolute value
if strcmp(vComp,'abs')
  velocity_abs = sqrt( velocity{1}.^2 + velocity{2}.^2 );
end


%% Plot options
% Hold on so that previous plots don't get deleted
hold on;
% Which velocity component to plot?
switch vComp
  case 'abs'
    velocity_plot = velocity_abs;
    myLableCol = '$|u|$~[m/s]';
  case 'v1'
    velocity_plot = velocity{1};
    myLableCol = '$u_1$~[m/s]';
  case 'v2'
    velocity_plot = velocity{2};
    myLableCol = '$u_2$~[m/s]';
end

% Compute some derived quantities
% Define ticks for pcolor
myTicks = unique(round( linspace(limitsVel(1),limitsVel(2),pcolorNTicks)*10 ) / 10);
% Define unit for lengths
myOMag = floor( log10(max(myGrid.vs{1})));
switch myOMag
  case {1,0,-1}
    % m
    myXlabel = '$x_1$~[m]';myYlabel = '$x_2$~[m]';
    myScaleFac = 1;
  case {-2,-3,-4}
    % mm
    myXlabel = '$x_1$~[mm]';myYlabel = '$x_2$~[mm]';
    myScaleFac = 1e3;
  otherwise
    % m
    myXlabel = '$x_1$~[m]';myYlabel = '$x_2$~[m]';
    myScaleFac = 1;
end
% Define ticks x1 and x2 axis
nTicks = 3;
myTicksX =  unique( round( linspace(0,max(myGrid.vs{1}),nTicks)*100*myScaleFac ) / 100/myScaleFac );
myTicksY =  unique( round( linspace(min(myGrid.vs{2}),max(myGrid.vs{2}),nTicks)*100*myScaleFac ) / 100/myScaleFac );
% Define labels x1 and x2 axis
myTicksXLab = cellstr(num2str(myTicksX'*myScaleFac))';
myTicksYLab = cellstr(num2str(myTicksY'*myScaleFac))';
% Get actual time + decide what unit to plot it
if ~strcmpi(tNow,'none')
  t_now = solverSetup.schemeData.innerData{indSD}.t_vec(1);
  myOMagT = floor( log10(t_now));
  switch myOMagT
    case {1,0,-1}
      % s
      myTUnit = '[s]';
      myTScaleFac = 1;
    case {-2,-3,-4}
      % ms
      myTUnit = '[ms]';
      myTScaleFac = 1e3;
    otherwise
      % s
      myTUnit = '[s]';
      myTScaleFac = 1;
  end
end

if solver.lablesOff
  myXlabel = ''; myYlabel = '';
  myTicksX = [min(myTicksX) max(myTicksX)]; myTicksY = [min(myTicksY) max(myTicksY)];
  myTicksXLab = '';myTicksYLab = '';
  myLableCol = '';
  myTicks = [];
end

%% Do plot pcolor
if doPcolor
  velocity_plot(data>=0)=nan; % only plot inside fresh region
  h = pcolor(myGrid.xs{1},myGrid.xs{2},velocity_plot);   % Axial Velocity Component
  set(h,'EdgeColor','none');
  shading(gca,'interp')
  
  % set(gca, 'CLim', [min(velocity_plot(:)), max(velocity_plot(:))]);
  set(gca, 'CLim', limitsVel);
  set(gca,'TickLabelInterpreter','Latex');
  
  cmap = getCustomColormap(1024,'BuRd');
  colormap(cmap);
  if ~solver.lablesOff
    cb = colorbar('Ticks',myTicks);
    zlab = get(cb,'ylabel');
    set(zlab,'String',myLableCol);set(cb,'TickLabelInterpreter','Latex');
    set( zlab ,'Interpreter', 'Latex' );
  end
  
  % scale axis
  xlabel(myXlabel,'Interpreter', 'Latex');
  ylabel(myYlabel,'Interpreter', 'Latex');
  set(gca,'XTick',myTicksX);set(gca,'YTick',myTicksY);
  set(gca,'XTickLabel',myTicksXLab);set(gca,'YTickLabel',myTicksYLab);
  set(gca,'FontSize',32);
  
  if ~strcmpi(tNow,'none')
    title([ 't = $', num2str(t_now*myTScaleFac,3),'$~',myTUnit ],'FontSize',32,'Interpreter', 'Latex' );
  end
  
  axis off
  grid off
  
end


%% Do plot quiver
if doQuiver
  dPlotX1 = quiverOptions(1);
  dPlotX2 = quiverOptions(2);
  scale = quiverOptions(3);
  h = quiver(myGrid.xs{1}(1:dPlotX1:end,1:dPlotX2:end),myGrid.xs{2}(1:dPlotX1:end,1:dPlotX2:end),...
      velocity{1}(1:dPlotX1:end,1:dPlotX2:end),velocity{2}(1:dPlotX1:end,1:dPlotX2:end) , scale ,'color','k' );
end


%% Plot burner
if doPlotBurner>0
  % Get general settings
  p = solverSetup.schemeData.innerData{indSD}.p;
  % Should confinement be plotted?
  doConfinement = (doPlotBurner==2 || strcmpi(p.geom,'V') || strcmpi(p.geom,'M'));
  % Increase xlim and y lim
  add2XLim = myTicksX(end)*0.1;
  xlim([myTicksX(1)-add2XLim myTicksX(end)])
  
  if doConfinement
    ylim([myTicksY(1) p.R_a])
  else
    add2YLim = myTicksY(end)*0.05;
    ylim([myTicksY(1) myTicksY(end)+add2YLim])
  end
  % Plot central axis
%   line([myTicksX(1)-add2XLim myTicksX(end)],[0 0],'Color','k','lineWidth',2,'LineStyle','-.')
  % Plot burner
  line([myTicksX(1)-add2XLim myTicksX(1)],[p.R_i p.R_i],'Color','k','lineWidth',6,'LineStyle','-')
  if doConfinement
    line([myTicksX(1) myTicksX(1)],[p.R_i p.R_a],'Color','k','lineWidth',6,'LineStyle','-')
    line([myTicksX(1) myTicksX(end)],[p.R_a p.R_a],'Color','k','lineWidth',6,'LineStyle','-')
    if strcmp(p.domainType,'full')
      line([myTicksX(1) myTicksX(1)],-[p.R_i p.R_a],'Color','k','lineWidth',6,'LineStyle','-')
      line([myTicksX(1) myTicksX(end)],-[p.R_a p.R_a],'Color','k','lineWidth',6,'LineStyle','-')
      line([myTicksX(1)-add2XLim myTicksX(1)],-[p.R_i p.R_i],'Color','k','lineWidth',6,'LineStyle','-')
    end
  else
    line([myTicksX(1) myTicksX(1)],[p.R_i myTicksY(end)+add2YLim],'Color','k','lineWidth',6,'LineStyle','-')
    if strcmp(p.domainType,'full')
      line([myTicksX(1) myTicksX(1)],-[p.R_i myTicksY(end)+add2YLim],'Color','k','lineWidth',6,'LineStyle','-')
      line([myTicksX(1)-add2XLim myTicksX(1)],-[p.R_i p.R_i],'Color','k','lineWidth',6,'LineStyle','-')
    end
  end
  % Plot central rod
  if strcmpi(p.geom,'V') || strcmpi(p.geom,'M')
    line([myTicksX(1)-add2XLim myTicksX(1)],[p.R_rod p.R_rod],'Color','k','lineWidth',6,'LineStyle','-')
    line([myTicksX(1) myTicksX(1)],[0 p.R_rod],'Color','k','lineWidth',6,'LineStyle','-')
    if strcmp(p.domainType,'full')
      line([myTicksX(1)-add2XLim myTicksX(1)],-[p.R_rod p.R_rod],'Color','k','lineWidth',6,'LineStyle','-')
      line([myTicksX(1) myTicksX(1)],-[0 p.R_rod],'Color','k','lineWidth',6,'LineStyle','-')
    end
  end
  
end


end

