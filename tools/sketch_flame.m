function [ myF ] = sketch_flame( p , varargin )
%sketch_flame Creates sketch of a specific geometry (use together with flameResponsesAnalytic() )

%% Parse varargin
% Line style
ind = find(strcmpi(varargin,'lstyle'),1);
if ~isempty(ind)
  % User
  lstyle = varargin{ind+1};
else
  % Default
  lstyle = {'-','-.',':','--','-o','-.o',':o','--o','-x'};
end

% Where should ticz file be exported to? Specify *.tex file!
ind = find(strcmpi(varargin,'ticz2'),1);
if ~isempty(ind)
  % Use user specified limits
  ticz2 = varargin{ind+1};
  doEx2Ticz = 1;
else
  % use default 
  ticz2 = [];
  doEx2Ticz = 0;
end

% Positioning text
ind = find(strcmpi(varargin,'texPos'),1);
if ~isempty(ind)
  % User
  texPos = varargin{ind+1};
  texPosMode = 'manual';
else
  % Automatic
  texPos = [];
  texPosMode = 'auto';
end

% Specify figure to plot in
ind = find(strcmpi(varargin,'fig'),1);
if ~isempty(ind)
  % User specified figure
  myF = varargin{ind+1};
else
  % new figure
  myF = [];
end

% Line width flame front
ind = find(strcmpi(varargin,'linWidthF'),1);
if ~isempty(ind)
  % User specified lineWidth
  linWidthF = varargin{ind+1};
else
  % Default
  linWidthF = 0.7;
end

% Plot full of half geometry?
ind = find(strcmpi(varargin,'half'),1);
if ~isempty(ind)
  % Plot half geometry
  doFull = 0;
else
  % Plot full geometry
  doFull = 1;
end

% Add text?
ind = find(strcmpi(varargin,'noText'),1);
if ~isempty(ind)
  % No
  noText = 0;
else
  % Yes
  noText = 1;
end

% Plot in subplot?
ind = find(strcmpi(varargin,'subplot'),1);
if ~isempty(ind)
  % Yes
  doSubplot = 1;
else
  % No
  doSubplot = 0;
end

% factor x2 direction
ind = find(strcmpi(varargin,'facX2'),1);
if ~isempty(ind)
  % User specified
  facX2 = varargin{ind+1};
else
  % Default
  facX2 = 1.2;
end

% factor -x1 direction
ind = find(strcmpi(varargin,'facmX1'),1);
if ~isempty(ind)
  % User specified
  facmX1 = varargin{ind+1};
else
  % Default
  facmX1 = 0.05;
end

% Get color order
myColorOrd = get(0, 'DefaultAxesColorOrder');

% Check if input is cell array object,and if not convert it to one
if ~iscell(p)
  tmp = p; p = cell(1,1); p{1} = tmp;
end

%% Plot
if isempty(myF)
  myF = figure('Position',[0 50 800 550],'Color','w');hold on;
end

% Plot flames
x1_maxPlot = 0;
x1_minPlot = 0;
for ii=1:length(p)
  plot([0 p{ii}.H_flame],[p{ii}.R_flame,0],'LineStyle',lstyle{ii},'LineWidth',linWidthF,'Color', myColorOrd(ii,:))
  if doFull
    plot([0 p{ii}.H_flame],[-p{ii}.R_flame,0],'LineStyle',lstyle{ii},'LineWidth',linWidthF,'Color', myColorOrd(ii,:))
  end
  x1_maxPlot = max(p{ii}.H_flame,x1_maxPlot);
  x1_minPlot = min(p{ii}.H_flame,x1_minPlot);
  % Add text
  if noText
    if strcmpi(texPosMode,'auto')
      text( 0.8*p{ii}.H_flame-0.1 , 0.2*p{ii}.R_flame , ['$\alpha=',num2str(p{ii}.alphaDegree),'^\circ$'])
    else
      text(texPos(ii,1),texPos(ii,2),['$\alpha=',num2str(p{ii}.alphaDegree),'^\circ$'])
    end
  end
end

delta_x1 = x1_maxPlot - x1_minPlot;
x1_maxPlot = x1_maxPlot+0.1*delta_x1;
x1_minPlot = x1_minPlot-facmX1*delta_x1;

x2maxPlot = p{1}.R_flame * facX2;
x2minPlot = -p{1}.R_flame * facX2;

% Plot center line
line([x1_minPlot x1_maxPlot],[0,0],'LineStyle','-.','Color','k','LineWidth',1)
xlim([x1_minPlot x1_maxPlot])

if doFull
  ylim([x2minPlot,x2maxPlot])
else
  ylim([0,x2maxPlot])
end

% Plot confinement
currYlims = get(gca,'yLim');
currXlims = get(gca,'xLim');
line([currXlims(1) 0],[p{1}.R_i,p{1}.R_i],'Color','k','LineWidth',2)
if doFull
  line([currXlims(1) 0],[-p{1}.R_i,-p{1}.R_i],'Color','k','LineWidth',2)
  line([0 0],[-p{1}.R_i,currYlims(1)],'Color','k','LineWidth',2)
end
line([0 0],[p{1}.R_i,currYlims(2)],'Color','k','LineWidth',2)

% extend center line
line([currXlims(1) 0],[0,0],'LineStyle','-.','Color','k','LineWidth',1)

% axis off
if ~doSubplot
  set(gca,'Position',[0 0 1 1])
  axis off
  axis equal
end



if doEx2Ticz
  axoptions={'legend style={font=\footnotesize},',...
  'scaled y ticks = true',...
  'y tick label style={/pgf/number format/1000 sep={}, fixed, fixed zerofill,precision=1}'};
  cleanfigure('handle',myF,'minimumPointsDistance',5e-2)
  matlab2tikz(ticz2,'figurehandle',myF,...
    'Height', '\figureheight', 'Width', '\figurewidth', 'ShowInfo', false,'automaticLabels',true,...
    'maxChunkLength',4000,'parseStrings',false,'interpretTickLabelsAsTex',true,...
    'extraaxisoptions',axoptions);
end


end

