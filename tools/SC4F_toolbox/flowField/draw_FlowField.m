function [ myF , plts ] = draw_FlowField( X_up_mesh , X_down_mesh , vel , x1L_flame , p , varargin )
%DRAW_FLOWFIELD Plots flow field including flame etc in L1 coordiantes
%
% Inputs:
%     X_up_mesh   : All points at which velocity data is available UPSREAM of the area jump/ flame (flat flame)
%                     Cell array { X1_up , X2_up } where [X1_up,X2_up] = meshgrid(x1_up,x2_up);
%     X_down_mesh : All points at which velocity data is available DOWNSTREAM of the area jump/ flame (flat flame)
%                     Cell array { X1_down , X2_down } where [X1_down,X2_down] = meshgrid(x1_down,x2_down);
%     vel         : velocity struct
%     x1L_flame   : Flame position as complex numbers (L1 coordinates)
%     p           : Struct with flame settings as returned from setUpPredefinedFlame()
%
% by Thomas Steinbacher Oct 2018

%% Parse varargin
% Specify figure to plot in
ind = find(strcmpi(varargin,'fig'),1);
if ~isempty(ind)
  % User specified figure
  myF = varargin{ind+1};
else
  % new figure
  myF = [];
end

% Specify mean flame position
ind = find(strcmpi(varargin,'meanF'),1);
if ~isempty(ind)
  % User specified mean flame position
  x1L_flame_mean = varargin{ind+1};
else
  % No mean flame position
  x1L_flame_mean = [];
end

% Specify Xi grid in orer to speed up velocity field computation (no mapping reqired)
ind = find(strcmpi(varargin,'Xi'),1);
if ~isempty(ind)
  % User specified scaling
  Xi_all = varargin{ind+1};
  Xi_up = Xi_all{1};
  Xi_down = Xi_all{2};
else
  % Unity scaling
  Xi_all = [];
end

% Plot Source locations realted to  gas expansion
ind = find(strcmpi(varargin,'sourceL'),1);
if ~isempty(ind)
  % Get indices of gas expansion sources
  ind_sources = varargin{ind+1};
else
  % Do not plot
  ind_sources = [];
end

% What to plot?   
%  abs     -> Absolute velocity    (Default)
%  u1      -> x1 component flow
%  u2      -> x2 component flow
ind = find(strcmpi(varargin,'plotField'),1);
if ~isempty(ind)
  % Get user choice
  plotField = varargin{ind+1};
else
  % Do not plot
  plotField = 'abs';
end

% Add quiver plot?
ind = find(strcmpi(varargin,'doQuiver'),1);
if ~isempty(ind)
  % yes
  doQuiver = 1;
  % Scaling factor for quiver plot
  myFac = varargin{ind+1};
else
  % N0
  doQuiver = 0;
end

% Add stream lines?
ind = find(strcmpi(varargin,'doStreamLines'),1);
if ~isempty(ind)
  % yes
  doStreamLines = 1;
  Ns = varargin{ind+1}; % number of stream lines
else
  % N0
  doStreamLines = 0;
end

% Color limits
ind = find(strcmpi(varargin,'cLim'),1);
if ~isempty(ind)
  % yes
  myClim = varargin{ind+1};
else
  % N0
  myClim = [];
end


% Plot Kutta panel?
ind = find(strcmpi(varargin,'Kutta'),1);
if ~isempty(ind)
  % yes
  doPlKutta = 1;
else
  % N0
  doPlKutta = 0;
end

%% Pre-process
% Convert mesh to vector in order to plot velocity (quiver)
X_up = X_up_mesh{1}(:) + 1i*X_up_mesh{2}(:);
X_down = X_down_mesh{1}(:) + 1i*X_down_mesh{2}(:);
% Combine up and down
X_all = [X_up;X_down];
  
% Get limits
myX1Lims = [ min(real(X_all)) , max(real(X_all))];
myX2Lims = [ min(imag(X_all)) , max(imag(X_all))];
% add puffer
myX1Lims = myX1Lims + [ -myX1Lims(1)*0.05 , myX1Lims(2)*0.05 ];
myX2Lims = myX2Lims + [ -myX2Lims(1)*0.03 , myX2Lims(2)*0.03 ];


% Set up Kutta condition
if strcmpi(p.CombType,'backwardFacingStep')
  % Enable Kutta condition
  doKutta = 'doKutta';
else
  % Disable Kutta condition
  doKutta = 'noKutta';
end

% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );
% Map mesh if no Xi is specified
if isempty(Xi_all)
  Xi_up =  myMap.xi_x_L1( X_up );
  Xi_down =  myMap.xi_x_L1( X_down );
end

% Evaluate velocity
[ u_up , vel_up ] = evalFlowField_physicalDomain( Xi_up , vel , myMap , doKutta , 1 ,...
    'xi' , 'L1' , 'noPanel2Point' );
[ u_down , vel_down ] = evalFlowField_physicalDomain( Xi_down , vel , myMap , doKutta , 1 ,...
    'xi' , 'L1' , 'noPanel2Point' );
u_all = [ u_up ; u_down ];  
  
% Evaluate absolute velocity
if strcmpi(plotField,'abs')
  absVel_up = sqrt(real(u_up(:)).^2 + imag(u_up(:)).^2);
  absVel_down = sqrt(real(u_down(:)).^2 + imag(u_down(:)).^2);
  absVel_lim = [ min(min(absVel_up),min(absVel_down)) , max(max(absVel_up),max(absVel_down)) ];
end

% Compute stream lines
u_up_mat = reshape(u_up,size(X_up_mesh{1}));
u_down_mat = reshape(u_down,size(X_down_mesh{1}));

if doStreamLines
  % Upper flow field
  ds_tmp = p.R_i/Ns;
  x1_s = zeros(Ns,1)+myX1Lims(1)+1e-4;x2_s = linspace(ds_tmp/2,p.R_i-ds_tmp/2,Ns).';
  myStrmlns_up = stream2(X_up_mesh{1},X_up_mesh{2},real(u_up_mat),imag(u_up_mat),x1_s,x2_s);
  
  % Downstream flow field
  % Starting point of downstream stream lines are ending points of upstream one
  x1_s = zeros(length(myStrmlns_up),1); x2_s = zeros(length(myStrmlns_up),1);
  for ii=1:length(myStrmlns_up)
    x1_s(ii) = min(real(X_down))+1e-4;
    x2_s(ii) = myStrmlns_up{ii}(end,2);
  end 
  myStrmlns_down = stream2(X_down_mesh{1},X_down_mesh{2},real(u_down_mat),imag(u_down_mat),x1_s,x2_s);
end


%% Prepare plot
if isempty(myF)
  myF = figure('Position',[0 50 800 550],'Color','w');hold on;
  figIni = 1;
else
  figIni = 0;
end

% axis off
set(gca,'Position',[0 0 1 1])
axis off
axis equal


%% Plot flow field
switch plotField
  case 'abs'
    % pcolor plot (absolute velocity)
    plts{1} = pcolor(X_up_mesh{1},X_up_mesh{2},reshape(absVel_up,size(X_up_mesh{1})));
    plts{2} = pcolor(X_down_mesh{1},X_down_mesh{2},reshape(absVel_down,size(X_down_mesh{1})));
    if isempty(myClim)
      myClim = [p.s_l_u*0.5 p.E*p.s_l_u*1.8];
    end
  case 'u2'
    plts{1} = pcolor(X_up_mesh{1},X_up_mesh{2},imag(u_up_mat));
    plts{2} = pcolor(X_down_mesh{1},X_down_mesh{2},imag(u_down_mat));
    if isempty(myClim)
      myClim = [-1 1]*0.1;
    end
  case 'u1'
    plts{1} = pcolor(X_up_mesh{1},X_up_mesh{2},real(u_up_mat));
    plts{2} = pcolor(X_down_mesh{1},X_down_mesh{2},real(u_down_mat));
    if isempty(myClim)
      myClim = [p.s_l_u*0.5 p.E*p.s_l_u*1.8];
    end
end
nPlts = length(plts);

% quiver plot (arrows)
if doQuiver
  nth = 2;
  plts{nPlts+1} = quiver(real(X_all(1:nth:end)),imag(X_all(1:nth:end)),...
    myFac*real(u_all(1:nth:end)),myFac*imag(u_all(1:nth:end)),0,'Color','k');
end
nPlts = length(plts);

% plot stream lines
if doStreamLines
  for ii=1:length(myStrmlns_up)
    mySL_tmp = [ myStrmlns_up{ii} ; myStrmlns_down{ii} ];
    plts{nPlts+ii} = plot(mySL_tmp(:,1),mySL_tmp(:,2),'k-','LineWidth',1);
  end
end
nPlts = length(plts);

% plot source locations
if ~isempty(ind_sources)
  sourceLoc = L2_to_L1(vel.sourceDat.x(ind_sources),p);
  plts{nPlts+1} = plot( real(sourceLoc) , imag(sourceLoc) , 'ro','MarkerFaceColor','r','markers',3 );
  nPlts = length(plts);
end

% Plot panel
if isfield(vel,'panel')
  for pp=2:length(vel.panel(1,:))
    if strcmpi(vel.panel{1,pp},'source')
      myLStyle = 'ro-';
      myCol = 'r';
    else
      myLStyle = 'go-';
      myCol = 'g';
    end
    plts{nPlts+1} = plot(L2_to_L1(vel.panel{3,pp},p),myLStyle,'MarkerFaceColor',myCol,'markers',3);
    nPlts = length(plts);
  end
end

% Plot Kutta Panel
if doPlKutta && isfield(vel_up.Kutta,'H')
  myBeta = vel_up.Kutta.beta;
  xi_pSL = linspace(1,1+vel_up.Kutta.H*cos(myBeta),30) + 1i*linspace(0,vel_up.Kutta.H*sin(myBeta),30);
  x_pSL = L2_to_L1(myMap.x_xi(xi_pSL),p);
  plot(x_pSL,'k-','LineWidth',2)
end

%% Plot flame
% Plot instantaneouse flame
plts{nPlts+1} = plot(real(x1L_flame),imag(x1L_flame),'LineStyle','-','LineWidth',2,'Color', 'g');

% Plot mean flame position if desired
if ~isempty(x1L_flame_mean)
  plts{nPlts+2} = plot(real(x1L_flame_mean),imag(x1L_flame_mean),'LineStyle','--','LineWidth',2,'Color', 'b');
end
nPlts = length(plts);


%% Plot combustor
% Plot center line
line(myX1Lims,[0,0],'LineStyle','-.','Color','k','LineWidth',2)

% Plot confinement
line([myX1Lims(1) 0],[p.R_i,p.R_i],'Color','k','LineWidth',4)
line([0 0],[p.R_i,p.R_a],'Color','k','LineWidth',4)
line([0 myX1Lims(2)],[1 1]*p.R_a,'Color','k','LineWidth',4)



%% Format
% Limits
xlim(myX1Lims)
ylim(myX2Lims)

% Color map etc
set(plts{1},'EdgeColor','none');
set(plts{2},'EdgeColor','none');
% shading(gca,'interp')
set(gca,'TickLabelInterpreter','Latex');

% color map
cmap = getCustomColormap(1024,'BuRd');
colormap(cmap);

% Color bar
% myTicks = floor(absVel_lim):0.5:absVel_lim;
myTicks = floor(myClim(1)*2)/2 : 1: ceil(myClim(2));
cb = colorbar('Ticks',myTicks);
set(gca, 'CLim', myClim);
zlab = get(cb,'ylabel');
set(cb,'TickLabelInterpreter','Latex');
set( zlab ,'Interpreter', 'Latex' );


% adjust figure
if figIni
  myRatio = (myX1Lims(2)-myX1Lims(1)) / (myX2Lims(2)-myX2Lims(1));
  if strcmpi(p.CombType,'duct')
    myHeight = 300;
  else
    myHeight = 700;
  end
  set(myF,'Position',[0 50 myHeight*myRatio myHeight*1.1])
end


end



