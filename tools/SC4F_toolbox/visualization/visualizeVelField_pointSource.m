function [ myF , sp ] = visualizeVelField_pointSource( Xi_grid , X_grid , vel , myMap , varargin )
%VISUALIZEVELFIELD_POINTSOURCE Visualizes velocity field for a given distribution of vortexes and sources
%
%
%
%
% by Thomas Steinbacher


%% Parse varargin
% Plot vortex paths?
ind = find(strcmpi(varargin,'VortPath'),1);
if ~isempty(ind)
  % Use user specified path
  VortPath = varargin{ind+1};
else
  % Default: No
  if ~isempty(vel.vortDat)
    VortPath = zeros( length(vel.vortDat.G) , 1 );
  else
    VortPath = zeros( 1 , 1 );
  end
end

% Plot vortex paths in image domain? (Only if 'doImage' is chosen)
ind = find(strcmpi(varargin,'VortPathIm'),1);
if ~isempty(ind)
  % Use user specified path
  VortPathIm = varargin{ind+1};
else
  % Default: No
  if ~isempty(vel.vortDat)
    VortPathIm = zeros( length(vel.vortDat.G) , 1 );
  else
    VortPathIm = zeros( 1 , 1 );
  end
end

% Also visualize image domain?
ind = find(strcmpi(varargin,'doImage'),1);
if ~isempty(ind)
  % Yes
  doImage = 1;
  plXi = varargin{ind+1};
else
  % Default: No
  doImage = 0;
end

% In which figure should be plotted?
ind = find(strcmpi(varargin,'fig'),1);
if ~isempty(ind)
  % Use user specified figure
  myF = varargin{ind+1};
  clf(myF)
  hold on;
else
  % Default: Create new figure
  if doImage
    myF = figure('Position',[0 50 1600 800],'Color','w');
  else
    myF = figure('Position',[0 50 800 800],'Color','w');
  end
  hold on;
end

% Add user specified velocity to each point in physical domain?
ind = find(strcmpi(varargin,'addVel'),1);
if ~isempty(ind)
  % Yes
  addVel = varargin{ind+1};
else
  % Default: No
  addVel = [];
end

%% Prepare
% Color order vor vortex paths
if ~isempty(vel.vortDat)
  if length(vel.vortDat.G)<=7
    myColorOrd = get(0, 'DefaultAxesColorOrder');
  else
    myColorOrd = repmat( [1 1 1]*0.3 , length(vel.vortDat.G) , 1 );
  end
end

% Get real valued grid data
X1_grid = real(X_grid);
X2_grid = imag(X_grid);

% Unwrap velocity field data
vortDat = vel.vortDat;
sourceDat = vel.sourceDat;

% Get limits for plotting
x1_lims = [ min(X1_grid(:)) , max(X1_grid(:)) ];
x2_lims = [ min(X2_grid(:)) , max(X2_grid(:)) ];
maxAbsXVal = max( max( abs(x1_lims) ) , max( abs(x2_lims) ) );
deltaX1 = x1_lims(2) - x1_lims(1);
deltaX2 = x2_lims(2) - x2_lims(1);
Nx1 = deltaX1/10; Nx2 = deltaX2/10;

if doImage
  xi1_lims = plXi.xiLims{1};
  xi2_lims = plXi.xiLims{2};
  deltaXi1 = xi1_lims(2) - xi1_lims(1);
  deltaXi2 = xi2_lims(2) - xi2_lims(1);
  Nxi1 = deltaXi1/10; Nxi2 = deltaXi2/10;
end

% Should Routh's correction be applied?
if ~isempty(vortDat)
  if isfield(myMap,'RouthsCorr') && any( vortDat.r0==0 )
    plotRouthsCorr = 1;
  else
    plotRouthsCorr = 0;
  end
else
  % No vortices
  plotRouthsCorr = 0;
end

%% Plot
figure(myF)

if doImage
  sp{1}=subplot(1,2,1); hold on;
else
  sp = {};
end
title('Physical Domain','Interpreter','Latex')

% (A) Physical domain
% Plot arrows
[ u_res_grid ] = evalFlowField_PointSource( Xi_grid(:) , vel , 'myMapping' , myMap )...
  .* conj( myMap.dxi_dx( Xi_grid(:) ) );
% Filter out too high velocity components
elimVel = abs(u_res_grid)>0.2*maxAbsXVal;
u_res_grid(elimVel) = 0;
% Add velocity specified by user (if desired)
if ~isempty(addVel)
  u_res_grid = u_res_grid + addVel;
end
% Now plot
quiver(X1_grid(:),X2_grid(:),real(u_res_grid),imag(u_res_grid),0,'Color','b');
% plot vortex/ source centers
if ~isempty(vortDat); plot(vortDat.x,'og','MarkerSize',12,'LineWidth',3); end
if ~isempty(sourceDat); plot(sourceDat.x,'xc','MarkerSize',12,'LineWidth',3); end
% plot vortex path
if ~isempty(vortDat)
  for nn=1:length(vortDat.G)
    plot(VortPath(nn,:),'Color',myColorOrd(nn,:),'LineWidth',2);
    % Visualize r0
    if vel.vortDat.r0(nn)>0
      [ x1_circ , x2_circ ] = myCircle( real(vortDat.x(nn)) , imag(vortDat.x(nn)) , vortDat.r0(nn) , 50 );
      plot(x1_circ,x2_circ,'--','Color',myColorOrd(nn,:),'LineWidth',1)
    end
  end
end

% Plot Routh's correction
if plotRouthsCorr && ~isempty(vortDat)
  myRouthCorr = myMap.RouthsCorr( vortDat.G , vortDat.xi );
  quiver( real(vortDat.x(:)) , imag(vortDat.x(:)) , real(myRouthCorr) , imag(myRouthCorr) , 0 , 'Color' , 'r' ,...
    'LineWidth',2)
end

% Plot settings
xlim(x1_lims);ylim(x2_lims); grid on; box on;
xlabel('$x_1$','Interpreter','Latex');ylabel('$x_2$','Interpreter','Latex');
ax = gca;
ax.XTick = x1_lims(1):Nx1:x1_lims(2);
ax.YTick = x2_lims(1):Nx2:x2_lims(2);
ax.XTickLabel = {};
ax.YTickLabel = {};
set(gca,'FontSize',24)

% Plot Vertexes
if isfield(myMap,'vertLables') && isfield(myMap,'vertexes') && isfield(myMap,'vertRelPos')
  for ii=1:length(myMap.vertexes)
    if ~isinf(myMap.vertexes(ii))
      myVertLabPos = myMap.vertexes(ii) + myMap.vertRelPos(ii);
      text( real(myVertLabPos) , imag(myVertLabPos) , myMap.vertLables{ii} , 'FontSize' , 24 ,...
        'HorizontalAlignment','center','Interpreter','Latex' )
    end
  end
end


% (B) Image domain
if doImage
  sp{2}=subplot(1,2,2); hold on;
  title('Image Domain','Interpreter','Latex')
  
  % remove u_p for plotting in image domain
  vel.u_p = 0;
  
  % Shade wall region in grey
  patchPoints = [ xi1_lims(1) , xi1_lims(2) , xi1_lims(2) , xi1_lims(1) ; 0 , 0 , xi2_lims(1) , xi2_lims(1) ];
  a2 = patch(patchPoints(1,:) , patchPoints(2,:),'k');
  set(a2,'FaceAlpha',0.5,'FaceColor',[1 1 1]*0.7,'EdgeColor','k','EdgeAlpha',0)
  
  % Plot arrows
  [ u_res_grid_Im ] = evalFlowField_PointSource( plXi.Xi1_grid(:)+1i*plXi.Xi2_grid(:) , vel , 'myMapping' , myMap );
  % Filter out too high velocity components
  elimVel = abs(u_res_grid_Im)>0.1*maxAbsXVal/plXi.myscal;
  u_res_grid_Im(elimVel) = 0;
  % Now plot
  quiver(plXi.Xi1_grid(:),plXi.Xi2_grid(:),real(u_res_grid_Im)*plXi.myscal,imag(u_res_grid_Im)*plXi.myscal,0,'Color','b');
  % plot vortexes/ sources centers
  if ~isempty(vortDat); plot(vortDat.xi,'og','MarkerSize',12,'LineWidth',3); end
  if ~isempty(sourceDat); plot(sourceDat.xi,'xc','MarkerSize',12,'LineWidth',3); end
  % Plot mirror vortrexes and sources
  if ~strcmpi(vel.doMirror,'noMirror') && ~isempty(vortDat)
    plot(conj(vortDat.xi),'or','MarkerSize',12,'LineWidth',3)
  end
  if ~strcmpi(vel.doMirror,'noMirror') && ~isempty(sourceDat)
    plot(conj(sourceDat.xi),'xy','MarkerSize',12,'LineWidth',3)
  end
  % plot vortex path
  if ~isempty(vortDat)
    for nn=1:length(vortDat.G)
      plot(VortPathIm(nn,:),'Color',myColorOrd(nn,:),'LineWidth',2);
      if ~strcmpi(vel.doMirror,'noMirror')
        % Plot path mirror vortex
        plot(conj(VortPathIm(nn,:)),'--','Color',myColorOrd(nn,:),'LineWidth',2);
      end
      % Visualize r0
      if vel.vortDat.r0(nn)>0
        [ xi1_circ , xi2_circ ] = myCircle( real(vortDat.xi(nn)) , imag(vortDat.xi(nn)) ,...
          vel.vortDat.r0(nn)*abs( myMap.dxi_dx(vortDat.xi(nn)) ) , 50 );
        plot(xi1_circ,xi2_circ,'--','Color',myColorOrd(nn,:),'LineWidth',1)
      end
    end
  end
  % Plot wall
  line(xi1_lims,[0 0],'LineWidth',5,'Color','k');
    
  % Plot settings
  xlim(xi1_lims);ylim(xi2_lims); grid on; box on;
  xlabel('$\xi_1$','Interpreter','Latex');ylabel('$\xi_2$','Interpreter','Latex');
  ax = gca;
  ax.XTick = xi1_lims(1):Nxi1:xi1_lims(2);
  ax.YTick = xi2_lims(1):Nxi2:xi2_lims(2);
  ax.XTickLabel = {};
  ax.YTickLabel = {};
  set(gca,'FontSize',24)
  
  sp{1}.Position=[0.03,sp{1}.Position(2),sp{1}.Position(3)*1.3,sp{1}.Position(4)];
  sp{2}.Position=[0.55,sp{2}.Position(2),sp{2}.Position(3)*1.3,sp{2}.Position(4)];
  
  % Plot prevertexes
  if isfield(myMap,'prevertexes') && isfield(myMap,'vertLables')
    for ii=1:length(myMap.prevertexes)
      if ~isinf(myMap.prevertexes(ii))
        plot( [myMap.prevertexes(ii) myMap.prevertexes(ii)] , [ 0 deltaXi2*0.02 ] , '-k' , 'LineWidth' , 2 )
        text( myMap.prevertexes(ii) , -deltaXi2*0.04 , myMap.vertLables{ii} , 'FontSize' , 24 ,...
          'HorizontalAlignment','center','Interpreter','Latex' )
      end
    end
  end
  
end

end

