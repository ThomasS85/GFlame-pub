function [ myF ] = visualize_velField( p , varargin )
%VISUALIZE_VELFIELD Visualizes geometry and specified velocity field. default is irrotational field due to
%  source which imposes a velocity of 1 at x1->-infinity. Solenoidal velocity can be chosen via varargin
%
% Inputs:
%   - p        :  Struct with flame settings as returned from setUpPredefinedFlame()
%
% Outputs:
%   - myF      : Function handle to figure
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 28.05.2015 as part of GFLAME 0.1          //
% // Last modified: 28.05.2015 by steinbacher           //
% ////////////////////////////////////////////////////////

%% Parse varargin
% Number of velocity field sample points upstream [nx1 nx2]
ind = find(strcmpi(varargin,'nPointsUP'),1);
if ~isempty(ind)
  % Use user specified value
  nPointsUP = varargin{ind+1};
else
  % use default value 
  nPointsUP = [1 6];
end

% Number of velocity field sample points downstream [nx1 nx2]
ind = find(strcmpi(varargin,'nPointsDOWN'),1);
if ~isempty(ind)
  % Use user specified value
  nPointsDOWN = varargin{ind+1};
else
  % use default value 
  nPointsDOWN = [20 13];
end

% Should, insteady of potential solution, field of a vortex be shown?
ind = find(strcmpi(varargin,'vort'),1);
if ~isempty(ind)
  % Yes visualize vortex field -> input of cell array { xs , Gamma }
  vel_param = varargin{ind+1};
  doVort = 1;
else
  % Default: Visualize potential velocity field with u=1 at x1->-infinity
  vel_param = {1};
  doVort = 0;
end

% Scaling factor for lengths of arrows
ind = find(strcmpi(varargin,'fac'),1);
if ~isempty(ind)
  % Use user specified value
  myFac = varargin{ind+1};
else
  % use default value 
  myFac = 1;
end


%% Plot flame and casing
% New figure
myF = figure('Position',[0 50 1000 800],'Color','w','Name','Physical Domain');hold on;
% Plot Flame and casings
tmp = cell(1,1);tmp{1}=p;
sketch_flame( tmp , 'noText' , 'half' , 'fig' , myF );hold on;


%% Plot velocity field
% Generate sample points upstream of area jump (L1)
x1_up = linspace(-0.1*p.H_flame,-0.05*p.H_flame,nPointsUP(1));
x2_up = linspace(0,p.R_i,nPointsUP(2));
[X1_up,X2_up] = meshgrid(x1_up,x2_up);
X_up = X1_up(:) + 1i*X2_up(:);
% Generate sample points downstream of area jump (L1)
x1_down = linspace(0,1.3*p.H_flame,nPointsDOWN(1));
x2_down = linspace(0,p.R_a,nPointsDOWN(2));
[X1_down,X2_down] = meshgrid(x1_down,x2_down);
X_down = X1_down(:) + 1i*X2_down(:);
% Evaluate velocity field at sample points
if doVort
  % Solenoidal velocity field
  u_up = evalVel_vortex( X_up , vel_param{1} , vel_param{2} , p , 'L1' );
  u_down = evalVel_vortex( X_down , vel_param{1} , vel_param{2} , p , 'L1' );
else
  % Irrotational velocity field
  u_up = evalVel_acoustic( X_up , vel_param{1} , p , 'L1' );
  u_down = evalVel_acoustic( X_down , vel_param{1} , p , 'L1' );
end
% Now plot
quiver(X1_up(:),X2_up(:),myFac*real(u_up(:)),myFac*imag(u_up(:)),0,'Color','k');
quiver(X1_down(:),X2_down(:),myFac*real(u_down(:)),myFac*imag(u_down(:)),0,'Color','k');


end

