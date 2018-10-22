function [ Xi_all , X_up_mesh , X_down_mesh ] = make_grid4plot( p )
%MAKE_GRID4PLOT Summary of this function goes here
%   Detailed explanation goes here

% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );

% Reference x1-length
if ~strcmp(p.flameType,'flat')
  x1_refL = p.H_flame;
else
  x1_refL = p.R_i*1;
end

% Define grid
if strcmpi(p.CombType,'duct')
  % Number of sample points
  nPointsUP = [100 80];
  nPointsDOWN = [100 80];
  % Generate sample points upstream of area jump (L1)
  x1_up = linspace(-1.3*x1_refL,-0.000001*x1_refL,nPointsUP(1));
  x2_up = linspace(0,p.R_i,nPointsUP(2));
else
  % Number of sample points
  nPointsUP = [5 30];
  nPointsDOWN = [40 70];
  % Generate sample points upstream of area jump (L1)
  x1_up = linspace(-0.1*x1_refL,-0.000001*x1_refL,nPointsUP(1));
  x2_up = linspace(0,p.R_i,nPointsUP(2));
end

[X1_up,X2_up] = meshgrid(x1_up,x2_up);
X_up = X1_up(:) + 1i*X2_up(:);
X_up_mesh = {X1_up,X2_up};
% Generate sample points downstream of area jump (L1)
x1_down = linspace(0.000001*x1_refL,1.3*x1_refL,nPointsDOWN(1));
x2_down = linspace(0,p.R_a,nPointsDOWN(2));
[X1_down,X2_down] = meshgrid(x1_down,x2_down);
X_down = X1_down(:) + 1i*X2_down(:);
X_down_mesh = {X1_down,X2_down};

% Map points to image domain
meshFile = ['output',filesep,'XI_mapped_',p.flameType,...
  num2str(nPointsUP(1)),num2str(nPointsUP(2)),num2str(nPointsDOWN(1)),num2str(nPointsDOWN(2)),'.mat'];
if exist(meshFile,'file')==2
  load(meshFile)
else
  Xi_all{1} =  myMap.xi_x_L1( X_up );
  Xi_all{2} =  myMap.xi_x_L1( X_down );
  save(meshFile,'Xi_all','X_up_mesh','X_down_mesh')
end

end

