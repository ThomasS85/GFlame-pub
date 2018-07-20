function [ val_expanded ] = expandVelocityFromZero( cFuns , valC , grid , data )
%EXPANDVELOCITYFROMZERO Extends velocity at level-set zero line to the grid such that signed distance property
% is maintained by computing convection
%
% Inputs:
%   cFuns - Iso-lines as, e.g., provided by extracIsoLines_SA()
%   valC  - Values of the scalar to be extended off the zero iso line
%   grid  - grid struct
%   data  - Implicit surface function
%
% Outputs:
%
%
% by Thomas Steinbacher, Jul 2018

%% Initialize close points
% Find points that are close to zero iso-line
closePoints_bool = isNearInterface(data , 0 , 0);
closePoints = zeros(sum(closePoints_bool(:)),2);
closePoints(:,1) = grid.xs{1}(closePoints_bool); % x1
closePoints(:,2) = grid.xs{2}(closePoints_bool); % x2

% sign of data for all close points
signData_CP = sign(data(closePoints_bool));

% Compute for all close points shortest distance to curve C + the corresponding closest point on C
cPoints = cFuns{1}';
[closeP_onC,dist2C,ta] = distance2curve(cPoints,closePoints,'linear');

% Define coordinate along flame front
myDs = sqrt(sum(diff(cPoints,[],1).^2,2));
myDs = [0; myDs]; % add starting point
myS = cumsum(myDs);
myS = myS/myS(end); % normalize

% Interpolate valC to closeP_onC 
valC_closeP = interp1(myS,valC,ta);


%% Extend scalar to close points + compute initial val field
val_CP = valC_closeP + dist2C .* signData_CP;
val_init = zeros(size(data));
val_init(closePoints_bool) = val_CP;


%% Evaluate sign(data) * (grad data) / |grad data|  -> effectiveVelocity
speed = sign(data);
magnitude = zeros(size(data));
effectiveVelocity = cell(2,1);
for ii=1:grid.dim
  % upwindFirstWENO5  or upwindFirstENO2  or upwindFirstFirst
  [ derivL, derivR ] = upwindFirstFirst(grid, data, ii);
  
  % Effective velocity in this dimension (scaled by \|\grad \phi\|).
  prodL = speed .* derivL;
  prodR = speed .* derivR;
  magL = abs(prodL);
  magR = abs(prodR);
  
  % Determine the upwind direction.
  %   Either both sides agree in sign (take direction in which they agree),
  %   or characteristics are converging (take larger magnitude direction).
  flowL = ((prodL >= 0) & (prodR >= 0)) | ...
    ((prodL >= 0) & (prodR <= 0) & (magL >= magR));
  flowR = ((prodL <= 0) & (prodR <= 0)) | ...
    ((prodL >= 0) & (prodR <= 0) & (magL < magR));
  
  % Now we know the upwind direction, add its contribution to \|\grad \phi\|.
  magnitude = magnitude + derivL.^2 .* flowL + derivR.^2 .* flowR;
  
  % CFL condition: sum of effective velocities from O&F (6.2).
  effectiveVelocity{ii} = prodL .* flowL + prodR .* flowR;

end

% Finally, calculate speed * \|\grad \phi\|
magnitude = sqrt(magnitude);

% Normalize fields
effectiveVelocity{1} = effectiveVelocity{1} ./ magnitude;
effectiveVelocity{2} = effectiveVelocity{2} ./ magnitude;


%% Now solve PDE to propagate initial velocity field
% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = effectiveVelocity;
schemeData.grid = grid;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.95, 'stats', 'off');

% Choose approximations at appropriate level of accuracy.
schemeData.derivFunc = @upwindFirstFirst;
integratorFunc = @odeCFL1;

% Choose time interval
tSpan = [ 0 max(grid.max - grid.min) ];

% Solve
[ ~ , val_expanded ] = feval(integratorFunc, schemeFunc, tSpan, val_init(:),...
  integratorOptions, schemeData);

val_expanded = reshape(val_expanded,size(data));


%% Debug plots
% % flame shape + closest points on flame
% figure;hold on; plot(cPoints(:,1),cPoints(:,2),'b.'); plot(closeP_onC(:,1),closeP_onC(:,2),'ro')
% % valC along flameshape and for closest points
% figure;hold on;
% plot3(cPoints(:,1),cPoints(:,2),valC,'b.');
% plot3(closeP_onC(:,1),closeP_onC(:,2),valC_closeP,'ro')
% xlabel('x1');ylabel('x2')
% zlabel('scalar')
% % plot gradient magnitude field of data
% figure;
% surf(grid.xs{1},grid.xs{2},magnitude)
% xlabel('x1');ylabel('x2')
% zlabel('magnitude')
% % plot gradient field of data
% figure;
% surf(grid.xs{1},grid.xs{2},effectiveVelocity{1})
% xlabel('x1');ylabel('x2')
% zlabel('eff. Vel. x1')
% figure;
% surf(grid.xs{1},grid.xs{2},effectiveVelocity{2})
% xlabel('x1');ylabel('x2')
% zlabel('eff. Vel. x2')
% figure
% myVelNorm = sqrt(effectiveVelocity{1}.^2 + effectiveVelocity{2}.^2);
% surf(grid.xs{1},grid.xs{2},myVelNorm)
% xlabel('x1');ylabel('x2')
% zlabel('||eff. Vel||_2')
% 
% a = 0;

end
