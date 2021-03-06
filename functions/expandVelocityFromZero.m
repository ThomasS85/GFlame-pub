function [ val_expanded ] = expandVelocityFromZero( cFuns , valC , grid , data )
%EXPANDVELOCITYFROMZERO Extends scalar at level-set zero line to the grid in normal direction of zero iso-ine
%   If velocity field shall be expanded off the zero iso-line, in order to maintain signed distance property
%   of G-field, this can either done by expanding the individual velocity components (u_1,u_2) separately or
%   by expanding the absolute value of the flame normal velocity at the zero iso-line! 
%
% Inputs:
%   cFuns - Iso-lines as, e.g., provided by extracIsoLines_SA()
%            { iso-line 1 , iso-line 2 , ... }
%   valC  - cell array, with each cell contains values of the scalar to be extended off the zero iso line
%            { scalar 1 over iso-line 1 , scalar 1 over iso-line 2 , ... ; 
%              scalar 2 over iso-line 1 , scalar 2 over iso-line 2 , ... ;
%              ... }
%   grid  - grid struct
%   data  - Implicit surface function (need not to be signed distance function since information is propagated
%           by normalized velocity:   sign(G) * grad G / |grad G|
%
% Outputs:
%   val_expanded  - scalar field that was expanded from the zero iso-line
%
% by Thomas Steinbacher, Jul 2018
%

%% Pre-process data
% Number of iso -lines
nIso = length(cFuns);
% Number of scalar fields
nScal = length(valC(:,1));

%% Initialize close points
% Find points that are close to zero iso-line
closePoints_bool = isNearInterface(data , 0 , 0);
closePoints = zeros(sum(closePoints_bool(:)),2);
closePoints(:,1) = grid.xs{1}(closePoints_bool); % x1
closePoints(:,2) = grid.xs{2}(closePoints_bool); % x2

% sign of data for all close points
% signData_CP = sign(data(closePoints_bool));

% Loop over all iso-lines to compute initial scalar field val_init
val_init = cell(nScal,1);
for ii=1:nIso
  % Compute for all close points shortest distance to curve C + the corresponding closest point on C
  cPoints = cFuns{ii}';
  [closeP_onC,dist2C,ta] = distance2curve(cPoints,closePoints,'linear');
  
  % Only take points close to the ii-th curve
  isClose_tmp = dist2C<max(grid.dx);
  
  % Define coordinate along flame front
  myDs = [ 0 ; sqrt(sum(diff(cPoints,[],1).^2,2)) ];
  myS = cumsum(myDs);
  [ myS , ia ] = unique( myS/myS(end) ); % normalize + ensure uniquess
  
  % What to add in order to extend scalar fields at slope 1
  % add2CP = dist2C(isClose_tmp) .* signData_CP(isClose_tmp);
  
  % Loop over all provided scalars
  for jj=1:nScal
    % Interpolate valC to closeP_onC
    val_CP = interp1(myS,valC{jj,ii}(ia),ta(isClose_tmp));
    
    % Extend scalar to close points + compute initial val field
    % val_CP = valC_closeP + add2CP;
    
    % Loop over all points -> inefficient! replace!?
    closP_tmp = closePoints_bool(:);
    count = 1;
    for kk=1:length(closP_tmp)
      if closP_tmp(kk)
        closP_tmp(kk) = isClose_tmp(count);
        count = count + 1;
      end
    end
    closP_tmp = reshape(closP_tmp,size(data));
    
    val_init{jj} = zeros(size(data));
    val_init{jj}(closP_tmp) = val_CP;
    
%     % Debug plot: valC along flameshape and for closest points
%     figure;hold on;
%     plot3(cPoints(:,1),cPoints(:,2),valC{jj,ii},'b.');
%     plot3(closeP_onC(isClose_tmp,1),closeP_onC(isClose_tmp,2),val_CP,'ro')
%     xlabel('x1');ylabel('x2')
%     zlabel('scalar')
    
  end
  
%   % Debug plot: flame shape + closest points on flame
%   figure;hold on; plot(cPoints(:,1),cPoints(:,2),'b.'); plot(closeP_onC(:,1),closeP_onC(:,2),'ro')
%   legend({'points of 0 iso-line','points that are closest to grid points'})
  
  
end


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

% % Debug plot: velocity field (should point away from zero level set
% figure; hold on;
% contour(grid.vs{1}, grid.vs{2}, data', [0,0])
% quiver(grid.xs{1}(:),grid.xs{2}(:),effectiveVelocity{1}(:),effectiveVelocity{2}(:))

% Debug computation: check nor of velocity
% norm_effVel = sqrt( effectiveVelocity{1}.^2 + effectiveVelocity{2}.^2 );


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

% Choose time interval: THIS SETTING AFFECTS COMPUTATIONAL EFFORT MOST 
%   -> If t_max set to 1/2 of max(grid.max - grid.min) this reduces computational costs by factor ~2!
tSpan = [ 0 max(grid.max - grid.min) / 2 ];

% Solve
val_expanded = cell(nScal,1);
for jj=1:nScal
  [ ~ , val_expanded{jj} ] = feval(integratorFunc, schemeFunc, tSpan, val_init{jj}(:),...
    integratorOptions, schemeData);
  
  val_expanded{jj} = reshape(val_expanded{jj},size(data));
end

%% Debug plots

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

