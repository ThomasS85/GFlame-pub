function [ vortDat ] = RK4_Vorticity_pointSources_step( vel , dt , varargin )
%RK4_POINTSOURCES_STEP Performes one single RK4 step for a given vortex/ source distribution

%% Parse varargin
% Should mapping be applied? Only relevant for r0!
ind = find(strcmpi(varargin,'myMapping'),1);
if ~isempty(ind)
  % Use user specified mapping
  myMap = varargin{ind+1};
else
  % Default: No mapping
  myMap.dxi_dx = @(xi) 1;
  myMap.xi_x = @(x) x;
  myMap.x_xi = @(xi) xi;
end

% Suppress application of Routh's correction
ind = find(strcmpi(varargin,'noRouth'),1);
if ~isempty(ind)
  % No Routh's correction
  doRouthsCorr = 0;
else
  % Default: Do, if available
  doRouthsCorr = 1;
end


%% Check
 % Vortex gets too close to wall! Prevent it from getting closer
if imag(vel.vortDat.xi) < 1e-7
  vel.vortDat.xi = real(vel.vortDat.xi) +1i*1e-6;
  vel.vortDat.x = myMap.x_xi(vel.vortDat.xi);
end

% Should Routh's correction be applied?
if doRouthsCorr
  if isfield(myMap,'RouthsCorr') && vel.vortDat.r0==0
    doRouthsCorr = 1;
  else
    doRouthsCorr = 0;
  end
end


%% Init
vortDat = vel.vortDat;


%% Classical Runge-Kutta 4th order
if doRouthsCorr 
  % Apply Routh's correction
  % 1
  k1 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) ...
    .* conj( myMap.dxi_dx(vel.vortDat.xi) + myMap.RouthsCorr(vel.vortDat.G,vel.vortDat.xi) );
  % 2
  vel.vortDat.x = vortDat.x + k1/2;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k2 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) ...
    .* conj( myMap.dxi_dx(vel.vortDat.xi) + myMap.RouthsCorr(vel.vortDat.G,vel.vortDat.xi) );
  % 3
  vel.vortDat.x = vortDat.x + k2/2;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k3 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) ...
    .* conj( myMap.dxi_dx(vel.vortDat.xi) + myMap.RouthsCorr(vel.vortDat.G,vel.vortDat.xi) );
  % 4
  vel.vortDat.x = vortDat.x + k3;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k4 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) ...
    .* conj( myMap.dxi_dx(vel.vortDat.xi) + myMap.RouthsCorr(vel.vortDat.G,vel.vortDat.xi) );
  
else
  % Do not apply Routh's correction
  % 1
  k1 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) .* conj( myMap.dxi_dx(vel.vortDat.xi) );
  % 2
  vel.vortDat.x = vortDat.x + k1/2;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k2 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) .* conj( myMap.dxi_dx(vel.vortDat.xi) );
  % 3
  vel.vortDat.x = vortDat.x + k2/2;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k3 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) .* conj( myMap.dxi_dx(vel.vortDat.xi) );
  % 4
  vel.vortDat.x = vortDat.x + k3;
  vel.vortDat.xi = myMap.xi_x( vel.vortDat.x );
  k4 = dt * evalFlowField_PointSource( vel.vortDat.xi , vel , 'myMapping' , myMap ) .* conj( myMap.dxi_dx(vel.vortDat.xi) );

end

% Step!
vortDat.x = vortDat.x + ( k1/6 + k2/3 + k3/3 + k4/6 );
vortDat.xi = myMap.xi_x( vortDat.x );

if any(isnan(vortDat.xi))
  error('Vortex got too close to vertex which is mapped to infinity! Matrix is singular.')
end

end

