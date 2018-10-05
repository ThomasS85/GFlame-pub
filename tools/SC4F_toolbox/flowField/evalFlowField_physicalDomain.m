function [ u_res , vel , dOmega_dxi_Kutta ] = evalFlowField_physicalDomain( myX , vel , myMap , varargin )
%EVALFLOWFIELD Evaluates the flow field, resulting from a number of vortexes at x_v, a number of sources
%   x_s and constant part u_p in the physical x plane. Allows application of a Kutta-condition by application of
%   a single panel method (see Disselhorst and van Wijngaarden 1980)
%
% Inputs:
%    myX            - Position at which velocity should be evaluated
%                     (L2 system; if provided in L1, set option accordingly! Also possible to directly provide Xi!)
%    vel.vortDat    - Struct which contains information about all vortexes
%                     .xi  : Nx1 complex array with positions of all vortexes [ 1+1i*2 ; 2+1i*3 ]
%                     .G   : Nx1 real array with circulations of each vortex [ 2 ; 4 ]
%                     .r0  : Nx1 real array with kernel width of each vortex [ 3 ; 4]
%                             If set to zero, a point vortex is assumed, otherwise a Lamb vortex
%    vel.sourceDat  - Struct which contains information about all (point) sources
%                     .xi  : Nx1 complex array with positions of all sources [ 1+1i*2 ; 2+1i*3 ]
%                     .G   : Nx1 real array with strengths of each source [ 2 ; 4 ]
%    vel.u_p        - Complex constant velocity which is added to each point
%    vel.Kutta      - Information about Kutta condition
%                     .G   : Circulation of vortex panel
%                     .H   : Length of the vortex panel (in physical domain)
%                     .beta: shear layer angle
%                     .u   : Velocity right at Kutta point
%    myMap          - Struct with mapping information as returned from return_SCmap_SCFT()
%
% outputs:
%    u_res          - Velocity (physical domain L2, if specified otherwise then L1) at all specified x-locations
%    vel            - Struct with velocity field information; in particular Kutta data might be changed
%                     compared to input!
%    dOmega_dxi_Kutta - Potential due to Kutta pannel
%
%
% by Thomas Steinbacher (03.2018)

%% Settings
% Choose vorticity distribution in vortex sheet: (visuallly case 2 looks better)
%   1: sqrt(xi-1)
%   2: sqrt(xi-1) / xi
KuttaCase = 2;

%% Parse varargin
% Is x provided in L1 coordinate system? If so, first transform to L2!
ind = find(strcmpi(varargin,'L1'),1);
if ~isempty(ind)
  % Do map from L1 to L2
  mapL1L2 = 1;
else
  % Don not map to L2; x already in L2!
  mapL1L2 = 0;
end

% Is x already provided in image domain? If so, no mapping is required
ind = find(strcmpi(varargin,'xi'),1);
if ~isempty(ind)
  % x/xs is xi/xi_s in image domain
  isXi = 1;
else
  % x/xs is x/xs in physical domain (L1 or L2 system)
  isXi = 0;
end

% Meet Kutta condition? If so, the prevertex coordinates of where the Kutta condition should be met has to be
%  specified
ind = find(strcmpi(varargin,'doKutta'),1);
if ~isempty(ind)
  % Meet Kutta condition
  doKutta = 1;
  xiKutta = varargin{ind+1};
else
  % No Kutta condition
  doKutta = 0;
end

% Get length of desired time step and mean flow jet velocity in order to compute length of shear layer panel.
%  If not specified, the length already stored in vel.Kutta.H is used
%  or, if not existent, a default panel length is assumed
ind = find(strcmpi(varargin,'KuttaJet'),1);
if ~isempty(ind)
  % Get length of time step
  dt = varargin{ind+1};
  % Get mean flow velocity of jet
  u_jet = varargin{ind+2};
else
  % Default: No time step length and no mean jet velocity -> Use default length
  dt = [];
  u_jet = [];
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


% Should vorticity of the panel NOT be forwarded to a point vortex?
ind = find(strcmpi(varargin,'noPanel2Point'),1);
if ~isempty(ind)
  % No!
  doPanel2Point = 0;
else
  % Default: Yes
  doPanel2Point = 1;
end


%% Preprocessing
% Should Routh's correction be applied?
if doRouthsCorr
  if isfield(myMap,'RouthsCorr') && ~isempty(vel.vortDat)
    if vel.vortDat.r0==0
      doRouthsCorr = 1;
    else
      doRouthsCorr = 0;
    end
  else
    doRouthsCorr = 0;
  end
end


%% Map points to image domain
if isXi
  myXi = myX;
else
  if mapL1L2
    % x provided in L1
    myXi =  myMap.xi_x_L1(myX);
  else
    % x provided in L2
    myXi =  myMap.xi_x(myX);
  end
end


%% Apply Kutta condition
if doKutta
  % If there is already an active Kutta panel, propagate this to a point vortex
  if isfield(vel,'Kutta') && doPanel2Point
    if vel.Kutta.active && vel.Kutta.H*vel.Kutta.G~=0
      if isfield(vel.vortDat,'G') && isfield(vel.vortDat,'xi')
        [ myCirc , myXi_s ] = panel2PointVortex( vel , vel.Kutta.beta , KuttaCase , xiKutta );
        vel.vortDat.G = [ myCirc ; vel.vortDat.G ]; % attach to position 1
        vel.vortDat.x = [ myMap.x_xi( myXi_s ) ; vel.vortDat.x ];
        vel.vortDat.xi = [ myMap.xi_x( vel.vortDat.x(1) ) ; vel.vortDat.xi ];
        vel.vortDat.r0 = [ 0 ; vel.vortDat.r0 ];
      else
        % Compute circulation of vortex street by assuming a square root vorticity distribution
        [ myCirc , myXi_s ] = panel2PointVortex( vel , vel.Kutta.beta , KuttaCase , xiKutta );
        vel.vortDat.G = myCirc;
        vel.vortDat.x = myMap.x_xi( myXi_s );
        vel.vortDat.xi = myMap.xi_x( vel.vortDat.x(1) );
        vel.vortDat.r0 = 0;
      end
      % Set panel inactive
      vel.Kutta.active = 0;
    end
  end
  
  % Evaluate velocity and potential at Kutta point xiKutta
  [ ~ , dOmega_dxi_xiKutta_sum ] = evalFlowField_PointSource( xiKutta , vel , 'myMapping' , myMap );
  
  % Angle functions (dOmega_dxi_xiKutta_sum should be real due to no flux BC!)
  % Compute corret velocity in order to evalaute angle
  if dOmega_dxi_xiKutta_sum > 0 
    % Separation downstream
    MyBeta = pi/3;
  else
    % Separation upstream
%     MyBeta = pi*2/3;
    MyBeta = pi/3;
  end
  
  
  % Set angle vortex sheet (image domain)
  vel.Kutta.beta = MyBeta;
  
  % Compute panel length H:
  if ( ~isempty(dt) && ~isempty(u_jet) )
    % Assume only mean flow velocity of jet responsible for advective transport of vorticity
    %  -> Vorticity is transported by half the jet velocity + local perturbation velocity upK=(U^+ + U^-)/2
    if isfield(vel.Kutta,'H')
      % If Kutta pannel exisits, evaluate limit of velocity imposed by pannel right at separation point
      limLeft = -1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
        1i*pi*sin(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) );
      limRight = 1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
        pi*cos(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) );
      limCent = ( limLeft + limRight ) /2;
      vel.Kutta.H = imag( myMap.xi_x( myMap.x_xi( xiKutta ) + (u_jet/2+limCent)*dt ) );
    else
      % If no pannel exists, use mean flow jet velocity to compute H (effect of sources/ vortexes on edge
      %  should then be 0 anyway!)
      vel.Kutta.H = imag( myMap.xi_x( myMap.x_xi( xiKutta ) + (u_jet/2)*dt ) );
    end
  else
    % If field Kutta exists, take length specified here. Otherwise, compute length from l_ref
    if ~isfield(vel.Kutta,'H')
      % Derive H from reference length
      vel.Kutta.H = abs( myMap.dxi_dx( myMap.xi_x( myMap.x_xi(xiKutta) + myMap.l_ref * 0.6 ) ) ) * myMap.l_ref*0.6;
    end
    
  end
  
  % Evaluate potential vortex sheets
  switch KuttaCase
    case 1
      % sqrt(xi-1)
      % (1) Compute Circulation vortex panel G (from dPhi_sum/dxi + dPhi_panel/dxi = 0 -> Kutta condition)
      vel.Kutta.G = dOmega_dxi_xiKutta_sum * pi /...
        ( 1i * sqrt(vel.Kutta.H ) * ( exp(-1i*vel.Kutta.beta) -  exp(1i*vel.Kutta.beta) ) );
      % (2) Now evaluate velocity due to vortex panel at all points myXi
      dOmega_dxi_Kutta_up   = @(xi)  1i / pi * vel.Kutta.G * ( sqrt(vel.Kutta.H) - sqrt((1-xi)*exp(-1i*vel.Kutta.beta)) .*...
        atan( sqrt( vel.Kutta.H ./ ((1-xi)*exp(-1i*vel.Kutta.beta)) ) ) ) * exp(-1i*vel.Kutta.beta);
      dOmega_dxi_Kutta_down   = @(xi)  -1i / pi * vel.Kutta.G * ( sqrt(vel.Kutta.H) - sqrt((1-xi)*exp(1i*vel.Kutta.beta)) .*...
        atan( sqrt( vel.Kutta.H ./ ((1-xi)*exp(1i*vel.Kutta.beta)) ) ) ) * exp(1i*vel.Kutta.beta);
      dOmega_dxi_Kutta = @(xi) dOmega_dxi_Kutta_up(xi) + dOmega_dxi_Kutta_down(xi);
    case 2
      % sqrt(xi-1) / xi
      % (1) Compute Circulation vortex panel G (from dPhi_sum/dxi + dPhi_panel/dxi = 0 -> Kutta condition)
      vel.Kutta.G = -dOmega_dxi_xiKutta_sum * pi /...
        ( 1i * atan(sqrt(vel.Kutta.H )) * ( exp(-1i*vel.Kutta.beta) -  exp(1i*vel.Kutta.beta) ) );
      % (2) Now evaluate velocity due to vortex panel at all points myXi
      dOmega_dxi_Kutta_up   = @(xi) 1i / pi * vel.Kutta.G ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) .* ( atan( sqrt(vel.Kutta.H) ) - ...
        sqrt( (1-xi)*exp(-1i*vel.Kutta.beta) ) .* atan( sqrt( vel.Kutta.H ./ ( (1-xi)*exp(-1i*vel.Kutta.beta) ) ) ) ) *exp(-1i*vel.Kutta.beta);
      dOmega_dxi_Kutta_down   = @(xi) -1i/pi * vel.Kutta.G ./( (xi-1)*exp(1i*vel.Kutta.beta)+1 ) .* ( atan( sqrt(vel.Kutta.H) ) - ...
        sqrt( (1-xi)*exp(1i*vel.Kutta.beta) ) .* atan( sqrt( vel.Kutta.H ./ ( (1-xi)*exp(1i*vel.Kutta.beta) ) ) ) ) *exp(1i*vel.Kutta.beta);
      dOmega_dxi_Kutta = @(xi) dOmega_dxi_Kutta_up(xi) + dOmega_dxi_Kutta_down(xi);
  end
  
  % Set Kutta panel active
  vel.Kutta.active = 1;
  
else
  % No contribution due to Kutta condition
  dOmega_dxi_Kutta = @(xi) 0;
  
end


%% Apply Routh's rule (only right at point vortexes)
u_Routh = zeros(size(myXi));
if doRouthsCorr
  for ii=1:length(myXi)
    indsTem = myXi(ii)==vel.vortDat.xi;
    if any(indsTem)
      u_Routh(ii) = myMap.RouthsCorr( vel.vortDat.G(indsTem) , vel.vortDat.xi(indsTem) );
    end
  end
end


%% Evaluate resulting velocity
% only evaluate vortices/ sources of non-zero strength    -> Implement?!
% nonZeroS = vel.sourceDat.G ~= 0;
% nonZeroV = vel.vortDat.G ~= 0;
% vel_tmp.sourceDat.G = vel.sourceDat.G(nonZeroV);
%  ...

% evaluate!
[ ~ , dOmega_dxi ] = evalFlowField_PointSource( myXi , vel , 'myMapping' , myMap );
u_res =  conj( ( dOmega_dxi + dOmega_dxi_Kutta(myXi) ) .* myMap.dxi_dx( myXi ) ) ...
  + conj( u_Routh );

if doKutta
  % Set velocity at Kutta condition point to analytical determined limit
  myXi(abs(xiKutta-myXi)<5e-6) = xiKutta;
  u_res(abs(xiKutta-myXi)<=1e-5 & real(xiKutta-myXi)>=0 ) = -1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
    1i*pi*sin(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) );
  u_res(abs(xiKutta-myXi)<=1e-5 & real(xiKutta-myXi)<0 ) = 1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
    pi*cos(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) );
end


% %% Debug limit
%
% % limit
% limLeft = real( -1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
%       1i*pi*sin(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) ) );
% limRight = real( 1/myMap.l_ref * sqrt(1-myMap.h)*dOmega_dxi_xiKutta_sum *...
%       pi*cos(3/2*vel.Kutta.beta) / ( 2*atan(sqrt(vel.Kutta.H)) * sin(vel.Kutta.beta) ) );
% limCent = ( limLeft + limRight ) /2;
%
% % Plot Kutta condition as used above slightly above real axis
% myXi2 = 1i*1e-5 + 1+linspace(-5e-1,5e-1,1000);
% [~,dOdxi_2] = evalFlowField_PointSource( myXi2 , vel , 'myMapping' , myMap );
% myV2 =  conj( ( dOdxi_2 + dOmega_dxi_Kutta(myXi2) ) .* myMap.dxi_dx( myXi2 ) );
% f1 = figure; hold on;
% plot(real(myXi2),real(myV2),':b');
% plot(1,limLeft,'or');plot(1,limRight,'or');plot(1,limCent,'xr');
% f2 = figure; hold on;
% plot(real(myXi2),imag(myV2),':b')
%
% % Implement formula Kutta condition 1
% dOdxi_Kut = @(xi) dOmega_dxi_xiKutta_sum * (...
%   1/(2*1i*sin(vel.Kutta.beta)) * ( 1./(xi-1+exp(1i*vel.Kutta.beta)) - 1./(xi-1+exp(-1i*vel.Kutta.beta)) ) ...
%   + sqrt(1-xi)./( atan(sqrt(vel.Kutta.H)) * 2*1i*sin(vel.Kutta.beta) ) .* (...
%      exp( 1i*vel.Kutta.beta/2)./(xi-1+exp(-1i*vel.Kutta.beta)) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp( 1i*vel.Kutta.beta)) ))...
%     -sign(real(1-xi)).*exp(-1i*vel.Kutta.beta/2)./(xi-1+exp( 1i*vel.Kutta.beta)) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp(-1i*vel.Kutta.beta)) ))...
%     )...
%   );
% myV_3 = conj( myMap.dxi_dx(myXi2) .* ( dOdxi_2 + dOdxi_Kut(myXi2) ) );
% figure(f1);plot(real(myXi2),real(myV_3),'-.g');
% figure(f2);plot(real(myXi2),imag(myV_3),'-.g');
%
% % Implement formula Kutta condition 2
% dOdxi_Kut2 = @(xi) dOmega_dxi_xiKutta_sum / ( atan(sqrt(vel.Kutta.H))*2*1i*sin(vel.Kutta.beta) ) .* (...
%   exp(-1i*vel.Kutta.beta) ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) .* (...
%     atan(sqrt(vel.Kutta.H)) - sqrt((1-xi)*exp(-1i*vel.Kutta.beta)) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp(-1i*vel.Kutta.beta)) )) )...
%  -exp( 1i*vel.Kutta.beta) ./ ( (xi-1)*exp( 1i*vel.Kutta.beta)+1 ) .* (...
%     atan(sqrt(vel.Kutta.H)) - sqrt((1-xi)*exp( 1i*vel.Kutta.beta)) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp( 1i*vel.Kutta.beta)) )) )...
%     );
% myV_4 = conj( myMap.dxi_dx(myXi2) .* ( dOdxi_2 + dOdxi_Kut2(myXi2) ) );
% figure(f1);plot(real(myXi2),real(myV_4),'--c');
% figure(f2);plot(real(myXi2),imag(myV_4),'--c');
%
% % Implement formula Kutta condition 3
% dOdxi_Kut3 = @(xi) dOmega_dxi_xiKutta_sum / ( atan(sqrt(vel.Kutta.H))*2*1i*sin(vel.Kutta.beta) ) .* (...
%   atan(sqrt(vel.Kutta.H)) * (  exp(-1i*vel.Kutta.beta) ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) - ...
%                                exp( 1i*vel.Kutta.beta) ./ ( (xi-1)*exp( 1i*vel.Kutta.beta)+1 )...
%                             )...
%  +sqrt(1-xi) .* ( ...
%     exp( 3/2*1i*vel.Kutta.beta) ./ ( (xi-1)*exp( 1i*vel.Kutta.beta)+1 ) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp( 1i*vel.Kutta.beta)) ))...
%    -sign(real(1-xi)).*exp(-3/2*1i*vel.Kutta.beta) ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp(-1i*vel.Kutta.beta)) ))...
%    )  );
% % dOdxi_Kut3 = @(xi) dOmega_dxi_xiKutta_sum / ( atan(sqrt(vel.Kutta.H))*2*1i*sin(vel.Kutta.beta) ) .* (...
% %   atan(sqrt(vel.Kutta.H)) * (  exp(-1i*vel.Kutta.beta) ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) - ...
% %                                exp( 1i*vel.Kutta.beta) ./ ( (xi-1)*exp( 1i*vel.Kutta.beta)+1 )...
% %                             )...
% %  +( ...
% %    sqrt((1-xi)*exp( 1i*vel.Kutta.beta)).*exp( 1i*vel.Kutta.beta) ./ ( (xi-1)*exp( 1i*vel.Kutta.beta)+1 ) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp( 1i*vel.Kutta.beta)) ))...
% %   -sqrt((1-xi)*exp(-1i*vel.Kutta.beta)).*exp(-1i*vel.Kutta.beta) ./ ( (xi-1)*exp(-1i*vel.Kutta.beta)+1 ) .* atan(sqrt(vel.Kutta.H ./ ((1-xi)*exp(-1i*vel.Kutta.beta)) ))...
% %    )  );
% myV_5 = conj( myMap.dxi_dx(myXi2) .* ( dOdxi_2 + dOdxi_Kut3(myXi2) ) );
% figure(f1);plot(real(myXi2),real(myV_5),':y','LineWidth',2);
% figure(f2);plot(real(myXi2),imag(myV_5),':y','LineWidth',2);
%
% figure(f1);legend({'Kutta original','limit','Kutta 1', 'Kutta 2' , 'Kutta 3'})


%% If input was L1, also map output back to L1
if mapL1L2
  u_res = conj(u_res);
end

end


function [ myCirc , myXi_s ] = panel2PointVortex( vel , myBeta , KuttaCase , xiKutta )
% Compute circulation of vortex street + center of vorticity
switch KuttaCase
  case 1
    % sqrt(xi-1)
    myCirc =  vel.Kutta.G * 2/3 *vel.Kutta.H^(3/2);
  case 2
    % sqrt(xi-1) / xi
    myCirc =  vel.Kutta.G * integral( @(t) sqrt( t ./ ( 1 + 2*cos(myBeta)*t + t.^2 )) , 0 , vel.Kutta.H );
end

myXi_s = xiKutta + exp(1i*myBeta)*vel.Kutta.H;

end
