function [ u_vort , dOmega_dxi ] = evalVel_vortex( x , xs , Gamma , p , varargin )
%EVALVEL_VORTEX evaluates the velocity due to a vortex of strength Gamma placed at the position
%  x_s (L2) at the points x (L2). Zero flux boundary conditions are met by introducing a mirror vortex in
%  image domain
%
% !!!!!!!!!!!!
% Function only required to be compatible with Symposium paper! 
% otherwise use evalFlowField_physicalDomain() instead!
% !!!!!!!!!!!!
%
% Inputs:
%   - x      : Position at which velocity should be evaluated 
%              (L2 system; if provided in L1, set option accordingly!)
%   - xs     : Position of vortexes in L1 system; can be a vector 
%              (L2 system; if provided in L1, set option accordingly!)
%   - Gamma  : Strength of vortex
%   - p         :  Struct with flame settings as returned from setUpPredefinedFlame()
%
% Outputs:
%   - u_vort      : Velocity field due to vortexes at points x
%                   (L2 system; if desired in L1, set option accordingly!)
%   - dOmega_dxi  : Complex conjugate of velocity field in image domain
%
%
% % Coordinate System L2:
%    ^
% x2 |__ . __ . __ . __ . __
%    |          o             x_L(1)
%    |       o                  :
%    |--- o                     :
%    |   |                    x_L(end)
%  0 ---------------------> 
%        0               x1
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 09.11.2017                                //
% // Last modified: 09.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////


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

% Use Lamb-Oseen Vortex instead of point vortex?  {r_0}
ind = find(strcmpi(varargin,'lamb'),1);
if ~isempty(ind)
  % Use lamb vortex
  vortparams = varargin{ind+1};
else
  % Default: Use point vortex
  vortparams = {0};
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
  doKutta = 'doKutta';
else
  % Default: No Kutta condition
  doKutta = 'noKutta';
end


%% Map points to image domain
if isXi
  xi = x;
  xi_s = xs;
else
  if mapL1L2
    % x provided in L1
    xi =  SCmapInv_SCFT( x , p , 'L1' );
    xi_s =  SCmapInv_SCFT( xs , p , 'L1' );
  else
    % x provided in L2
    xi =  SCmapInv_SCFT( x , p );
    xi_s =  SCmapInv_SCFT( xs , p );
  end
end


%% Evaluate velocity field in image domain
% Get mapping infos
[ s ] = return_SCmap_SCFT( p );

% Wrap data to fit evalFlowField_PointSource() inputs
vel.vortDat.G = Gamma;
vel.vortDat.xi = xi_s;
vel.vortDat.r0 = vortparams{1} * ones(size(Gamma));
vel.sourceDat = [];
vel.u_p = 0;
vel.doMirror = 'mirror';


u_vort = evalFlowField_physicalDomain( xi , vel , s , doKutta , 1 , 'xi' );
dOmega_dxi = conj(u_vort);

% Velocity field due to vortexes in physical domain (L1)
if mapL1L2
   u_vort = dOmega_dxi;
end



% % Evaluate velocity field in image domain
% [ ~ , dOmega_dxi ] = evalFlowField_PointSource( xi , vel , 'myMapping' , s );
% 
% % Routh's correction
% RouthsCorr = zeros(size(xi));
% if isfield(s,'RouthsCorr') && vortparams{1}==0
%   for ii=1:length(Gamma)
%     myInd_xi_s = xi==xi_s(ii);
%     RouthsCorr(myInd_xi_s) = s.RouthsCorr( Gamma(ii) ,xi(myInd_xi_s) );
%   end
% end
%       
% % Velocity field due to vortexes in physical domain (L1)
% if mapL1L2
%   % If input was L1, then also return in L1
%   u_vort = s.dxi_dx(xi) .* dOmega_dxi + RouthsCorr;
% else
%   % If input was L2, then also return in L2 (complex conjugate)
%   u_vort = conj( s.dxi_dx(xi) .* dOmega_dxi + RouthsCorr );
% end


% doLamb = 0;
% if doLamb
%   % Remove points which are too close at trailing edge
%   vort2Remove = abs(xi_s-1)<1e-1;
%   xi_s(vort2Remove) = 2;
%   Gamma(vort2Remove) = 0;
%   % Radius depends on xi value since length are not preserved in image domain! We assume dxi_dx is
%   % approximately const. everywhere inside the vortexcore radius
%   r_0 = (vortparams{1} * abs(s.dxi_dx(xi_s))/s.l_ref);
%   dOmega_dxi = @(xi,xi_s,Gamma,myR_0) myLambOseenVortex(xi,xi_s,Gamma,myR_0);
%   %   dOmega_dxi = @(xi,xi_s,Gamma) 1i * Gamma/(2*pi) * ( ( xi_s - conj(xi_s) ) ./ ...
%   %     ( ( xi - xi_s ) .* ( xi - conj(xi_s) ) ) ...
%   %     + ( exp(-abs(xi-conj(xi_s)).^2/r_0^2)./(xi-conj(xi_s)) - exp(-abs(xi-xi_s).^2/r_0^2)./(xi-xi_s) ) );
%   %   dOmega_dxi = @(xi,xi_s,Gamma) 1i * Gamma/(2*pi) * 1 ./ ( xi - xi_s ) .* ...
%   %      ( 1 - exp(-abs(xi-xi_s).^2/r_0^2) );
% else
%   % Velocity field due to point vortex and mirror vortex in image domain
%   dOmega_dxi = @(xi,xi_s,Gamma,myR_0) 1i * Gamma/(2*pi) * ( xi_s - conj(xi_s) ) ./ ...
%     ( ( xi - xi_s ) .* ( xi - conj(xi_s) ) );
%   r_0 = zeros(size(xi_s));
% end
% 
% % Velocity field due to vortexes in physical domain (L1)
% if mapL1L2
%   % If input was L1, then also return in L1
%   u_vort_fun = @(xi,xi_s,Gamma,myR_0) s.dxi_dx(xi) .* dOmega_dxi(xi,xi_s,Gamma,myR_0);
% else
%   % If input was L2, then also return in L2 (complex conjugate)
%   u_vort_fun = @(xi,xi_s,Gamma,myR_0) conj( s.dxi_dx(xi) .* dOmega_dxi(xi,xi_s,Gamma,myR_0) );
% end
% 
% % Compute velocity due to all vortexes
% u_vort = zeros(size(x));
% for ii=1:length(xi_s)
%   if length(Gamma) == length(xi_s)
%     u_vort = u_vort + u_vort_fun( xi , xi_s(ii) , Gamma(ii) , r_0(ii) );
%   else
%     u_vort = u_vort + u_vort_fun( xi , xi_s(ii) , Gamma(1) , r_0(ii) );
%   end
% end


end

