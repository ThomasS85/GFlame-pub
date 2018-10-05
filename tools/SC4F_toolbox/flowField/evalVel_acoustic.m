function [ u_ac , u_ac_fun ] = evalVel_acoustic( x , u_in , p , varargin )
%EVALVEL_ACOUSTIC evaluates compact acoustic velocity field (irrotational) at points x (L2)
%
% Inputs:
%   - x      : Position at which velocity should be evaluated 
%              (L2 system; if provided in L1, set option accordingly!)
%   - u_in   : Velocity far upstream of the flame in the feed channel
%   - p      :  Struct with flame settings as returned from setUpPredefinedFlame()
%
% Outputs:
%   - u_ac      : Velocity field due to source at -infinity
%                   (L2 system; if desired in L1, set option accordingly!)
%   - u_ac_fun  : Function handle which evaluates velocity field due to source at -infinity
%                   (L2 system; if desired in L1, set option accordingly!)
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

% Is x already provided in image domain? If so, no mapping is required
ind = find(strcmpi(varargin,'xi'),1);
if ~isempty(ind)
  % x is xi in image domain. Output L2 if not chosen differently!
  isXi = 1;
else
  % x is x in physical domain (L1 or L2 system)
  isXi = 0;
end


%% Map points to image domain
if isXi
  xi = x;
else
  if mapL1L2
    % x provided in L1
    xi =  SCmapInv_SCFT( x , p , 'L1' );
  else
    % x provided in L2
    xi =  SCmapInv_SCFT( x , p );
  end
end

%% Evaluate velocity field in image domain
% Get mapping infos
[ s ] = return_SCmap_SCFT( p );

% Velocity field due to vortexes in physical domain (L1)
if mapL1L2
  % If input was L1, then also return in L1
  u_ac_fun = @(xi) s.ccVelIrr(xi,u_in);
else
  % If input was L2, then also return in L2 (complex conjugate)
  u_ac_fun = @(xi) conj( s.ccVelIrr(xi,u_in) );
end

u_ac = u_ac_fun(xi);


end

