function [ Q_int ] = Q_exact( x1L_flame_s , s_L_flame_s , p )
%Q_FROM_XI Evaluates the heat releas rate of a perturbed flame normalized by density rho and heat of reaction
%delta h_r
%
% Inputs:
%   - x1L_flame_s      : Perturbed flame position in laboratory coordinates (L1)
%   - s_L_flame_s      : Flame speed distribution over flame
%   - p                : Flame struct
%
%   Outputs:
%           - Q_int         : Integral heat release rate devided by (rho* delta h_r)
%                               Density times heat of reaction
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 10.11.2017                                //
% // Last modified: 10.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////

%% Parse varargin

%% Flame aligned coordinates
% Check input dimeniosn and correct if required
[z,s] = size(x1L_flame_s);
if z>s; x1L_flame_s = x1L_flame_s.'; end;
[z,s] = size(s_L_flame_s);
if z>s; s_L_flame_s = s_L_flame_s.'; end;

% Compute flame coordinate
myCurve = [ real(x1L_flame_s) ; imag(x1L_flame_s) ];
myDs = sqrt(sum(diff(myCurve,[],2).^2,1));
myDs = [0, myDs]; % add starting point
myS = cumsum(myDs);

%% Evaluuate heat release
if strcmp(p.Fdim,'2D')
  % Slit flames
  Q_int = trapz(myS,s_L_flame_s);
  
elseif strcmp(p.Fdim,'3D')
  % Conical flames
  
  % local radius
  my_rad = imag(x1L_flame_s);
  % Evaluate surface integral
  Q_int = 2*pi*trapz( myS , my_rad.*s_L_flame_s );
end

end

