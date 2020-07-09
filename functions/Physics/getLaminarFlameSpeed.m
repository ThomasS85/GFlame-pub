function [ s_l_u ] = getLaminarFlameSpeed( mixture , Phi , T_in , p_in )
%GETLAMINARFLAMESPEED Function computes laminar flame speed based on equivalence ratio, Temperature of the
% fresh gas and pressure
%
% Inputs:
%   mixture     : Name of the mixture
%   Phi         : Equivalence ratio (scalar or matrix) [-]
%   T_in        : Temperature of unburnt mixture [K]
%   p_in        : pressure of unburnt mixture [Pa] (if p_in==1 then 1 atm is assumed)
%
% Outputs:
%   s_l_u       : Laminar flame speed [m/s]
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.04.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% If p_in=1 then assume pressure in atm
if p_in == 1
  p_in = 101325;
end

% Compute s_L
if strcmpi(mixture,'MethaneAir')
  % Methan Air
  
  % Modelling flame speed fluctuations 
  %   according to Lieuwen 2003, Modeling Premixed Combustion–Acoustic Wave Interactions: A ReviewEq. (42))
  A = 0.6079; B = -2.554; C = 7.31; D = 1.23;
  s_l_u = A*Phi.^B .* exp(-C*(Phi-D).^2); % Mean flame speed at T_in=300K and p=101325Pa
  
  % Influence of Temperature and pressure 
  %   according to Poinsot, Veynante: Theoretical and Numerical Combustion, 2005, p.55
  % Interpolate coefficients in table
  alpha_T = interp1( [0.8 1 1.2] , [2.105 1.612 2.0] , Phi , 'linear','extrap');
  alpha_p = interp1( [0.8 1 1.2] , [-0.504 -0.374 -0.438] , Phi , 'linear','extrap');
  % Now compute flame speed
  s_l_u = s_l_u .* ( T_in / 300 ).^alpha_T .* ( p_in / 101325 ).^alpha_p;

elseif strcmpi(mixture,'myMixture')
  % Insert computation of s_l_u for you mixture here
  s_l_u = nan;
  
else
  error('Flame speed calculation for chosen mixture not implemented!')
end

end

