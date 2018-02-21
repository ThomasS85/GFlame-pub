function [ marksteinLength ] = getMarksteinLength( mixture , s_l_u , E , T_in )
%GETMARKSTEINLENGTH Function computes the Markstein length for the laminar flame speed (consumption rate) for
%a given mixture
%
% Inputs:
%   mixture     : Name of the mixture
%   s_l_u       : Laminar flame speed [m/s]
%   E           : Expansion ratio:  E = rho_u / rho_b [-]
%   T_in        : Temperature of unburnt mixture [K]
%
% Outputs:
%   marksteinLength  : Markstein length
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.04.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


if strcmpi(mixture,'MethaneAir')
  % Methan Air
  
  % Quantities for the unburnt gas
  lambda_u = 0.0272;                      % [W/m/K] : From cantera (Phi=1,T_in=300K)
  cp_u = 1080;                            % [J/kg/K] : From cantera (Phi=1,T_in=300K)
  rho_u = 1.12;                           % [kg/m^3] : From cantera (Phi=1,T_in=300K)
  a = lambda_u / (rho_u*cp_u);            % Thermal diffusivity unburnt mixture
  D_CH4_air = diffusivityMethanAir( T_in );  % Diffusion coeffiocient for methane in air
  
  % Lewis Number Le
  Le = a / D_CH4_air;
  
  % Flame thickness
  delta_f = a / s_l_u;
  
  
  % Markstein number 
  %   according to Poinsot, Veynante: Theoretical and Numerical Combustion, 2005
  % alpha = 1 - 1 / E;
  beta = 18.4 ;           % according to p.41, table 2.2 (Flame 2)
  fac = E - 1;            % fac = (T_b-T_u)/T_u
  Ma = 0.5 * beta *( Le - 1 ) / fac * ( polylog(2,0) - polylog(2,-fac) );   % p.71, Eq.(2.110)
  
  
  % According to: Tseng et. al., 1993. Laminar Burning Velocities and Markstein Numbers of 
  %   Hydrocarbon/Air Flames. Combustion and Flame 95, 410–426.
  % Ma = 10.2 * ( phi - 0.74 ); % yields very similar results!  CHECK!

  
  % Compute Markstein length
  marksteinLength = Ma * delta_f;
  
  
elseif strcmpi(mixture,'myMixture')
  % Insert computation of marksteinLength for you mixture here  
  marksteinLength = nan;
  
else
  error('Markstein calculation for chosen mixture not implemented!')
end


end



