function p = setDerivedParamters(p)
% Calculates all parameters from user settings
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Derived parameters physics/ geometry
% Flame speed: If flame speed was set to a value <=0 then calculate flame speed from Phi_m
if p.s_l_u <= 0
  [ p.s_l_u ] = getLaminarFlameSpeed( p.mixture , p.Phi_m , p.T_in , p.p_in );
end

% Markstein Number: If markstein length was set to a value <=0 then calculate it
if p.marksteinLength <=0
  [ p.marksteinLength ] = getMarksteinLength( p.mixture , p.s_l_u , p.E , p.T_in );
end

% Flame
% Flame angle
p.alpha = asin( p.s_l_u / p.u_1_L );
p.alphaDegree = p.alpha / pi * 180;
% Flame Radius
if strcmpi(p.geom,'Vinv') ||  strcmpi(p.geom,'M')
  p.R_flame = p.R_i + p.lateralOffset;
elseif strcmpi(p.geom,'V')
  p.R_flame = p.R_a;
end
% Flame height
p.H_flame = p.R_flame / ( tan(p.alpha) );
% Flame length
p.L_flame = p.R_flame / ( sin(p.alpha) );
% Density expansion raio
p.gamma = 1 - 1/p.E;
% Convective time scale
p.tau_c = p.L_flame / ( p.u_1_L * cos(p.alpha) );

% Physical quantities
% Density and molar mass
[ p.M_sum , p.rho_in ] = getMoleMassFromPhi( p.mixture , p.Phi_m );

% Confinement
% Confinement Ratio Cr
p.Cr = p.R_i / p.R_a;
% Gamma* for velocity model with confinement
if strcmpi(p.geom,'Vinv')
  if strcmp(p.flameType,'conical')
    % (See Diss Cuquel p. 133)
    p.gamma_s = sqrt( 1 - 2*p.Cr^2 *( 1 + sin(p.alpha)^2 * (p.E-1) -p.E*p.Cr^2 ) ) - 1;
  elseif strcmp(p.flameType,'slit')
    % See derivation MUPAD
    p.gamma_s = sqrt( 1 - 2*p.Cr*cos(p.alpha)^2 - 2*p.Cr^2*p.E*( sin(p.alpha)^2/p.Cr - 1 ) ) - 1;
  end
else
  % Nothing derived yet -> Ignore influence!
  p.gamma_s = 0;
end


end