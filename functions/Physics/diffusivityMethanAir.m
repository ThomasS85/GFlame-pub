function [ D_CH4_air ] = diffusivityMethanAir( T )
%DIFFUSIVITYMETHANAIR Calculates Diffusion coeffiocient for methane in air depending on Temperature
%
% From https://chrisbharding.wordpress.com/tag/chapman-and-enskog/   (search for methane)
% And http://www.chemie.de/lexikon/Diffusionskoeffizient.html
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.04.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% molar masses
% Air
M_air = 29; % g/mole
M_CH4 = 16; % g/mole
M_AB = 2 / ( 1/M_air + 1/M_CH4 ); % g/mole


%% Sigma/ eps_kappa
sigma_AB = 3.735; % A
eps_kappa_AB = 108.1; % K
T_s = T / eps_kappa_AB;

%% Collusion integral
A = 1.06036;
B = 0.15610;
C = 0.193;
D = 0.47635;
E = 1.03587;
F = 1.52996;
G = 1.76474;
H = 3.89411;

Omega_D = A ./ T_s.^B + C ./ exp(D*T_s) + E ./ exp(F*T_s) + G ./ exp(H*T_s);

% Diffusivity Methan in air (ideal gas law)
D_CH4_air = 0.00266 * T.^(3/2) ./ ( M_AB^(1/2)*sigma_AB^2*Omega_D ) * 1e-4; % [m^2/s]

end

