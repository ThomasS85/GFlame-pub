function [ M_sum , rho_sum , x , y ] = getMoleMassFromPhi( mixture , Phi )
%GETMOLEFRACTIONFROMPHI Calculates the mole fractions x_i=n_i/n_sum or the mass fractions y_i = m_i/m_sum
% for a given equivalence ratio Phi for methan mixed with air which consists of 21% O_2 and 79% N
%
% Inputs:
%   mixture : Name of the mixture
%   Phi     : Desired equivalence ratio
%
% Outputs:
%   M_sum   : Molare mass of the mixture [g/mol]
%   rho_sum : Density of the mixture [kg/m^3]
%   x       : mole fractions
%   y       : mass fractions
%
% Global reaction (only with oxygen): CH4 + 2 * O2 -> C02 + 2 * H20
%
% by Thomas Steinbacher (copyed from TFDtools: getMoleMassFractionsFromPhi() )



%% Settings
% Reference values (for determining density)
p_ref = 101325; %[Pa]
T_ref = 293; % [K]
% Universal gas constant
R_m = 8.3144598; %[J/mol/K]
% Precision
prec1 = 8;

%% Calculate mass and mole fractions
% Definition air
x_O2_air = 0.21;
x_N2_air = 0.79;


if strcmpi(mixture,'MethaneAir')
  % Molar masses [kg/kmol]
  M_O2 = 31.9988;
  M_N2 = 14.0067*2;
  M_CH4 = 16.043;
  
  
  % Calculate m_fuel2air for stoichiometric conditions
  m_fuel2air_st = M_CH4 / ( 2*M_O2 + 2*x_N2_air/x_O2_air* M_N2 );
  
  % Calculate m_fuel2air
  m_fuel2air = Phi * m_fuel2air_st;
  
  % Calculate x_CH4
  x_CH4 = m_fuel2air * ( M_N2*x_N2_air + M_O2*x_O2_air ) ./ ( M_CH4*(x_N2_air+x_O2_air) +...
    m_fuel2air * ( M_N2*x_N2_air + M_O2*x_O2_air ) );
  
  % Calculate x_O2
  x_O2 = x_O2_air * ( 1 - x_CH4 ) ./ ( x_O2_air + x_N2_air );
  
  % Calculate x_N
  x_N2 = x_N2_air * x_O2 ./ x_O2_air;
  
  % Calculate Mass fractions from mole fractions
  M_sum =  M_CH4*x_CH4 + M_O2*x_O2 + M_N2*x_N2 ;
  y_CH4 = M_CH4 ./ M_sum .* x_CH4;
  y_O2 = M_O2 ./ M_sum .* x_O2;
  y_N2 =  M_N2 ./ M_sum .* x_N2;
  
  
  
  % Compute density
  rho_sum = p_ref * M_sum*1e-3 / ( R_m*T_ref );
  
  % Compute cv
  
else
  error('Flame speed calculation for chosen mixture not implemented!')
end

% Output
% mass fractions
y.CH4 = y_CH4;
y.N2 = y_N2;
y.O2 = y_O2;

% mole fractions
x.CH4 = x_CH4;
x.N2 = x_N2;
x.O2 = x_O2;

% Only display if no output is requested
if nargout==0
  disp(['Phi: ',num2str(Phi,prec1)])
  disp(['y_CH4: ',num2str(y_CH4,prec1)])
  disp(['y_O2: ',num2str(y_O2,prec1)])
  disp(['y_N2: ',num2str(y_N2,prec1)])
  disp(['x_CH4: ',num2str(x_CH4,prec1)])
  disp(['x_O2: ',num2str(x_O2,prec1)])
  disp(['x_N2: ',num2str(x_N2,prec1)])
  disp(['Resulting Molare Mass: ',num2str(M_sum,prec1),' [g/mol]'])
  disp(['Resulting Density @ Reference Conditions (Ideal Gas): ',num2str(rho_sum,prec1),' [kg/m^3]'])
end

%% Debug
% if length(Phi) == 1
%   % Exact version
%   Phi_test = M_CH4 * x_CH4 / ( M_O2*x_O2 + M_N2*x_N2 )  / m_fuel2air_st;
%   disp(['Sum (=1): ',num2str(x_CH4 + x_O2 + x_N2)]);
%   disp(['Phi test: ',num2str(Phi_test,6)])
%   % Version alp
%   Phi_test_alp = M_CH4 * x_CH4_alp / ( M_O2*x_O2_alp + M_N2*x_N2_alp )  / m_fuel2air_st;
%   disp(['Sum Alp (=1): ',num2str(x_CH4_alp + x_O2_alp + x_N2_alp)]);
%   disp(['Phi test Alp: ',num2str(Phi_test_alp,6)])
% end


end

