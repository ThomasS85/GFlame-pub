function [ p ] = settings_Example( )
%settingsFTF Sets parameters for computation FTF

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% General Settings
% Flame geometry: Vinv (inverse V) or V
p.geom = 'Vinv';
% Flame type: 'slit' or 'conical'
p.flameType = 'conical';
% Simulate full geometry or use symmetry? ( 'full' or 'sym' )
p.domainType = 'sym';
% name for output folder for output/ input data ('.' for case root dir)
p.caseFolder = '.';


%% Unburnt Mixture
% Mixture: MethaneAir
%     if another mixture is desired, define it in all functions stored at functions/Physics, e.g.
%      getLaminarFlameSpeed() -> see example mixture 'myMixture'!
p.mixture = 'MethaneAir';
% Mean equivallence ratio
p.Phi_m = 0.8;
% Temperature of unburnt mixture
p.T_in = 300;  % [K]
% Pressure of unburnt mixture
p.p_in = 101325;  % [Pa]
% kinematic viscosity
p.nu = 15.82e-6;      % [m^2/s]

% Laminar flame speed (Methane at Phi~0.86 -> See Diss. Cuquel p. 124 )
% if s_l_u is set to a value <= 0 it will be computet from mean equivallence ratio Phi_m
p.s_l_u = -1;  % [m/s]

% Axial flow speed at inlet
p.u_1_L = 1.12;   % [m/s]


%% Other Parameters Physics 
% Markstein length (curvature correction: activate curvature in solver
% settings!). If set to a value <=0 it will be estimated automatically!
p.marksteinLength = -1;
% burnt to unburnt gas volumetric expansion ratio: E = rho_u / rho_b
p.E = 6.5;


%% Parameters Geometry
% half diameter feed channel (=width of flame region!)
p.R_i = 11e-3;  % Cuquel
% p.R_i = 1.0e-3;   % Kornilov
% half diameter combustion chamber (not needed yet!)
p.R_a = 25e-3;    % Cuquel
% p.R_a = 2.5e-3;   % Kornilov
% For V flames: Radius of the central rod
p.R_rod = 0.5e-3;
% Flame lift off (Flame is anchored at x1 = x1Lim(1) + liftOff )
p.liftOff = 0;
% Flame lateral offset (Flame is anchored at x2 =  R_i + lateralOffset )
p.lateralOffset = 0;


%% Sketch of the geometry:
% Vinv-flame
%     ^
% x_2 |
%     |
% R_a |_    _____________________
%     |    |
%     |    |_   _   _   _   _   _   _
%     |    |    o                    ^
% R_i |____|_   _  o_   _   _   _   _v lateralOffset
%     |               o
%     |                  o
%     |                     o
%   0 |----|-------------------o--------> 
%          0                      x_1
%          |< >|
%         liftOff

% V-flame
%         ^
%     x_2 |
%         |
%     R_a |_    __________________________
%         |    |                      o
%         |    |_   _   _   _   _  o_   _
%         |    |                o            
%     R_i |____|_   _   _  _ o _   _   _   _
%         |              o     
%         |           o          
%  R_rod _|____    o                    
%       0 |----|-----------------------------> 
%              0                   x_1
%              |< >|
%             liftOff

% M-flame
%         ^
%     x_2 |
%         |
%     R_a |_    __________________________
%         |    |                      
%         |    |_   _   _   _   _  _   _
%         |    |                            
%     R_i |____|_   _   _  _  _   _   _   _
%         |        o           
%         |           o          
%  R_rod _|____    o                    
%       0 |----|-----------------------------> 
%              0                   x_1
%              |< >|
%             liftOff

end

