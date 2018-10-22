function [ p ] = setROMParameters( Fdim , type , R_i , R_a , s_l_u , u_1_bulkFeed , varargin )
%SETROMPARAMETERS Calculates and sets all necessary parameters for the ROM
%   for a given flame and checks user input 
%
% A flame is defined by:
%     (1) Dimension:      2D or 3D (axissymmetric)
%     (2) Shape:          inverseV or V (not tested!)
%     (3) Flame base radius/ width R_i (V flames: R_i is rod radius)
%     (4) Outer radius of combustion chamber R_a
%     (5) Laminar flame speed s_l_u (defined @ unburnt side)
%     (6) Bulk flow velocity in feed channel u_1_bulkFeed (important for mass flux)
% 
% Optional Input:
%     - name        : A name for the flame can be supplied 
%                     (default = type_Fdim)
%     - E           : Expansion ratio for steady state correction 
%                     E = rho_u / rho_b  (default = 1, no expansion)
%     - u1CenterIn  : Velocity at center axis at inlet (important for no-slip walls)
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Set combustor type (required to find appropriate SC mapping + to parse input)
if R_i<R_a
  p.CombType = 'backwardFacingStep';
elseif R_i==R_a
  p.CombType = 'duct';
else
  p.CombType = 'unknown';
end


%% Parse input
% Dimension
if ~ ( strcmp(Fdim,'3D') || strcmp(Fdim,'2D') )
  error('Unknown flame dimension chosen! (Use 2D or 3D!)')
end
% Shape
if ~ (strcmp(type,'inverseV') || strcmp(type,'V') || strcmp(type,'flat') )
  error('Unknown flame shape chosen! (Use inverseV or V!)')
elseif strcmp(type,'V')
  warning('V-shaped flames are not valiadated yet! Be carefull!)')
end
% radii
if R_i <= 0 || R_a <= 0
  error('Please choose a radius greater than zero!')
end

if R_i>R_a
  if strcmp(type,'V')
    error('Rod radius cannot be greater than combustion chamber radius!')
  else
    error('Combustion chamber radius cannot be greater than flame base radius!')
  end
end
% speeds
if u_1_bulkFeed / s_l_u > 10  
  % Flame angle would be smaller than 5.7 degrees!
  error('Flame angle too small. Such flames are not supported!')
elseif s_l_u >= u_1_bulkFeed && ~strcmpi(p.CombType,'duct')
  error('Flame speed has to be smaller than bulk flow velocity in feed channel!')
end


%% Parse varargin
% Set name
ind = find(strcmpi(varargin,'name'),1);
if ~isempty(ind)
  % Take flame name from varargin
  p.FlameName = varargin{ind+1};
else
  % Default flame name
   p.FlameName = [type,'_',Fdim];
end

% Set expansion ratio
ind = find(strcmpi(varargin,'E'),1);
if ~isempty(ind)
  % Take expansion ratio from varargin
  E = varargin{ind+1};
  if E<1 || E>10
    error('Check provided expansion ratio! Value seems to be unphysical.')
  end
else
  % Default expansion ratio
   E = 1;
end

% Set inlet velocity at flame foot (on the center axis)
ind = find(strcmpi(varargin,'u1CenterIn'),1);
if ~isempty(ind)
  % Take inlet velocity at center axis from varargin
  u_1_centerIn = varargin{ind+1};
  if u_1_centerIn<u_1_bulkFeed || u_1_centerIn>2*u_1_bulkFeed
    error('Check provided inlet velocity at center axis! Value seems to be unphysical.')
  end
else
  % Default: Constant inlet velocity profile
   u_1_centerIn = u_1_bulkFeed;
end

%% Set fixed parameters
% Flame dimension
p.Fdim = Fdim;
% Flame type
p.flameType = type;

% Combustion Chamber Radius
p.R_a = R_a;
% Flame speed
p.s_l_u = s_l_u;
% Expansion ratio
p.E = E;
% Flame Radius / rod radius
if strcmp(type,'V')
  p.R_r = R_i;
else
  p.R_i = R_i;
end


%% Empiric flame parameters 
% Lift off height 
p.H_liftOff = 0;
% Flame maximum radius (the close the flame is anchored at edge the worse
% is irrot part of IR)
if strcmp(type,'V')
  p.R_flame = p.R_a - p.R_r;
else
  p.R_flame = p.R_i;
end


%% Set derived parameters
% Combustion Chamber: Confinement Ratio
if strcmp(type,'V')
  p.r_tilde = p.R_r / p.R_a;
  p.Cr = 1 - p.r_tilde;
else
  p.Cr = p.R_i / p.R_a;
end

% Flame shape
% (1) Flame shape based on flow speed and flame speed (no confinement)
% Flame angle  
p.alpha0 = asin( p.s_l_u / u_1_centerIn );
p.alpha0Degree = p.alpha0 / pi * 180;
% Flame height
if ~strcmp(type,'flat')
  p.H_flame0 = p.R_flame / ( tan(p.alpha0) );
else
  p.H_flame0 = 0;
end
% Flame length
p.L_flame0 = p.R_flame / ( sin(p.alpha0) );
% Bulk flow velocity
p.u_1_bulkFeed = u_1_bulkFeed;
% Velocity at inlet (center axis)
p.u_1_centerIn0 = u_1_centerIn;

% (2) Correction due to confinement 
% Non-dimensional acceleration as function of confinement number Cn
p.acc_fun = @(Cn) ( 1 - 2*Cn*cos(p.alpha0)^2*(1-p.E) - 2*p.E*Cn.*(1-Cn)  ).^0.5 - 1;
% Critical confinement number (from Remie et al. 2006 and simplified with WolframAlpha)
p.Cn_crit = 1 - cos(p.alpha0) * (p.E-1)/p.E;

if strcmp(p.Fdim , '3D') && strcmp(p.flameType , 'inverseV')
  % Confinement number 
  p.Cn = p.Cr^2;
  if p.Cn>=p.Cn_crit
    % Acceleration 
    acc = p.acc_fun(p.Cn);
    % Flame Shape
    p.flameShape = @(r) 1/acc * ( ( 1-(1-r/p.R_flame)/(1+2/acc) ).^(-2) - 1 ) * p.H_flame0;
    % Flame Height: eq. 6.23
    p.H_flame = p.flameShape(0);
    % Flame angle
    p.alpha = atan( p.R_flame / p.H_flame );
    p.alphaDegree = p.alpha / pi * 180;
    % Flame length
    p.L_flame = p.R_flame / ( sin(p.alpha) );
    
    % if set to 1: inlet velocity is taken, 0: velocity at flame tip
    %   (2/3 since center of gravity of cone is 1/3 H_f away from base, also see Diss Cuquel p. 142ff)
    % fac_u0 = 2/3;
    % p.u_1_bulkFeed = ( p.u_1_bulkFeed * ( 1 + acc*p.H_flame/p.H_flame0 )*(1-fac_u0) + p.u_1_bulkFeed*fac_u0 ) ;
    p.u_1_centerIn = p.u_1_centerIn0 * acc / ( log( 1 + acc ) );
  else
    p.alpha = p.alpha0;
    p.alphaDegree = p.alpha0Degree;
    p.H_flame = p.H_flame0;
    p.L_flame = p.L_flame0;
    p.u_1_centerIn = p.u_1_centerIn0;
    % Flame shape (straight line)
    p.flameShape = @(r) p.H_flame - cot(p.alpha)*r;
  end
  
elseif strcmp(p.Fdim , '2D') && strcmp(p.flameType , 'inverseV')
  % Confinement number 
  p.Cn = p.Cr;
  if p.Cn>=p.Cn_crit
    % Acceleration 
    acc = p.acc_fun(p.Cn);
    % Flame Shape
    p.flameShape = @(r) ( ( 1 - r/p.R_flame ) ./ ( acc*r/p.R_flame + 1 ) ) * p.H_flame0;
    % Flame Height: eq. 6.23
    p.H_flame = p.flameShape(0);
    % Flame angle
    p.alpha = atan( p.R_flame / p.H_flame );
    p.alphaDegree = p.alpha / pi * 180;
    % Flame length
    p.L_flame = p.R_flame / ( sin(p.alpha) );
    
    % if set to 1: inlet velocity is taken, 0: velocity at flame tip
    %  (1/2 since center of gravity of cone is 1/2 H_f away from base)
    %fac_u0 = 1/2;
    %p.u_1_bulkFeed = ( p.u_1_bulkFeed0 * ( 1 + acc*p.H_flame/p.H_flame0 )*(1-fac_u0) + p.u_1_bulkFeed0*fac_u0 ) ;
    
    p.u_1_centerIn = p.u_1_centerIn0 * acc / ( log( 1 + acc ) );
    
    %myTau_r = p.R_flame / (acc^2*p.u_1_bulkFeed0) * ( acc - log(acc+1) );
    %p.u_1_bulkFeed = p.L_flame / (cos(p.alpha)*myTau_r);
  else
    p.alpha = p.alpha0;
    p.alphaDegree = p.alpha0Degree;
    p.H_flame = p.H_flame0;
    p.L_flame = p.L_flame0;
    p.u_1_centerIn = p.u_1_centerIn0;
    % Flame shape (straight line)
    p.flameShape = @(r) p.H_flame - cot(p.alpha)*r;
  end

elseif strcmp(p.Fdim , '3D') && strcmp(p.flameType , 'V')
  % Confinement number 
  p.Cn = 1 - p.r_tilde^2;
  if p.Cn>=p.Cn_crit
    % Acceleration 
    acc = p.acc_fun(p.Cn);
    % Radial velocity
    % u_hat = p.u_1_centerIn0 * ( (p.R_a^2-p.R_r^2)^2/(2*p.R_a^2) - (1+acc)*((p.R_a-p.R_r)^2/3 + (p.R_a-p.R_r)/2*p.R_r)) /...
    %   ( cot(p.alpha0)*( (p.R_a-p.R_r)^2/(3*sin(p.alpha0)) + p.R_r*(p.R_a-p.R_r)/(2*sin(p.alpha0)) - p.R_a*(p.R_a-p.R_r)/2-p.R_r*p.R_a ) );
    % Flame Shape Wedge
    %p.flameShape = @(r) 1/acc * ( exp( acc/cos(p.alpha0)/p.R_flame*(r-p.R_r) ) - 1 ) * p.H_flame0;
    p.flameShape = @(r) cot(p.alpha0) * ( r - p.R_r );
    % Flame Height: eq. 6.23
    p.H_flame = p.flameShape( p.R_a );
    % Flame angle
    p.alpha = atan( p.R_flame / p.H_flame );
    p.alphaDegree = p.alpha / pi * 180;
    % Flame length
    p.L_flame = p.R_flame / ( sin(p.alpha) );
    
    % if set to 1: inlet velocity is taken, 0: velocity at flame tip
    %   (2/3 since center of gravity of cone is 1/3 H_f away from base, also see Diss Cuquel p. 142ff)
    fac_u0 = 2/3 - p.R_r^2 / ( 2 * ( p.R_r + p.R_a )^2 );
    %fac_u0 = 0.56;
    p.u_1_centerIn = ( p.u_1_centerIn0 * ( 1 + acc*p.H_flame/p.H_flame0 )*(1-fac_u0) + p.u_1_centerIn0*fac_u0 ) ;
  else
    p.alpha = p.alpha0;
    p.alphaDegree = p.alpha0Degree;
    p.H_flame = p.H_flame0;
    p.L_flame = p.L_flame0;
    p.u_1_centerIn = p.u_1_centerIn0;
  end
  
else
  % Don't use correction due to confinement
  warning('No correction due to confinement performed since there exists no model for the chosen flame yet.')
  p.alpha = p.alpha0;
  p.alphaDegree = p.alpha0Degree;
  p.H_flame = p.H_flame0;
  p.L_flame = p.L_flame0;
  p.u_1_centerIn = p.u_1_centerIn0;

end


