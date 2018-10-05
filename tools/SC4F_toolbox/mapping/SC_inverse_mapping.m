function [ xi_SC ] = SC_inverse_mapping( x , s , varargin )
%SC_MAPPING Maps a given complex vector with locations in the physical 
%       domain to the image domain using inverse Schwarz-Christoffel mapping
%
%   Mapping x -> xi
%
%   Input:  - x         : Complex vector with coordinates in physical
%                         domain (L2 system)
%           - s         : Solver settings struct
%   Output: - xi_SC     : Complex vector with coordinates in image
%                         domain xi = xi1 + i*xi2
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
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Initial value for root finding point
ind = find(strcmpi(varargin,'xi0'),1);
if ~isempty(ind)
  % Use user specified value
  xi0 = [ real(varargin{ind+1}) , imag(varargin{ind+1})  ];
else
  % use default value from p struct
  xi0 = [ real(s.xi_start) , imag(s.xi_start) ];
end

% Should as initial value for root finding inside the for loop the value in
% the step before be used (faster for mappings of continuous functions whose values are close
% next to each other -> NOT SUITABLE for points on CENTERLINE!)
ind = find(strcmpi(varargin,'updateInit'),1);
if ~isempty(ind)
  % Yes
  updateInit = 1;
else
  % No ,use always the same initial value
  updateInit = 0;
end



%% Define Parameters
% initialise output
xi_SC_vec = zeros(length(x),2);


%% Normalize lengths
x = x / s.l_ref;


%% Check if input is row or column vector
[r,c] = size(x);


%% Perform inverse SC-Mapping solving the complex implicit transformation equation
%Set options for fsolve and perform initial Step
fun = @(xi) SC_inverse_mapping_Residuum( xi , x(1) , s );
opts = optimoptions(@fsolve,'Display','off','TolFun',1e-10,'TolX',1e-10);
xi_SC_vec(1,:) = fsolve(fun,xi0,opts);

if updateInit
  % Initial value is flexible
  solveFun = @(fun,xi0n) fsolve(fun,xi0n,opts);
else
  % Always use same initial value (xi0)
  solveFun = @(fun,xi0n) fsolve(fun,xi0,opts);
end

% Loop over x-vector
for ii = 2:length(x) 
  % Solve inverse mapping problem (always using xi0 as initial guess)
  fun = @(xi) SC_inverse_mapping_Residuum( xi , x(ii) , s );
  xi_SC_vec(ii,:) = solveFun( fun , xi_SC_vec(ii-1,:) );
  
end

% Convert vector to imaginary number 
xi_SC = xi_SC_vec(:,1) + 1i*xi_SC_vec(:,2);


%% Output vector should be rotated like input
if c>r
  xi_SC = xi_SC.';  %% Attention: don't use only ' since this calcs conjugate transpose
end


end

