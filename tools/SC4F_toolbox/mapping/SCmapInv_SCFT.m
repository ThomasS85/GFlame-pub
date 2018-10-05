function [ xi_SC ] = SCmapInv_SCFT( x , p , varargin )
%SC_MAPPING Maps a given complex vector with locations in the physical 
%       domain to the image domain using inverse Schwarz-Christoffel mapping
%
%   Mapping x -> xi
%
%   Input:  - x         : Complex vector with coordinates in physical
%                         domain (L2 system; if provided in L1, set option accordingly!)
%           - p         :  Struct with flame settings as returned from setUpPredefinedFlame()
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
% To do:
%   - remove updateInit as default
%   - compute optimal xi= for each given x
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 09.11.2017                                //
% // Last modified: 09.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Initial value for root finding point
ind = find(strcmpi(varargin,'xi0'),1);
if ~isempty(ind)
  % Use user specified value
  xi0 = [ real(varargin{ind+1}) , imag(varargin{ind+1})  ];
else
  % use default value from p struct 
  xi0( real(x)>=0 , : ) = repmat( [ 1.05 , 0 ] , length(x(real(x)>=0)) , 1 );
  xi0( real(x)<0  , : ) =  repmat( [ 0.95 , 0 ] , length(x(real(x)<0))  , 1 );
end
% Get number of initial values
nXi0 = length(xi0(:,1));

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

% Is x provided in L1 coordinate system? If so, first transform to L2!
ind = find(strcmpi(varargin,'L1'),1);
if ~isempty(ind)
  % Do map from L1 to L2
  mapL1L2 = 1;
else
  % Don not map to L2; x already in L2!
  mapL1L2 = 0;
end


%% Map from L1 to L2 system if desired
if mapL1L2
  x = L1_to_L2(x,p);
end


%% Initialization
% initialise output
xi_SC_vec = zeros(length(x),2);
% Get mapping infos
[ s ] = return_SCmap_SCFT( p );
% % Mapping does not work properly far away from area jump -> warning!
% if any(real(x)>3*(p.R_a-p.R_i)) || any(real(x)<-(p.R_a-p.R_i))
%   warning('Inverse mapping only works relyable  close to area jump. Desired point, however, is far off!')
% end


%% If x is infinity, check mapping if this is mapped to a finite point. Otherwise skip point with warning
indInf = find( isinf(x) );
% Find lowest non-zero index
ind1 = find( ~isinf(x) , 1 , 'first' );
% Use mapping information from s if available
if isfield(s,'vertexes') && isfield(s,'prevertexes')
  for ii=1:length(indInf)
    indMatch = real(x(indInf(ii))) == real(s.vertexes);
    if ~isempty(indMatch)
      xi_SC_vec(indInf(ii),1) = real( s.prevertexes(indMatch) );
      xi_SC_vec(indInf(ii),2) = imag( s.prevertexes(indMatch) );
    else
      xi_SC_vec(indInf(ii),:) = [nan ,nan];
      warning(['Could not find inverse mapping for ',num2str(x(indInf(ii))),'! NaN returned'])
    end
  end
end

% Set indInf to 0 if it is empty (avoids problems later in for loop)
if isempty(indInf); indInf = 0; end

% Return if no finite entry was found
if isempty(ind1)
  % Convert vector to imaginary number 
  xi_SC = xi_SC_vec(:,1) + 1i*xi_SC_vec(:,2);
  return
end

%% Normalize lengths
x = x / s.l_ref;


%% Check if input is row or column vector
[r,c] = size(x);


%% Perform inverse SC-Mapping solving the complex implicit transformation equation
%Set options for fsolve and solve for first x-vector entry
fun = @(xi) SC_myResiduum_SCFT( xi , x(ind1) , s );
opts = optimoptions(@fsolve,'Display','off','TolFun',1e-12,'TolX',1e-12);
xi_SC_vec(ind1,:) = fsolve(fun,xi0(ind1,:),opts);

if updateInit
  % Initial value is flexible
  solveFun = @(fun,xi0n,ii) fsolve(fun,xi0n,opts);
else
  if nXi0 == length(x)
    % Use specific initial value from xi0 vector
    solveFun = @(fun,xi0n,ii) fsolve(fun,xi0(ii,:),opts);
  else
    % Always use same initial value (xi0)
    solveFun = @(fun,xi0n,ii) fsolve(fun,xi0(ind1,:),opts);
  end
end

% Loop over rest of x-vector
for ii = ind1+1:length(x) 
  % Skip if mapped from +-infinity
  if ii~=indInf
    % Solve inverse mapping problem 
    fun = @(xi) SC_myResiduum_SCFT( xi , x(ii) , s );
    xi_SC_vec(ii,:) = solveFun( fun , xi_SC_vec(ii-1,:) , ii );
  end
end

% Convert vector to imaginary number 
xi_SC = xi_SC_vec(:,1) + 1i*xi_SC_vec(:,2);



%% Output vector should be rotated like input
if c>r
  xi_SC = xi_SC.';  %% Attention: don't use only ' since this computes conjugate transpose
end


end


%% Residuum function
function [ res ] = SC_myResiduum_SCFT( xi , x , s )
%SC_INVMAPRES returns the residuum of the implicit inverse SC mapping 
%     function for a given xi and x
%
%  Inputs:  - xi       : Guess for solution as vector [real , imag]
%           - x        : Right hand side of equation as complex number
%
%  Outputs: - res      : Residuum of complex inverse mapping equation 
%                         SCmapping(xi)-x
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////

% Convert xi vector to xi as complex number:
xi_compl = xi(1) + 1i*xi(2);
% real part
res(1) =  real(s.x_xi(xi_compl)/s.l_ref - x);
% Imaginary part
res(2) =  imag(s.x_xi(xi_compl)/s.l_ref - x);

end