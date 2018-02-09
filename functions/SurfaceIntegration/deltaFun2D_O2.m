function [ delta_sum, stencil_2nd , dataOut ] = deltaFun2D_O2( grid , data  , s )
% deltaFun2D_O1: calculate delta function with second order accuracy 
% according to Smereka 2006
%
% The delta function is only evaluated in a small band around the zero
% level iso line. This is archieved by setting all values of data outside
% the band to zero and subsequentely creating a sparse.
%
% This version doesn't include Boundary Conditions!
%
% DO NOT USE THIS VERSION! DEBUGGING SHOWS IT DOES NOT WORK PROPERLY!
%  Seems to have problems when the resolution is not so high
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Method only applicable if spacing in x1 and x2 direction is equal
% In order to fix this one has to treat the regular nodes differently
% (=which are axactely 0).
if grid.dx(1) ~= grid.dx(2)
  error('Method only applicible for equal x1 and x2 spacing!')
end

% Cut cells from every dimension according to stencil.
stencil_2nd = 2;    % 2nd order derivative!
stencil_1st = 1;    % central/ forward/ backward differences (1st order)

% define epsilon which prevents gradient for becoming 0 and thus devision
% by 0
epsilon = 1e-13;

% This factor defines how many cells should be linearily interpolated at maximum. Important in the case of a
% slope close to zero!
facHmax = 1;


%% (0) Collect some data
% Number of grid points in each dimension (without BC)
N_x1 = grid.N(1) ;
N_x2 = grid.N(2) ;
% Spacing in each dimension
h_x1 = grid.dx(1);
h_x2 = grid.dx(2);


%% (1) Prepare data fields
% For reinitialisation tempoaly change BC to extrapolation (no influence of BC!)
grid.bdry{1,1} = @addGhostExtrapolate;
grid.bdry{2,1} = @addGhostExtrapolate;
grid.bdryData = cell(2,1);
% Reinitialize data (output this field)
dataOut = signedDistanceIterativeSubCellFix(grid, data, s.deltaFunReinit.accuracy ,...
  s.deltaFunReinit.tMax, s.deltaFunReinit.errorMax);

  
% Set every entry whose absolute value is smaller eps to zero (otherwise
% this could lead to nan)
data = dataOut;
nodesTrue = sparse( abs(data)<=eps );
data( nodesTrue ) = 0;

%  Find all Regular nodes (=0); important to do that before cutting out
%  band!
deltaReg = sparse( grid.N(1) , grid.N(2) );
deltaReg(nodesTrue) = 1 /  h_x1;

% Only perform calculation in a small band around zero line
bandWindow = stencil_2nd + 1;
band = ( identifyInterfaceBand( data , bandWindow ) );
data(~band) = 0;
data = sparse(data);

% Cut cells from every dimension according to stencil.
dataCut = data( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );
deltaReg = deltaReg( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );



%% (2) 1st order terms
%Calculation of derivatives (as matrix-matrix multiplication with sparse!)

%%%%%%%%%%%%%%%%%%%%%
% (2.1) forward difference 
%%%%%%%%%%%%%%%%%%%%%
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[0 1],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[0 1],N_x2,N_x2);
Dx1_fwd = Dx1*data / h_x1;
Dx2_fwd = data*Dx2 / h_x2;
Dx1_fwd = Dx1_fwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
Dx2_fwd = Dx2_fwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );

%%%%%%%%%%%%%%%%%%%%%
% (2.2) backward difference
%%%%%%%%%%%%%%%%%%%%%
% Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[-1 0],N_x1,N_x1);
% Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[-1 0],N_x2,N_x2);
% Dx1_bwd = Dx1*data / h_x1;
% Dx2_bwd = data*Dx2 / h_x2;
% Dx1_bwd = Dx1_bwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
% Dx2_bwd = Dx2_bwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );

%%%%%%%%%%%%%%%%%%%%%
% (2.3) central difference
%%%%%%%%%%%%%%%%%%%%%
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[-1 1],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[-1 1],N_x2,N_x2);
Dx1_ced = Dx1*data / ( 2*h_x1 );
Dx2_ced = data*Dx2 / ( 2*h_x2 );
Dx1_ced = Dx1_ced( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
Dx2_ced = Dx2_ced( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );

%%%%%%%%%%%%%%%%%%%%%
% (2.4) Absolute value gradient
%%%%%%%%%%%%%%%%%%%%%
gradient_abs = sqrt( Dx1_ced.^2 + Dx2_ced.^2 + epsilon );

%%%%%%%%%%%%%%%%%%%%%
% (2.5) Normal vector
%%%%%%%%%%%%%%%%%%%%%
n_x1 = Dx1_ced ./ gradient_abs;
n_x2 = Dx2_ced ./ gradient_abs;

% All 1st order fields are now shorter by 2*stencil:


%% (3) 2nd order terms for jump in second derivative

%%%%%%%%%%%%%%%%%%%%%
% (3.1) Derivation of the normal vector (central difference)
% Also shorten derivative operator matrices
%%%%%%%%%%%%%%%%%%%%%
Dx1 = Dx1( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st ); % Dx1 still from centraf diff.
Dx2 = Dx2( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st ); % Dx2 still from centraf diff.
% nx1
Dnx1Dx1 = Dx1*n_x1 / ( 2*h_x1 );
Dnx1Dx1 = Dnx1Dx1( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
Dnx1Dx2 = n_x1*Dx2 / ( 2*h_x2 );
Dnx1Dx2 = Dnx1Dx2( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
% nx2
Dnx2Dx1 = Dx1*n_x2 / ( 2*h_x1 );
Dnx2Dx1 = Dnx2Dx1( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
Dnx2Dx2 = n_x2*Dx2 / ( 2*h_x2 );
Dnx2Dx2 = Dnx2Dx2( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
% cut also n_x1/ nx_2
n_x1 = n_x1( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
n_x2 = n_x2( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );

%%%%%%%%%%%%%%%%%%%%%
% (3.2) Calculate D
%%%%%%%%%%%%%%%%%%%%%
D_mat = ( n_x1.^2 - n_x2.^2 ) .* ( n_x1 .* ( n_x2.*Dnx2Dx1 - n_x1.*Dnx2Dx2 ) + ...
                                   n_x2 .* ( n_x1.*Dnx1Dx2 - n_x2.*Dnx1Dx1 )    );

%%%%%%%%%%%%%%%%%%%%%
% (3.3) calculate second derivatives
%%%%%%%%%%%%%%%%%%%%%
D2x1x1gx1 = -D_mat .* sign(n_x1);
D2x2x2gx2 =  D_mat .* sign(n_x2);

% Again, fields D2xxgx/D2yygy are now shorter by 2*stencil (absolute: 4*stencil)


% Cut all first order terms one more time so that matrix sizes are consistent
Dx1_fwd = Dx1_fwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
Dx2_fwd = Dx2_fwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
% Dx1_bwd = Dx1_bwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );
% Dx2_bwd = Dx2_bwd( stencil_1st+1:end-stencil_1st , stencil_1st+1:end-stencil_1st );



%% (4) Approximate h_x1/x2^+ to 3rd order

%%%%%%%%%%%%%%%%%%%%%
% (4.1) h_x1^+
%%%%%%%%%%%%%%%%%%%%%
% Phi_c
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)*9 ones(N_x1,1)*9 -ones(N_x1,1)],[-1 0 1 2],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)*9 ones(N_x2,1)*9 -ones(N_x2,1)],[-1 0 1 2],N_x2,N_x2);
Phi_c_x1 = Dx1*data / 16;
Phi_c_x2 = data*Dx2 / 16;
Phi_c_x1 = Phi_c_x1( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );
Phi_c_x2 = Phi_c_x2( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );

% Phi_c'
Phi_c_prime_x1 = Dx1_fwd;
Phi_c_prime_x2 = Dx2_fwd;

% Phi_c''
Dx1 = spdiags([ones(N_x1,1) -ones(N_x1,1) -ones(N_x1,1) ones(N_x1,1)],[-1 0 1 2],N_x1,N_x1);
Dx2 = spdiags([ones(N_x2,1) -ones(N_x2,1) -ones(N_x2,1) ones(N_x2,1)],[-1 0 1 2],N_x2,N_x2);
Phi_c_prime2_x1 = Dx1*data / ( 2*h_x1^2 );
Phi_c_prime2_x2 = data*Dx2 / ( 2*h_x2^2 );
Phi_c_prime2_x1 = Phi_c_prime2_x1( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );
Phi_c_prime2_x2 = Phi_c_prime2_x2( stencil_2nd+1:end-stencil_2nd , stencil_2nd+1:end-stencil_2nd );


u_I_x1 = - Phi_c_x1 ./ Phi_c_prime_x1 - Phi_c_prime2_x1.*Phi_c_x1.^2 ./ ( 2*Phi_c_prime_x1.^3 );
% Remove NaN which come since at some points there is 0/0 ( reduces number of elements in sparse matrix)
u_I_x1(isnan(u_I_x1)) = 0;
% Now compute h_x1Plus (absolute value)
h_x1Plus = sparse( abs( h_x1 / 2 - u_I_x1 ) );
h_x1Plus(u_I_x1==0) = 0; % remove non-zero elements outside the band
% First order approx if interface not well resolved 
node_data = data(  stencil_2nd+2:end-stencil_2nd+1 , stencil_2nd+1:end-stencil_2nd ); % Phi_(i+1,j)
nodes_true = abs(u_I_x1) >= h_x1/2;
h_x1Plus(nodes_true) = abs( node_data(nodes_true) ./ Dx1_fwd(nodes_true) ); 
h_x1Plus(h_x1Plus>=h_x1*facHmax) = 0; % important if slope is close to zero!


% Compute  h_x1^- from h_x1^+ (at boundaries do not extrapolate!)
h_x1Minus = sparse( [ zeros(1,N_x2-2*stencil_2nd) ; h_x1 - h_x1Plus(1:end-1,:) ] );
h_x1Minus(h_x1Minus==h_x1) = 0;
h_x1Plus = [ h_x1Plus(1:end-1,:) ; zeros(1,N_x2-2*stencil_2nd) ]; % also set boundary for h^+ to zero: do not extrapolate!

% % Compute  h_x1^- (absolute value)
% h_x1Minus = abs( - h_x1 / 2 + u_I_x1 );
% h_x1Minus(u_I_x1==0) = 0; % remove non-zero elements outside the band
% % First order approx if interface not well resolved 
% node_data = data( stencil_2nd:end-stencil_2nd-1 , stencil_2nd+1:end-stencil_2nd ); % Phi_(i-1,j)
% h_x1Minus(nodes_true) = abs( node_data(nodes_true) ./ Dx1_bwd(nodes_true) ); 
% h_x1Minus(h_x1Minus>=h_x1*facHmax) = 0; % important if slope is close to zero!


u_I_x2 = - Phi_c_x2 ./ Phi_c_prime_x2 - Phi_c_prime2_x2.*Phi_c_x2.^2 ./ ( 2*Phi_c_prime_x2.^3 );
% Remove NaN which come since at some points there is 0/0 ( reduces number of elements in sparse matrix)
u_I_x2(isnan(u_I_x2)) = 0;
% Now compute h_x2Plus
h_x2Plus = sparse( h_x2 / 2 - u_I_x2 );
h_x2Plus(u_I_x2==0) = 0; % remove non-zero elements far away from band
% First order approx if interface not well resolved 
node_data = data(  stencil_2nd+1:end-stencil_2nd , stencil_2nd+2:end-stencil_2nd+1 ); % Phi_(i,j+1)
nodes_true = abs(u_I_x2) >= h_x2/2;
h_x2Plus(nodes_true) = abs( node_data(nodes_true) ./ Dx2_fwd(nodes_true) ); 
h_x2Plus(h_x2Plus>=h_x2*facHmax) = 0; % important if slope is close to zero!

% Compute  h_x1^- from h_x1^+ (at boundaries do not extrapolate!)
h_x2Minus = sparse( [ zeros(N_x1-2*stencil_2nd,1) , h_x2 - h_x2Plus(:,1:end-1) ] );
h_x2Minus(h_x2Minus==h_x1) = 0;
h_x2Plus = [ h_x2Plus(:,1:end-1) , zeros(N_x1-2*stencil_2nd,1) ]; % do not extrapolate!

% % Compute  h_x2^-  (absolute value)
% h_x2Minus = abs( - h_x2 / 2 + u_I_x2 );
% h_x2Minus(u_I_x2==0) = 0; % remove non-zero elements outside the band
% % First order approx if interface not well resolved 
% node_data = data(  stencil_2nd+1:end-stencil_2nd , stencil_2nd:end-stencil_2nd-1); % Phi_(i,j-1)
% h_x2Minus(nodes_true) = abs( node_data(nodes_true) ./ Dx2_bwd(nodes_true) ); 
% h_x2Minus(h_x2Minus>=h_x2*facHmax) = 0; % important if slope is close to zero!


%% (4) Calculate the 4 components of the delta function 

% (4.0) Initialize delta functions as sparse
delta_x1_up = deltaReg * 0;
delta_x1_down = delta_x1_up;
delta_x2_up = delta_x1_up;
delta_x2_down = delta_x1_up;


% (4.1) x1 direction
% (4.1.1) delta x1_up (+x1)
% include BC for shift (but matrix should have same size as data!)
node_data = data( stencil_2nd+2:end-stencil_2nd+1 , stencil_2nd+1:end-stencil_2nd ); 
% don't allow data and node_data to be zero at the same time!
% nodes_true =  data .* node_data < 0 | ( data==0 & node_data~=0) ;
nodes_true = sparse( dataCut .* node_data < 0 );  % In paper this is <= !! Here the =0 case is included in deltaReg!
delta_x1_up(nodes_true)   = ( abs( h_x1Plus(nodes_true) .* n_x1(nodes_true) ) +...
  0.5 * h_x1Plus(nodes_true).^2 .* D2x1x1gx1(nodes_true) ) / h_x1^2;

% (4.1.2) delta x1_down (-x1)
node_data = data( stencil_2nd:end-stencil_2nd-1 , stencil_2nd+1:end-stencil_2nd );
nodes_true = sparse( dataCut .* node_data < 0 );
delta_x1_down(nodes_true) = ( abs( h_x1Minus(nodes_true) .* n_x1(nodes_true) ) -...
  0.5 * h_x1Minus(nodes_true).^2 .* D2x1x1gx1(nodes_true) ) / h_x1^2;

% (4.2) x2 direction
% (4.2.1) delta x2_up (+x2)
node_data = data( stencil_2nd+1:end-stencil_2nd , stencil_2nd+2:end-stencil_2nd+1 );
% don't allow data and node_data to be zero at the same time!
% nodes_true =  data .* node_data < 0 | ( data==0 & node_data~=0) ;
nodes_true = sparse( dataCut .* node_data < 0 );  % In paper this is <= !! Here the =0 case is included in deltaReg!
delta_x2_up(nodes_true) = ( abs( h_x2Plus(nodes_true) .* n_x2(nodes_true) ) +...
  0.5 * h_x2Plus(nodes_true).^2 .* D2x2x2gx2(nodes_true) ) / h_x2^2;

% (4.2.2) delta x2_down (-x2)
node_data = data( stencil_2nd+1:end-stencil_2nd , stencil_2nd:end-stencil_2nd-1 ); 
nodes_true = sparse( dataCut .* node_data < 0 );
delta_x2_down(nodes_true) = ( abs( h_x2Minus(nodes_true) .* n_x2(nodes_true) ) -...
  0.5 * h_x2Minus(nodes_true).^2 .* D2x2x2gx2(nodes_true) ) / h_x2^2;
                            
                            

%% Summation of all terms
delta_sum = delta_x1_up + delta_x1_down + delta_x2_up + delta_x2_down + deltaReg;

% At boundaries the forward or backward difference can become 0 and thus
% the delta function might show some inf entries. Replace them by 0:
delta_sum( delta_sum == inf ) = 0;

end

