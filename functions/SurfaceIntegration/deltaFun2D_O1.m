function [ delta_sum , stencil , dataOut ] = deltaFun2D_O1( grid , data , s  )
% deltaFun2D_O1: calculate delta function with forst order accuracy 
% according to Smereka 2006
% 
% The delta function is only evaluated in a small band around the zero
% level iso line. This is archieved by setting all values of data outside
% the band to zero and subsequentely creating a sparse.
%
% This version doesn't use the defined Boundary Conditions!
%
% Improvements: We need only to calculate h_x1Plus and from there we could calculate h_x1Minus. This way we do
% not need the backward differences. This would speed up calculations, however, we cannot use the nodeTrue
% thing since we have to deal with vectors!
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

% Define stencil (also defines how many cells are cut at the edges)
stencil = 1; % central/ forward/ backward differences

% define epsilon which prevents gradient for becoming 0 and thus devision
% by 0
epsilon = 1e-12;

% This factor defines how many cells should be linearily interpolated at maximum. Important in the case of a
% slope close to zero!
% facHmax = max( floor(s.deltaFunReinit.tMax / grid.dx(1)) - 1 , 1);
facHmax = 1;


% First order accuracy

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
% Reinitialize data
dataOut = signedDistanceIterativeSubCellFix(grid, data, s.deltaFunReinit.accuracy ,...
  s.deltaFunReinit.tMax, s.deltaFunReinit.errorMax );

  
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
bandWindow = stencil + 1;
band = ( identifyInterfaceBand( data , bandWindow ) );
data(~band) = 0;
data = sparse(data);

% Cut cells from every dimension according to stencil.
dataCut = data( stencil+1:end-stencil , stencil+1:end-stencil );
deltaReg = deltaReg( stencil+1:end-stencil , stencil+1:end-stencil );


%% (2) Calculation of derivatives (as matrix-matrix multiplication with sparse!)

% (2.1) forward difference 
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[0 1],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[0 1],N_x2,N_x2);
Dx1_fwd = Dx1*data / h_x1;
Dx2_fwd = data*Dx2 / h_x2;
Dx1_fwd = Dx1_fwd( stencil+1:end-stencil , stencil+1:end-stencil );
Dx2_fwd = Dx2_fwd( stencil+1:end-stencil , stencil+1:end-stencil );

% (2.2) backward difference
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[-1 0],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[-1 0],N_x2,N_x2);
Dx1_bwd = Dx1*data / h_x1;
Dx2_bwd = data*Dx2 / h_x2;
Dx1_bwd = Dx1_bwd( stencil+1:end-stencil , stencil+1:end-stencil );
Dx2_bwd = Dx2_bwd( stencil+1:end-stencil , stencil+1:end-stencil );

% (2.3) central difference
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[-1 1],N_x1,N_x1);
Dx2 = spdiags([-ones(N_x2,1) ones(N_x2,1)],[-1 1],N_x2,N_x2);
Dx1_ced = Dx1*data / ( 2*h_x1 );
Dx2_ced = data*Dx2 / ( 2*h_x2 );
Dx1_ced = Dx1_ced( stencil+1:end-stencil , stencil+1:end-stencil );
Dx2_ced = Dx2_ced( stencil+1:end-stencil , stencil+1:end-stencil );

% (2.4) Absolute value gradient
gradient_abs = sqrt( Dx1_ced.^2 + Dx2_ced.^2 + epsilon );

% (2.5) Normal vector
n_x1 = Dx1_ced ./ gradient_abs;
n_x2 = Dx2_ced ./ gradient_abs;

% All 1st order fields are now shorter by 2*stencil:


%% (3) Calculate the 4 components of the delta function 

% (3.0) Initialize delta functions as sparse
delta_x1_up = deltaReg * 0;
delta_x1_down = delta_x1_up;
delta_x2_up = delta_x1_up;
delta_x2_down = delta_x1_up;


% (3.1) x1 direction
% (3.1.1) delta x1_up (+x1)
% include BC for shift (but matrix should have same size as data!)
node_data = data( stencil+2:end-stencil+1 , stencil+1:end-stencil ); 
% don't allow data and node_data to be zero at the same time!
% nodes_true =  data .* node_data < 0 | ( data==0 & node_data~=0) ;
nodes_true = sparse( dataCut .* node_data < 0 );  % In paper this is <= !! Here the =0 case is included in deltaReg!
h_x1Plus =  abs( node_data(nodes_true) ./ Dx1_fwd(nodes_true) ); 
h_x1Plus(h_x1Plus>=h_x1*facHmax) = 0; % if h_x1Plus is bigger than a cell, then do not consider it!
delta_x1_up(nodes_true) = abs( h_x1Plus .* n_x1(nodes_true) ) / h_x1^2;

% (3.1.2) delta x1_down (-x1)
node_data = data( stencil:end-stencil-1 , stencil+1:end-stencil );
nodes_true = sparse( dataCut .* node_data < 0 );
h_x1Minus =  abs( node_data(nodes_true) ./ Dx1_bwd(nodes_true) ); 
h_x1Minus(h_x1Minus>=h_x1*facHmax) = 0; % if h_x1Minus is bigger than a cell, then do not consider it!
delta_x1_down(nodes_true) = abs( h_x1Minus .* n_x1(nodes_true) ) / h_x1^2;


% (3.2) x2 direction
% (3.2.1) delta x2_up (+x2)
node_data = data( stencil+1:end-stencil , stencil+2:end-stencil+1 );
% don't allow data and node_data to be zero at the same time!
% nodes_true =  data .* node_data < 0 | ( data==0 & node_data~=0) ;
nodes_true = sparse( dataCut .* node_data < 0 );  % In paper this is <= !! Here the =0 case is included in deltaReg!
h_x2Plus =  abs( node_data(nodes_true) ./ Dx2_fwd(nodes_true) ); 
h_x2Plus(h_x2Plus>=h_x2*facHmax) = 0; % if h_x2Plus is bigger than a cell, then do not consider it!
delta_x2_up(nodes_true) = abs( h_x2Plus .* n_x2(nodes_true) )  / h_x2^2;


% (3.2.2) delta x2_down (-x2)
node_data = data( stencil+1:end-stencil , stencil:end-stencil-1 ); 
nodes_true = sparse( dataCut .* node_data < 0 );
h_x2Minus =  abs( node_data(nodes_true) ./ Dx2_bwd(nodes_true) ); 
h_x2Minus(h_x2Minus>=h_x2*facHmax) = 0; % if h_x1Minus is bigger than a cell, then do not consider it!
delta_x2_down(nodes_true) = abs( h_x2Minus .* n_x2(nodes_true) )  / h_x2^2;



%% (4) Summation of all terms
delta_sum = delta_x1_up + delta_x1_down + delta_x2_up + delta_x2_down + deltaReg;

% At boundaries the forward or backward difference can become 0 and thus
% the delta function might show some inf entries. Replace them by 0:
delta_sum( delta_sum == inf ) = 0;


% For debuggung
% e = full(delta_sum);
% sum1 = sum(e,1)*h_x1;
% sum2 = sum(e,2)*h_x2;

end

