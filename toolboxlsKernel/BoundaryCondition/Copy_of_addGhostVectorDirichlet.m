function dataOut = addGhostVectorDirichlet(dataIn, dim, width, ghostData)
% addGhostDirichletFlameHolder: add ghost cells with Dirichlet boundary conditions
%       and a fixed zero value at a desired position which holds the flame
%
%   dataOut = addGhostDirichletFlameHolder(dataIn, dim, width, ghostData)
%
% Creates ghost cells to manage the boundary conditions for the array dataIn.
%
% This m-file fills the ghost cells with constant data
%   (ie Dirichlet boundary conditions).
%
% At present, this code can only handle state (and time) independent
%   Dirichlet data, although the value can be different on the upper
%   and lower boundaries of the grid.
%
% Notice that the indexing is shifted by the ghost cell width in output array.
%   So in 2D with dim == 1, the first data in the original array will be at
%          dataOut(width+1,1) == dataIn(1,1)
%
% parameters:
%   dataIn	Input data array.
%   dim		Dimension in which to add ghost cells.
%   width	Number of ghost cells to add on each side (default = 1).
%   ghostData	A structure (see below).
%
%   dataOut	Output data array.
%
% ghostData is a structure containing data specific to this type of
%   ghost cell.  For this function it contains the field(s)
%
%   .lowerVector Vector which contains the values for the lower boundary. 
%     Row or column doesn't matter (no default).
%   .upperVector Vector which contains the values for the upper boundary. 
%     Row or column doesn't matter (no default).
%
% This boundary condition can only be applied to 2 dimensional problems
% Best practice to calculate the lower-/ upperVector is to derive it from the 
%   initial condition

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing
%   agreement contained in the file LICENSE in the top directory of
%   the distribution.
%
% Ian Mitchell, 1/13/04
%
% Edited by Tomas Steinbacher (tsteinbacher@td.mw.tum.de)
% 11/20/14

%% parse input
% Get dimension of the problem and check applicability of BC
dims = ndims(dataIn);
if dims ~= 2
  error('Boundary condition addGhostDirichletFlameHolder requires a 2 dimensional geometry!');
end

if(nargin < 3)
  width = 1;
end

if((width < 0) || (width > size(dataIn, dim)))
  error('Illegal width parameter');
end

if((nargin == 4) && isstruct(ghostData))
  % ghostdata structure exists: Read in provided data
  
  % lower flame position: Has to be provided!
  if(isfield(ghostData, 'lowerVector'))
    lowerVector = ghostData.lowerVector;
  else
    error('ghostData structure must contain field lowerFlamePosIndx');
  end
  
  % upper flame position: if not provided same as lower one
  if(isfield(ghostData, 'upperVector'))
    upperVector = ghostData.upperVector;
  else
    upperVector = lowerVector;
  end
  
else
  % ghostdata structure does not exist: Set default values
  error('No data provided for vector dirichlet boundary condition!')
  
end

% calculate dimension in which the vector is inserted to
if dim == 1
  % take next dimension
  insertDim = 2;
else
  % take preceding dimension
  insertDim = 1;
end


%% Create output and fill it with input data inside the domain
% create cell array with array size
sizeIn = size(dataIn);
indicesOut = cell(dims, 1);
for i = 1 : dims
  indicesOut{i} = 1:sizeIn(i);
end
% indicesIn = indicesOut;

% create appropriately sized output array
sizeOut = sizeIn;
sizeOut(dim) = sizeOut(dim) + 2 * width;
dataOut = zeros(sizeOut);

% fill output array with input data
indicesOut{dim} = width + 1 : sizeOut(dim) - width;
dataOut(indicesOut{:}) = dataIn;

%% fill in the ghost cells
% create vector for boundary condition: Check right dimensions!
[r_l, c_l] = size(lowerVector);
[r_u, c_u] = size(upperVector);

if insertDim == 1
  % A row vector is needed
  if r_l < c_l
    lowerVector = lowerVector';
  end
  if r_u < c_u
    upperVector = upperVector';
  end
  lowerBoundary_values = repmat( lowerVector , 1 , width );
  upperBoundary_values = repmat( upperVector , 1 , width );
  
else
  % A column vector is needed
  if r_l > c_l
    lowerVector = lowerVector';
  end
  if r_u > c_u
    upperVector = upperVector';
  end
  lowerBoundary_values = repmat( lowerVector , width , 1 );
  upperBoundary_values = repmat( upperVector , width , 1 );
  
end


% (a) lower
indicesOut{dim} = 1 : width;
% cells to add along boundary (neccessary when curvature is activated)
n_domain = length(indicesOut{insertDim});                 % beta!
delta = ( n_domain - length(lowerBoundary_values) ) / 2;  % beta!
indicesOut{insertDim} = 1+delta : n_domain-delta;         % beta!
dataOut(indicesOut{:}) = lowerBoundary_values;

% (b) upper
indicesOut{dim} = sizeOut(dim) - width + 1 : sizeOut(dim);
% cells to add along boundary (neccessary when curvature is activated)
n_domain = length(indicesOut{insertDim});                 % beta!
delta = ( n_domain - length(upperBoundary_values) ) / 2;  % beta!
indicesOut{insertDim} = 1+delta : n_domain-delta;         % beta!
dataOut(indicesOut{:}) = upperBoundary_values;
