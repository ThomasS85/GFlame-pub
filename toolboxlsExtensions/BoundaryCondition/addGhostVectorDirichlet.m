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
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


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
  
  if(ghostData.towardZero)
    slopeMultiplier = -1;
  else
    slopeMultiplier = +1;
  end
  
else
  % ghostdata structure does not exist: Set default values
  error('No data provided for vector dirichlet boundary condition!')
  
end


%% Create output and fill it with input data inside the domain
% create cell array with array size
sizeIn = size(dataIn);
indicesOut = cell(dims, 1);
for i = 1 : dims
  indicesOut{i} = 1:sizeIn(i);
end
indicesIn = indicesOut;

% create appropriately sized output array
sizeOut = sizeIn;
sizeOut(dim) = sizeOut(dim) + 2 * width;
dataOut = zeros(sizeOut);

% fill output array with input data
indicesOut{dim} = width + 1 : sizeOut(dim) - width;
dataOut(indicesOut{:}) = dataIn;


%% fill in the ghost cells
% (a) lower: Fixed value (Sirichlet) + Extrapolation
% compute slopes
indicesIn{dim} = 1;
slopeBot = lowerVector - dataIn(indicesIn{:});
% adjust slope sign to correspond with sign of data at array edge
% slopeBot = slopeMultiplier * abs(slopeBot) .* sign(dataIn(indicesIn{:}));
% Set first values to dirichlet vector
indicesOut{dim} = width;
dataOut(indicesOut{:}) = lowerVector;
% now extrapolate
for i = 1 : width-1
  indicesOut{dim} = i;
  dataOut(indicesOut{:}) = (lowerVector + (width - i + 1) * slopeBot);
end


% (b) upper: Fixed value (Sirichlet) + Extrapolation
% compute slopes
indicesOut{dim} = sizeIn(dim);
slopeTop = upperVector - dataIn(indicesOut{:});
% adjust slope sign to correspond with sign of data at array edge
% indicesIn{dim} = sizeIn(dim);
% slopeTop = slopeMultiplier * abs(slopeTop) .* sign(dataIn(indicesIn{:}));
% Set first values to dirichlet vector
indicesOut{dim} = sizeIn(dim) + width + 1;
dataOut(indicesOut{:}) = upperVector;
% now extrapolate
for i = 1 : width - 1
  indicesOut{dim} = sizeOut(dim) - i + 1;
  dataOut(indicesOut{:}) = (upperVector + (width - i + 1) * slopeTop);
end


end