function dataOut = addGhostVectorDirichletSym(dataIn, dim, width, ghostData)
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
%   .symmetryPlane String which defines at what boundary symmetry can be 
%     assumde ('lower or 'upper'). The opposite boundary is a Dirichlet 
%     one (default = 'upper').
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
%
% Jul 2018: Instead of extrapolating for ghost cells at Dirichlet side, just copy the boundary vector
%           -> This would correspond to a fixed Neumann BC
%


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
  
  % Position of symmetry plane: Has to be provided!
  if(isfield(ghostData, 'symmetryPlane'))
    symmetryPlane = ghostData.symmetryPlane;
  else
    error('ghostData structure must contain field symmetryPlane');
  end
  
  % lower flame position: Has to be provided if symmetry plane is upper
  if(isfield(ghostData, 'lowerVector'))
    lowerVector = ghostData.lowerVector;
  elseif strcmp(symmetryPlane,'upper')
    error('ghostData structure must contain field lowerFlamePosIndx (or no valid symmetry plane chosen)!');
  end
  
  % upper flame position: Has to be provided if symmetry plane is lower
  if(isfield(ghostData, 'upperVector'))
    upperVector = ghostData.upperVector;
  elseif strcmp(symmetryPlane,'lower')
    error('ghostData structure must contain field lowerFlamePosIndx (or no valid symmetry plane chosen)!');
  end
  
  if(ghostData.towardZero)
    slopeMultiplier = -1;
  else
    slopeMultiplier = +1;
  end
  
  
else
  % ghostdata structure does not exist: Set default values
  error('No data provided for vector dirichlet boundary condition (defaults not available)!')
  
end


% Function only works for x2 direction!
if dim~=2
  error('Boundary Condition only implemented for x2 direction!')
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
% create vector for boundary condition: Check right dimensions!

if strcmp(symmetryPlane,'upper')
  
  % In case both boundary conditions are called one after another, cells have also been added in x1 direction.
  % This means the dirichlet vector would have to be extended by those number of cells. (e.g. by hessianSecond())
  n_lowerV = length(lowerVector);
  n_dataIn = size(dataIn,1);
  if n_dataIn > n_lowerV
    % add additional cells by linear extrapolation
    addCellsN = abs(n_lowerV - n_dataIn) / 2;
    lowerVector = [ interp1(lowerVector,0:-1:-addCellsN+1,'linear','extrap') ; lowerVector ; ...
      interp1(lowerVector,n_lowerV+1:n_lowerV+addCellsN,'linear','extrap') ];
  end
  
  
  % (a) lower: Fixed value (Dirichlet) + Extrapolation/ NOW NEUMANN!
  % compute slopes
%   indicesIn{dim} = 1;
%   slopeBot = lowerVector - dataIn(indicesIn{:});
  % adjust slope sign to correspond with sign of data at array edge
  % slopeBot = slopeMultiplier * abs(slopeBot) .* sign(dataIn(indicesIn{:}));
  % Set first values to dirichlet vector
  indicesOut{dim} = width;
  dataOut(indicesOut{:}) = lowerVector;
  % now extrapolate
  for i = 1 : width-1
    indicesOut{dim} = i;
%     dataOut(indicesOut{:}) = (lowerVector + ( width - i ) * slopeBot);
    dataOut(indicesOut{:}) = lowerVector; % no extrapolation, but just copy Dirichlet BC (=fixed Neumann!)
  end
  
  % (b) upper: Neumann (extrapolate)
  for i = 1 : width
    indicesOut{dim} = sizeOut(dim) - i + 1;
    indicesIn{dim} = sizeIn(dim);
    dataOut(indicesOut{:}) = (dataIn(indicesIn{:}));
  end
  
  
else
  
  % In case both boundary conditions are called one after another, cells have also been added in x1 direction.
  % This means the dirichlet vector would have to be extended by those number of cells. (e.g. by hessianSecond())
  n_upperV = length(upperVector);
  n_dataIn = size(dataIn,1);
  if n_dataIn > n_upperV
    % add additional cells by linear extrapolation
    addCellsN = abs(n_upperV - n_dataIn) / 2;
    upperVector = [ interp1(upperVector,0:-1:-addCellsN+1,'linear','extrap') ; upperVector ; ...
      interp1(upperVector,n_upperV+1:n_upperV+addCellsN,'linear','extrap') ];
  end
  
  % (a) lower: : Neumann (extrapolate)
  for i = 1 : width
    indicesOut{dim} = i;
    indicesIn{dim} = 1;
    dataOut(indicesOut{:}) = (dataIn(indicesIn{:}));
  end
  
  % (b) upper: Fixed value (Dirichlet) + Extrapolation/ NOW NEUMANN!
  % compute slopes
%   indicesIn{dim} = sizeIn(dim);
%   slopeTop = upperVector - dataIn(indicesIn{:});
  % adjust slope sign to correspond with sign of data at array edge
  % indicesIn{dim} = sizeIn(dim);
  % slopeTop = slopeMultiplier * abs(slopeTop) .* sign(dataIn(indicesIn{:}));
  % Set first values to dirichlet vector
  indicesOut{dim} = sizeIn(dim) + width + 1;
  dataOut(indicesOut{:}) = upperVector;
  % now extrapolate
  for i = 2 : width
    indicesOut{dim} = sizeIn(dim) + width + i;
%     dataOut(indicesOut{:}) = (upperVector + ( i - 1 ) * slopeTop);
    dataOut(indicesOut{:}) = upperVector; % no extrapolation, but just copy Dirichlet BC (=fixed Neumann!)
  end
  
end
