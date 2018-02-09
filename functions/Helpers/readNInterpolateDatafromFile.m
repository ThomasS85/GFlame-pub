function [ y_interp ] = readNInterpolateDatafromFile( x , myFile , varargin )
%READNINTERPOLATEDATAFROMFILE Reads and interpolates data from a file
%
% Reads in data from a 2 column file (first line is skipped). 1st column
% are x-values 8in rising order), 2nd are y-values:
%
%     tilte x   ,  title y
%     0.1       ,  0.354
%     0.2       ,  0.456
%     0.3       ,  0.25
%     0.4       ,  3.354
%     ...       ,  ...
%
% Default deleminter is ','. if you want a different one specify
%     ('delim','.')  (here '.' is specified)
%
% Based on this it interpolates the x value from the input linearly
% No extrapolation is performed unless the option 'extrap' is specified!
%
% Function can handle vector or matrix input for x! If x is a matrix it is assumed that all rows are the same
% and thus only the first column is taken into account. However, the output is again a matrix (using repmat)
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Set delimiter
ind = find(strcmp(varargin,'delim'),1);
if ~isempty(ind)
  delim = varargin{ind+1};
else
  delim = ',';
end
% Perform extrapolation: All data outside domain is set to a specific value
ind = find(strcmp(varargin,'extrap'),1);
if ~isempty(ind)
  extrap = 1;
  extrapVal = varargin{ind+1};
else
  extrap = 0;
  extrapVal = 0;
end


%% Check if input data x is matrix or vector. If matrix only take the first row since this is relevant for time delay
[~,s] = size(x);
if s>1
  % Only take first column
  x = x(:,1);
end

%% Read and Interpolate data
% Read data to M
try
  M = dlmread(myFile, delim , 1, 0);
catch
  M = [0 , extrapVal ; 1 , extrapVal ];
end

if extrap
  % Check length of M. If it contains only one line add a second in order
  % for interpolation to work
  if size(M,1)<=1
    M = [-1 , extrapVal ; M ];
  end
  % interpolate (and extrapolate to 0)
  y_interp = interp1( M(:,1) , M(:,2), x ,'linear',extrapVal);
  
else
  % check if x is inside the specidied region (no extrapoltion!)
  if min(x) < M(1,1) || max(x)>M(end,1)
    error('No extrapolation supported!')
  end
  
  % interpolate
  y_interp = interp1( M(:,1) , M(:,2), x );
end

%% Output data correctely
if s>1
  % If input was matrix also return matrix of same size
  y_interp = repmat(y_interp,1,s);
end

end

