function [ SR ] = IR2SR( IR , t_vec , varargin )
%IR2STEP Calculates the step response SR from a give impulse response IR
%
%
% - Function interpretes rows as time series
% - Matrix input:
%            x
%      ------>
%     |
%     |
%   t V
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2015 as part of GFLAME 0.1          //
% // Last modified: 08.05.2015 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Set filter parameter
ind = find(strcmpi(varargin,'filter'),1);
if ~isempty(ind)
  % Take filter parameter from varargin
  applyFilter = true;
  % Filter Parameter
  filterParam = varargin{ind+1};
  % Polynomial degree
  degreePoly = varargin{ind+2};
  % Frame size
  frameSize = 2*filterParam + 1; % Must be odd!
else
  % Default filter parameter
  applyFilter = false;
end



%% Check Input (column vector needed)
[l,c] = size(IR);
if c>l && l==1
  % If input is a vector also interprete line vector as time series
  IR = IR.';
  inputTransposed = true;
  l = c;
  c = 1;
else
  inputTransposed = false;
end


%% Iterate over columns and smooth, fit spline and derive
SR = zeros(l,c);
for ii = 1:c
  % Apply filter for smoothing if desired
  if applyFilter
    IR(:,ii) = savitzkyGolayFilt( IR(:,ii) , degreePoly , 0 , frameSize , [] , 1 ) ;
  end
  
  % Fit spline and integrate it
  IR_spline = spline( t_vec , IR(:,ii) );
  SR_spline = fnint(IR_spline);
  
  % Evaluate Integral
  SR(:,ii) = fnval(SR_spline,t_vec);
  
end

%% Output should be same format as input
if inputTransposed
  SR = SR.';
end

end

