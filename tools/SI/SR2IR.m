function [ IR ] = SR2IR( SR , dt , varargin )
%STEP2IR Calculates the impulse response IR from a give step response SR
%
% Function needs savitzkyGolayFilt()
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
  % Filter Parameter
  filterParam = varargin{ind+1};
  % Polynomial degree
  degreePoly = varargin{ind+2};
else
  % Default filter parameter
  % Filter Parameter
  filterParam = 10;
  % Polynomial degree
  degreePoly = 3;
end
% Frame size
frameSize = 2*filterParam + 1; % Must be odd!


%% Check Input (column vector needed)
[l,c] = size(SR);
if c>l && l==1
  % If input is a vector also interprete line vector as time series
  SR = SR.';
  inputTransposed = true;
  l = c;
  c = 1;
else
  inputTransposed = false;
end


%% Derive SR using Savitzky Golay Filter
IR = zeros(l,c); 
for ii = 1:c
  IR(:,ii) = -savitzkyGolayFilt( SR(:,ii) , degreePoly , 1 , frameSize , [] , 1 ) / dt^1;
end

%% Output should be same format as input
if inputTransposed
  IR = IR.';
end

end

