function [ y ] = myHeaviside( x )
%MYHEAVISIDE Modified Heaviside function wich is exactely the heaviside
%function except at x=0. Here, a y of 0 is returned:
% 
%   Inputs:
%      - x  : Scalar, vector or matrix
%
%   Outputs:
%      - y  : Scalar, vector or matrix (according to input)
%
%        |- x < 0   : y = 0
%   y = -|- x = 0   : y = 0
%        |- x > 0   : y = 1
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////

y = zeros(size(x));
y(x>0) = 1;

end

