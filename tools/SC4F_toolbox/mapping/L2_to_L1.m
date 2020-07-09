function [ x_L1 ] = L2_to_L1( x_L2 , p )
%L1_to_L2 function maps from L2 to L1 coordinate system
%
% Inputs:
%   - x_L2   :  Complex vector in L2 coordinate system
%   - p      :  Struct with flame settings as returned from setUpPredefinedFlame()
%
% Outputs:
%   - x_L1   :  Complex vector in L1 coordinate system
%
% Coordinate System L1:
%    ^
%    |  |
% x2 |-- o                     x_L(1)
%    |      o                    :
%    |         o                 :
%    |            o              :
%  0 |__ . __ . __ . __ . __   x_L(end)
%    |
%    ---------------------> 
%       0                x1
%
% % Coordinate System L2:
%    ^
% x2 |__ . __ . __ . __ . __  x_L(1)
%    |          o               :
%    |       o                  :
%    |--- o                     :
%  0 |   |_________________    x_L(end)
%    |
%    ---------------------> 
%        0               x1
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 09.11.2017                                //
% // Last modified: 09.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////

if strcmpi(p.CombType,'backwardFacingStep')
  x_L1 = real(x_L2) + 1i * ( p.R_a - imag(x_L2) );
else
  x_L1 = x_L2;
end

end

