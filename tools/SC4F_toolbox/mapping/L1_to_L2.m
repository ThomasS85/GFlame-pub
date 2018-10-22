function [ x_L2 ] = L1_to_L2( x_L1 , p )
%L1_to_L2 function maps from L1 to L2 coordinate system
%
% Inputs:
%   - x_L1   :  Complex vector in L1 coordinate system
%   - p      :  Struct with flame settings as returned from setUpPredefinedFlame()
%
% Outputs:
%   - x_L2   :  Complex vector in L2 coordinate system
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
  x_L2 = real(x_L1) + 1i * ( p.R_a - imag(x_L1) );
else
  x_L2 = x_L1;
end

end

