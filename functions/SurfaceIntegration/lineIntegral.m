function [ integralLine , data ] = lineIntegral( grid, data , myFun , s  )
%LINEINTEGRAL calculates the line integral of function fun along the zero
%level iso line of data
%
% Only 2D!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Calculate delta function 
if s.orderDelta == 1
  %(first order approximation)
  [ delta_sum , stencil , data ] = deltaFun2D_O1( grid , data , s );
else
  %(second order approximation)
  [ delta_sum , stencil , data ] = deltaFun2D_O2( grid , data , s );
end

% Multiply delta function with desired function fun 
if length(myFun)>1
  % Fun is matrix
  integralLine = delta_sum .* myFun( stencil+1:end-stencil , stencil+1:end-stencil ) ;
else
  % Fun is scalar
  integralLine = delta_sum .* myFun;
end

% Sum up and convert to full (otherwise it would be sparse)
integralLine = full( sum( integralLine(:) ) * grid.dx(1)*grid.dx(2) );

% For debugging spikes:
% if integralLine >0.5
%   disp('Wait...')
% end

end



