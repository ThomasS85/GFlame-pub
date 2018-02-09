function [ vel , schemeData , u_ref ] = velocityFieldExternal( t, ~, schemeData,varargin  )
%VELOCITYFIELDEXTERNAL Uses an external constant velocity field as input
%   
% Use the provided tools to create a conformal velociy field which can be
% used by this function!
% 
% function [ vel ] = velocityField( t, data, schemeData )
%
% schemeData:  
%       .dim              : Dimension
%       .externalVelData  : External velocity data as matrix: 
%                             (x1,x2,velocity component) 
%       .fileV          : Handle to output file
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (tsteinbacher@td.mw.de).        //
% // Created, 15.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////



%% Define complete velocity field
% Constant velocity field
vel{1} = schemeData.externalVelData(:,:,1);
vel{2} = schemeData.externalVelData(:,:,2);

% Writes reference velocity 
u_ref = vel{1}(1,1);

end

