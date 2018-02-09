function [ vel , schemeData , u_ref ] = velocityFieldUniformTransversal( t, data, schemeData )
%VELOCITYFIELD Definition of a velocity field; settings in separate file!
%
% function [ vel ] = velocityField( t, data, schemeData )
%
% Inputs:
%   - t             : Actual time
%   - data          : G-field (matrix)
%   - schemeData    : convectiveData struct from termConvection
%       .dim            : Dimension
%       .type           : Type (transient or constant)
%       .grid           : Grid (tollboxls)
%       .K              : Speed of convected distortions
%       .transientDef   : Definition for transient data (file, harmonic or
%                         step)
%       .vAmp           : Amplitude of distortion
%       .meanFlowSpeed  : Mean flow spead
%
%     Only for specific settings:
%       .t_transient0   : Time when transient part starts
%       .path2File      : Path to input file with velocity data
%       .omega          : angular frequency of harmonic distortion
%       .sharpness      : Sharpness of the step function
%
% Outputs:
%   - vel           : velocity field
%   - schemeData    : Modified version of input
%   - u_ref         : Reference Velocity: Scalar
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.de).        //
% // Created, 09.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////
  

%% Define complete velocity field
% initialise 2D cell vector
vel = zeros(schemeData.dim,1);
vel = num2cell(vel);

% Set constant velocity in x1 direction
vel{1} =  schemeData.meanFlowSpeed;

% Calculate velocity in x2 direction
[ u_s_rel,schemeData] = getTransientInput( schemeData , t , data );
vel{2} = schemeData.vAmp * u_s_rel;

% Writes reference velocity 
u_ref = vel{2}(1,1);

end

