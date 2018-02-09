function [ vel , schemeData , u_ref ] = velocityFieldConvective( t, data, schemeData )
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

if strcmp(schemeData.type,'transient')
  % transient velocity
  
  % Calculate distributed time delay for time t (with respect to area jump)
  tau_vec = schemeData.grid.xs{1} / schemeData.K;
  
  % Now create a time vector where actual time t is at reference position and corresponding to the computed
  % time delays past times for the other x-coordinates
  tau_vec = t - tau_vec;
  
  % Evaluate velocity in x1 direction with desired Kernel
  [ u_1s , u1s_ref ] = applyKernel_convectiveMods( data , schemeData , tau_vec );
  
  % Evaluate velocity field
  vel_field = schemeData.meanFlowSpeed + u_1s;

  % Copy velocity field to output
  vel{1} = vel_field;
  
  % Writes reference velocity 
  u_ref = u1s_ref + schemeData.meanFlowSpeed;
  
else
  % constant velocity
  vel{1} = schemeData.meanFlowSpeed;
  
  % Writes reference velocity 
  u_ref = schemeData.meanFlowSpeed;
  
end


end

