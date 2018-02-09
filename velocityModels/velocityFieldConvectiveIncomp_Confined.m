function [ vel , schemeData , u_ref ] = velocityFieldConvectiveIncomp_Confined( t, data, schemeData )
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

p = schemeData.p;
d_x1 = schemeData.grid.dx(1);

if strcmpi(p.flameType,'conical')
  % CONICAL: Mean velocity field from Cuquel Diss p. 133 (conical)
  vel{1} = schemeData.meanFlowSpeed * ( 1 + p.gamma_s/p.H_flame * schemeData.grid.xs{1} );
  vel{2} = -schemeData.meanFlowSpeed * p.gamma_s / (2*p.H_flame) * schemeData.grid.xs{2};
  fg = 0.5;
else
  % SLIT: Mean velocity field from MUPAD
  vel{1} = schemeData.meanFlowSpeed * ( 1 + p.gamma_s/p.H_flame * schemeData.grid.xs{1} );
  vel{2} = -schemeData.meanFlowSpeed * p.gamma_s / (p.H_flame) * schemeData.grid.xs{2};
  fg = 1;
end


if strcmp(schemeData.type,'transient')
   % transient velocity
  
  % Calculate distributed time delay for time t (with respect to reference position)
  tau_vec = schemeData.grid.xs{1} / schemeData.K;
  
  % Now create a time vector where actual time t is at reference position and corresponding to the computet
  % time delays past times for the other x-coordinates
  tau_vec = t - tau_vec;
  
  % Evaluate velocity in x1 direction with desired Kernel
  [ u_1s , u1s_ref ] = applyKernel_convectiveMods( data , schemeData , tau_vec );
  
  % Compute spatial derivative of u_1s using Savitzky Golay Filter
  degreePoly = 3; frameSize = 2*degreePoly+1; % Settings see savitzkyGolayFilt
  du_1s_dx1 = -savitzkyGolayFilt( u_1s(:,1) , degreePoly , 1 , frameSize , [] , 1 ) / d_x1^1;
  
  % Debug plot:
%   f_tmp=figure; yyaxis left; plot(schemeData.grid.vs{1},u_1s(:,1),'b','DisplayName','u');
%   hold on; yyaxis right;  plot(schemeData.grid.vs{1},du_1s_dx1,'r:','DisplayName','du/dx');
%   xlabel('x [m]');legend('show')
%   close(f_tmp);
  
  du_1s_dx1 = repmat( du_1s_dx1 , 1 , schemeData.grid.N(2) );
  
  % Evaluate radial velocity field: u_2
  
  u_2s = -fg * schemeData.grid.xs{2} .* du_1s_dx1;
  
  % Writes reference velocity 
  u_ref = u1s_ref + vel{1}(1,1);
  
  % Copy velocity field to output
  vel{1} = u_1s + vel{1};
  vel{2} = u_2s + vel{2};
  
else
  % Writes reference velocity 
  u_ref = vel{1}(1,1);
  
end


end

