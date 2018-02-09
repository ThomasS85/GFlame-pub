function [ t_vec , uref_vec , A_surf_vec ] = getLastNTimeStepsFromFile( fileName_A , fileName_u , t0 , t1 , N )
%GETLASTNTIMESTEPSFROMFILE Function opens and reads from output file the flame surface area and reference
% velocity and outputs time series of both for t0 to t1 with N elements (linearily interpolated)
%
% Inputs:
%   fileName_A  : Filename with flame surface data A_surf
%   fileName_u  : Filename with reference velocity data uref
%   t0          : Starting time of output time series
%   t1          : Ending time of output time series
%   N           : Number of elements of output time series
%
% Outputs:
%   t_vec       : Output time series, starting with the latest time 
%   uref_vec    : Output reference velocity (interpolated to t_vec)
%   uref_vec    : Output flame surface area (interpolated to t_vec)
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 30.03.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% (1) Flame surface
% Open file
f_A_surf = fopen(fileName_A,'r');
% read header
fgetl(f_A_surf);
% Read data from file
data = cell2mat( textscan(f_A_surf,'%f , %f','Headerlines',0) );
t_vec_A = data(:,1);
A_surf_vec = data(:,2);
% Close file
fclose(f_A_surf);


%% (2) reference velocity
% Open file
f_uref = fopen(fileName_u,'r');
% read header
fgetl(f_uref);
% Read data from file
data = cell2mat( textscan(f_uref,'%f , %f','Headerlines',0) );
t_vec_u = data(:,1);
uref_vec = data(:,2);
% Close file
fclose(f_uref);


%% Interpolate to same time step
% Construct time vector
t_vec = linspace( t1 , t0 , N )';
% Interpolate data to time vector (extrapolation to zero )
A_surf_vec = interp1( t_vec_A , A_surf_vec , t_vec , 'linear' , 0 );
uref_vec = interp1( t_vec_u , uref_vec , t_vec , 'linear' , 0 );

end

