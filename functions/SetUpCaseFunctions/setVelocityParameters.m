function [ vel ] = setVelocityParameters( vel , solver )
%SETVELOCITYPARAMETERS Calculates all parameters from user settings
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 10.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Store user settings in vel_tmp
vel_tmp = vel;
% Empty vel
vel = [];

%% Fill in general and basic settings which are used in each model
vel.velModel = vel_tmp.velModel;
vel.type = vel_tmp.type;
vel.transientDef = vel_tmp.transientDef;
vel.meanFlowSpeed = vel_tmp.meanFlowSpeed;
vel.kernel = vel_tmp.kernel;

% in case of constant velocity set transientDef to 'none'
if strcmp(vel.type,'constant')
  vel.transientDef = 'none';
end

% Check if Kernel is defined correctly
if ~ ( strcmpi(vel.kernel,'Dirac') ||  strcmpi(vel.kernel,'Diffusion') )
  error('The chosen Kernel is not defined! Please choose <Diffusion> or <Dirac> in the velocity settings file. Set v.mu to zero if you want a constant kernel width (no diffusion)')
end


if strcmpi(vel.kernel,'Diffusion')
  
  vel.mu = vel_tmp.mu;
  
  % Compute x1 length of domain
  L_x1 = ( solver.x1Lim(2) - solver.x1Lim(1) );
  if vel_tmp.kernelWidth >= L_x1
    error('Kernel Width is bigger than length of domain! Invalid settings.')
  elseif ~(vel_tmp.kernelWidth+eps>=solver.dxi*3)
    error(['Kernel width is smaller than twice the x-spacing dx1! Choose at least ',num2str(solver.dxi*3)])
  else
    vel.kernelWidth = vel_tmp.kernelWidth;
  end
end


%% Fill in specific settings for the different velocity models
% In case the confined velocity model is chosen check valididy
if strcmp(vel.velModel ,'convective') 
  % Convective velocity model
  vel.velocity = @velocityFieldConvective;
  vel.K = vel_tmp.K;
  % Postpone start of transient part by delay over the whole domain (x1)
%   vel_tmp.t_transient0 = vel_tmp.t_transient0 + solver.x1Lim(2) / vel.K;
  % Do not postpone in order to set transient starting point precisely
%   vel_tmp.t_transient0 = vel_tmp.t_transient0 ;

elseif strcmp(vel.velModel ,'convectiveIncomp')
  % Convective incompressible velocity model
  vel.velocity = @velocityFieldConvectiveIncomp;
  vel.K = vel_tmp.K;
  
elseif strcmp(vel.velModel ,'convectiveConfined')
  % Convective velocity model
  vel.velocity = @velocityFieldConvective_Confined;
  vel.K = vel_tmp.K;
  
elseif strcmp(vel.velModel ,'convectiveIncompConfined')
  % Convective velocity model
  vel.velocity = @velocityFieldConvectiveIncomp_Confined;
  vel.K = vel_tmp.K;
  
elseif strcmp(vel.velModel ,'uniform')
  % Uniform velocity model
  vel.velocity = @velocityFieldUniform;
  
elseif strcmp(vel.velModel ,'uniform_transversal')
  % Uniform velocity model
  vel.velocity = @velocityFieldUniformTransversal;
  
elseif strcmp(vel.velModel ,'SCcomputed')
  % Convective velocity model
  vel.velocity = @velocityFieldSCcomputed;
  vel.K = vel_tmp.K;
  
% elseif strcmp(vel.velModel ,'CFD_FOAM')
%   % Velocity from CFD: FOAM
%    vel.velocity = @velocityFieldExternal;
%    vel.transientDef = 'none';
%    % Pass settings
%    vel.path2Case = fileparts(vel_tmp.path2Case); % Extract only path
%    vel.fields2Load = vel_tmp.fields2Load;
%    vel.CFDorigin = vel_tmp.CFDorigin;
%    vel.dSave = vel_tmp.dSave;
%    vel.reloadData = vel_tmp.reloadData;
%    vel.loadTime = vel_tmp.loadTime;
%    
% elseif strcmp(vel.velModel ,'CFD_VTK')
%   % Velocity from CFD: VTK
%   vel.velocity = @velocityFieldExternal;
%   vel.transientDef = 'none';
%   % Pass settings
%   vel.path2Case = fileparts(vel_tmp.path2Case); % Extract only path
%   vel.fields2Load = vel_tmp.fields2Load;
%   vel.CFDorigin = vel_tmp.CFDorigin;
%   vel.dSave = vel_tmp.dSave;
%   vel.reloadData = vel_tmp.reloadData;
%   vel.loadTime = vel_tmp.loadTime;
  
% elseif strcmp(vel.velModel ,'VORTEX')
%   % Velocity from VORTEX code (vortex tracking)
%   vel.velocity = @velocityFieldVORTEX;
  
else
  error('Unknown velocity model chosen!')
  
end


%% Fill in specific parameters for the transient definitions
if strcmp(vel.transientDef , 'file' )
  % File as input
  vel.vAmp = vel_tmp.vAmp;
  [pathstr,name,ext] = fileparts(vel_tmp.path2File);
  vel.path2File = [pathstr,filesep,name,ext];
  vel.fileFOAM = vel_tmp.fileFOAM;
  vel.useFOAM_list = vel_tmp.useFOAM_list;
  vel.dx = vel_tmp.dx;
  vel.t_transient0 = vel_tmp.t_transient0;
  
elseif strcmp(vel.transientDef , 'harmonic' )
  % Transient behaviour: Harmonic function
  vel.vAmp = vel_tmp.vAmp;
  vel.frequ = vel_tmp.frequ;
  vel.omega = 2 * pi * vel.frequ;
  vel.t_transient0 = vel_tmp.t_transient0;
  vel.useFOAM_list = 0;

elseif strcmp(vel.transientDef , 'sweep' )
  % Transient behaviour: Sweep input
  vel.t_transient0 = vel_tmp.t_transient0;
  vel.vAmp = vel_tmp.vAmp;
  vel.n_sweep = vel_tmp.n_sweep;
  vel.nPeriods = vel_tmp.nPeriods;
  vel.f_min = vel_tmp.f_min;
  vel.f_max = vel_tmp.f_max;
  [pathstr,name,ext] = fileparts(vel_tmp.path2File);
  vel.path2File = [pathstr,filesep,name,ext];
  vel.useFOAM_list = 0;

elseif strcmp(vel.transientDef , 'broad' )
  % Transient behaviour: Sweep input
  vel.t_transient0 = vel_tmp.t_transient0;
  vel.vAmp = vel_tmp.vAmp;
  vel.n_IR = vel_tmp.n_IR;
  vel.FIR_ord = vel_tmp.FIR_ord;
  vel.f_maxBB = vel_tmp.f_maxBB;
  vel.tauBB = vel_tmp.tauBB;
  [pathstr,name,ext] = fileparts(vel_tmp.path2File);
  vel.path2File = [pathstr,filesep,name,ext];
  vel.useFOAM_list = 0;
  
elseif strcmp(vel.transientDef , 'step' )
  % Transient behaviour: Step function
  vel.vAmp = vel_tmp.vAmp;
  vel.t_transient0 = vel_tmp.t_transient0;
  vel.sharpness = vel_tmp.sharpness;
  vel.useFOAM_list = 0;

elseif strcmpi(vel.transientDef , 'Tax' )
  % Tax
  vel.useFOAM_list = 0;
  vel.A_mean = vel_tmp.A_mean;
  vel.path2sys = vel_tmp.path2sys;
  vel.CFL_ac = vel_tmp.CFL_ac;
  vel.vAmp = 1;               % Amplitude is 1 (no scaling of transient input data from TaX)  
  
elseif strcmp(vel.transientDef , 'none' )
  % constant
  vel.useFOAM_list = 0;
  
else
  error('Unknown definition of transient velocity field!')
end



end

