function [ v ] = settings_VelocityModel_Example( p )
%settings_VelocityModell Settings for the desired velocity model
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 10.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Basic settings
% Type?  'convectiveConfined' or 'convective' or 'convectiveIncomp' or 'uniform' or 'uniform_transversal'
%         or 'convectiveIncompConfined'
v.velModel = 'convective';
% 'transient' or 'constant'?
v.type = 'transient';
% if transient: temporal evolution defined by 
%   'file'
%   'harmonic'
%   'sweep'
%   'broad'
%   'step'
%   'TaX'
% You should not change transientDef when resuming a case (past excitation data might be wrong)!
v.transientDef = 'harmonic';


%% General settings
% Speed of mean flow (see settings_Example for definition of u_1_L)
v.meanFlowSpeed = p.u_1_L;        % [m/s] 
% Convective velocity of the distortion (only for convective velocity model)
v.K = v.meanFlowSpeed;            % [m/s] 
% Amplitude of transient part of the signal
v.vAmp = 0.1;
% When should transient part of the signal start? Only for step, harmonic, sweep and broad
v.t_transient0 = 0;
% Definition of Kernel for convective velocity models: Dirac or Diffusion
v.kernel = 'Dirac';


%% Kernel Settings
% Set sigma of kernel (Diffusion Kernel only)
v.kernelWidth = 0.002;
% Set kinematic viscosity (Diffusion Kernel only; for v.mu=0 the kernelWidth inst constant)
v.mu = 1.54e-5; % [m^2/s] Air-methane mixture at 293K (cantera, Gri-30)
v.mu = 3.93e-4; % [m^2/s] Air-methane mixture at 2008K (cantera, Gri-30)


%% Settings for velocity from file
% Name of file which contains a list of transient behaviour
v.path2File = [p.run_output_folder,'/uInlet'];
% Should an openFOAM list be used as input? (1 or 0)
v.useFOAM_list = 1;
% file name of openFoam list
v.fileFOAM = 'uInlet_FOAM';
% resample this list? dx=0 -> No   ; dx>0 : specifies resample spacing
v.dx = 1e-4;


%% Settings for velocity harmonic
% Frequency of the signal
v.frequ = 80;       % [Hz]


%% Settings for velocity sweep
% Creates a sine sweep with n_sweep frequencies
v.n_sweep = 60;
% number of periods per frequency
v.nPeriods = 4;
% Specify maximum/ minimum frequency [Hz]. If negative value is set the value is estimated
v.f_min = -10;
v.f_max = -1500;


%% Settings for velocity broad
% Length of the broad band signal is defined in n_IR * IR length
v.n_IR = 15;
% Maximum frequency [Hz] (negative value if this should be estimated)
v.f_maxBB = -1;
% Identification: Number of FIR coefficients? (vector if several should be tried)
v.FIR_ord = [20 15];
% Define time delay [s] (if negative this is estimated)
v.tauBB = -1;


%% Settings for velocity step
v.sharpness = 1e5;


%% Settings for velocity TaX
% mean flame surface area (take steady state mean flame area)
v.A_mean = 0.001; % [m^2]
% Path to *.mat file which contains the acoustic SISO system
v.path2sys = '/home/steinbacher/dataLocal/data/MATLAB_Projects/TaX_projects/Kornilov_setup/kornilovSys.mat';
% CFL number for acoustics (TaX model)
v.CFL_ac = 0.5;

% Amplitude v.vAmp is automaticaly set to 1 (no scaling of transient input data from TaX)



%% Settings for velocity from CFD data
% Path to desired CFD case
v.path2Case = '/home/steinbacher/dataLocal/CFD/Duchaine/rhoReactingFoam/MonoFrequent/000DrRF00MF107';
% define which fields to load
v.fields2Load = {'U','T'};
% Define what time step to load (in case v.type is set to constant). 
%   Set to a negative value to load the first available time step
v.loadTime = 0.142461;
% v.loadTime = -1;
% Origin of CFD data coordinates in the coordinate system used here (use
% [0,0] if identical!)
v.CFDorigin = [0 0];
% If CFD data set is huge only save every nth value 
v.dSave = 1 ;
% Reload CFD data if case folder is already existent in output folder?
v.reloadData = 1;


end

