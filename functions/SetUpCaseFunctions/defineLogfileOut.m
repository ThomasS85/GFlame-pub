function [ p,solver,vel ] = defineLogfileOut( p,solver,vel )
%DEFINELOGFILEOUT Defines information which is written 2 logfile
%

% write fields to log file (only numeric or string!)

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% General settings
p.logfile = { 'flameType' , 'domainType' , 's_l_u' , 'marksteinLength', 'nu' , ...
  'u_1_L','E','gamma','R_i','R_a','Cr','R_flame','H_flame','L_flame','alpha','alphaDegree',...
  'liftOff','lateralOffset','logFileName'};

%% Solver settings
solver.logfile = { 'dim' , 'curvature' , 'initial' ,'purgeWrite', 't0', 'tMax' , ...
  'accuracy','dxi','x1Lim','factorCFL','sigmaG','SurfInt_method','stats','outSteps','singleStep','doPlot',...
  'intervalWriteOI','sigmaG','doWriteData'};

%% Reinitialization settings (solver)
solver.reinit.logfile = {'do','intervalStep','accuracy','tMax','errorMax'};

%% Velocity settings
vel.logfile = { 'velModel' , 'type' , 'transientDef' , 'meanFlowSpeed' , 'kernel' , 'kernelWidth' , 'mu' };

% dependent on chosen type:
if strfind(vel.velModel,'convective') 
  vel.logfile = {vel.logfile{:} , 'K'};
  
elseif strcmp(vel.velModel,'CFD_FOAM')
  vel.logfile = {vel.logfile{:} , 'path2Case','CFDorigin','dSave'};
  
end

% dependent on chosen model:
if strcmp(vel.transientDef , 'file' )
  vel.logfile = { vel.logfile{:} , 'vAmp','path2File','useFOAM_list','fileFOAM'};
  
elseif strcmp(vel.transientDef , 'harmonic' )
  vel.logfile = {vel.logfile{:} , 'vAmp','frequ','omega','t_transient0'};
  
elseif strcmp(vel.transientDef , 'sweep' )
  vel.logfile = {vel.logfile{:} , 'vAmp','t_transient0','n_sweep','nPeriods','f_min','f_max'};
  
elseif strcmp(vel.transientDef , 'broad' )
  vel.logfile = {vel.logfile{:} , 'vAmp','t_transient0','n_IR','f_maxBB'};
  
elseif strcmp(vel.transientDef , 'step' )
  vel.logfile = {vel.logfile{:} , 'vAmp','t_transient0','sharpness'};
  
elseif strcmp(vel.transientDef , 'TaX' )
  vel.logfile = {vel.logfile{:} ,'A_mean','path2sys','CFL_ac'};
  
end

end

