function [f_vec , t_vec , f_orig ] = generateSpecialInputSignalGFLAME( vel , p , varargin )
%GENERATESINESWEEPGFLAME Function generates a sine sweep or broad band signal depending on velocity being
% either transientDef=sweep or transientDef=broad. Signal creation is based on user settings and maximum time lag
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.de).        //
% // Created, 07.04.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Estimate some important features for Transfer Funtion
f_upperLimit = 1e4;
% Estimate Cut off frequency
if strcmp(p.flameType,'conical')
  % f_cutOff = 1 / convective time scale (restoration)
  f_cutOff = 1 / p.tau_c;
  
elseif strcmp(p.flameType,'slit') && strfind(vel.velModel,'convective')
  % Difference between time of restoration and distortion
  tau_d = p.H_flame / vel.K;
  f_cutOff = min( 1 / abs(tau_d - p.tau_c) ,f_upperLimit ) / 2;  % the /2 is empirical
  
else
  % Otherwise take conical value as default
  f_cutOff = 1 / p.tau_c;
  
end

% f_max is 3x higher than cut off
f_max = 3 * f_cutOff;

% f_min is at 2% of f_max
f_min = max( 0.02 * f_max , 5 );



%% Generate signal
if strcmpi(vel.transientDef,'sweep')
  % Sweep
  
  % if user specified min/ max frequency than take this value
  if vel.f_min > 0
    f_min = vel.f_min;
  end
  if vel.f_max > 0
    f_max = vel.f_max;
  end
  
  % Create frequencies
  f_orig = linspace( f_min , f_max , vel.n_sweep );
  % Define Ts so that for highest frequency we still have 15 points / period
  Ts = 1 / (15*f_max);
  % Estimated length of IR
  tTransient = p.tau_c * 1.2;
  % Use TFD tools function to generate sine sweep
  [data, tEnd] = siCreateSineSweep(f_orig,Ts,vel.nPeriods,tTransient,'uref');
  
  % Export signal to vectors
  t_vec = data.SamplingInstants';
  f_vec = data.u';
  
  % Save sweep data
  if ~(exist(p.run_output_folder,'dir')==7)
    mkdir(p.run_output_folder) 
  end
  save([p.run_output_folder,filesep,'sweep_data.mat'] , 'f_vec' , 't_vec' , 'tEnd','tTransient','f_orig')

  
elseif strcmpi(vel.transientDef,'broad')
  % Broad
  
  % Desired maximum frequency
  if vel.f_maxBB > 0
    f_max = vel.f_maxBB;
  end
  
  % Desired length of input signal (impulse length times vel.n_IR times security factor 1.2)
  t_max = p.tau_c * 1.2 * vel.n_IR;
  % Nyquist frequency for desired maximum frequency of interest
  f_NQ = 2 * f_max;
  % Therfor we need a sample time of Ts
  Ts = 1 / f_NQ;
  % Thus the length of the signal vector must be
  N_BB = ceil( t_max / Ts );
  % Create time vector
  t_vec = linspace( 0 , t_max , N_BB );
  % Create signal
  f_vec = idinput( N_BB , 'prbs' , [0 1] , [-1 1] )';
  
  % Return maximum frequency 
  f_orig = f_max;
  
  
end


%% Get right starting time and write signal to files
% shift time vetor so that it starts at t0
t_vec = ( t_vec - t_vec(1) ) + vel.t_transient0;

% Only write data when no additional input was made
if nargin == 2
  % write to uInlet-file
  f1 = fopen([p.run_output_folder,filesep,'uInlet'],'w');
  fprintf(f1,'xData , yData');
  fprintf(f1,'\n%.8f , %.8f',[t_vec;f_vec]);
  fclose(f1);

end

% debug plot
% f = figure;
% plot(t_vec, f_vec);
% xlabel('t [s]')
% title('Sweep Input Signal')
% close(f)


end

