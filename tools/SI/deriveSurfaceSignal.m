function [ dAdt , t_new , A_new ] = deriveSurfaceSignal( GFcase , varargin)
%DERIVESURFACESIGNAL Derives signal of the velocity fluctuation applying a
% continuouse wavelet transform
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Should be plotted?
ind = find(strcmp(varargin,'plot'),1);
if ~isempty(ind)
  doPlot = 1;
else
  % Don't plot
  doPlot = 0;
end

% Smoothing parameter
ind = find(strcmp(varargin,'smoothing'),1);
if ~isempty(ind)
  scaleSmoothing = varargin{ind+1};
else
  % default value
  scaleSmoothing = 40;
end

% Number of samples for resampling
ind = find(strcmp(varargin,'resample'),1);
if ~isempty(ind)
  samples = varargin{ind+1};
else
  % Don't resample
  samples = 2000;
end

% Should reference velocity be plotted?
ind = find(strcmp(varargin,'plotRef'),1);
if ~isempty(ind)
  plotRef = 1;
  numSubplots = 3;
else
  % Don't plot
  plotRef = 0;
  numSubplots = 2;
end

% Should case be change?
ind = find(strcmp(varargin,'case'),1);
if ~isempty(ind)
  % Change run output folder
  GFcase.solver.resume2case = varargin{ind+1};
else
  % Load case specified by case number
  GFcase.solver.resume2case = GFcase.p.caseNumber;
end

% Should root folder be changed? 
ind = find(strcmp(varargin,'root'),1);
if ~isempty(ind)
  % Change root folder
  GFcase.p.caseRootDir = varargin{ind+1};
else
  % Leave everything as it is
end

% Limits
ind = find(strcmp(varargin,'lim'),1);
if ~isempty(ind)
  myLims = varargin{ind+1};
else
  % default value
  myLims = NaN;
end

% Norm derived time series to 1 so that they can be compared to IR?
ind = find(strcmp(varargin,'norm'),1);
if ~isempty(ind)
  doNorm = 1;
else
  % Don't plot
  doNorm = 0;
end

% Change dir to case root
change2caseRootDir( GFcase )


%% Read data
% define where to load data from
GFcase.p.run_output_folder = [GFcase.p.caseFolder,'/',num2str( GFcase.solver.resume2case )];

% Read surface data from file
f1 = fopen([GFcase.p.run_output_folder,'/Surface.out'],'r');
data_A = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
fclose(f1);

if plotRef
  % Read reference velocity data from file
  f1 = fopen([GFcase.p.run_output_folder,'/V_ref.out'],'r');
  data_V = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
  fclose(f1);
end

%% resample
t_new = linspace(data_A(1,1),data_A(end,1),samples);
dt = t_new(2)-t_new(1);
A_new = interp1(data_A(:,1),data_A(:,2),t_new);

%% Cut data
if ~isnan(myLims)
  if length(myLims)==2
    % Find indeces where to cut time series
    lowerInd = find( abs(t_new-myLims(1))<0.5*dt , 1 , 'first' );
    upperInd = find( abs(t_new-myLims(2))<0.5*dt , 1 , 'first' );
    % Now cut
    if lowerInd<upperInd
      t_new = t_new(lowerInd:upperInd);
      A_new = A_new(lowerInd:upperInd);
    else
      warning('There was a problem respecting the specified limits, so they were ignored!')
    end
  else
    warning('Limits must be specified as 2x1 vector like [-1 2]. Limits were ignored!')
  end
end

%% derive
% (a) Wavelet
% wavelet = 'bior1.3';
% dAdt_wavelet = derivative_cwt(A_new,wavelet,scaleSmoothing,dt,1);

% (b) Savitzky Golay Filter
degreePoly = 3;
derivative = 1;
frameSize = 2*scaleSmoothing + 1; % Must be odd!
dAdt_SavGol = -savitzkyGolayFilt(A_new,degreePoly,derivative,frameSize,[],2) / dt^derivative;

% (c) Moving average and central differences
% c.1 Moving average accounting for time delay
frameSize = 2*scaleSmoothing + 1;           % Window size
fDelay = (frameSize-1)/2;                   % Delay of the filter (indices)
tDelay = ( t_new(fDelay) - t_new(1) );          % Delay of the filter (time)
filterCoeff = ones(1, frameSize)/frameSize; % Filter coefficients
A_filtered = filter(filterCoeff, 1, A_new); % Moving average filter
% c.2 Cental differences
N_x1 = length(A_filtered);
Dx1 = spdiags([-ones(N_x1,1) ones(N_x1,1)],[-1 1],N_x1,N_x1);
dAdt_MACD = Dx1*A_filtered' / ( 2*dt );
% Cut off ends
t_new_noDelay = t_new(2:end-1) - tDelay;
A_filtered = A_filtered(2:end-1);
dAdt_MACD = dAdt_MACD(2:end-1);

%% Output
dAdt = dAdt_SavGol;

%% Norm time series if desired
if doNorm
  normFac = trapz(t_new,dAdt);
  dAdt = dAdt / normFac;
  t_new = t_new - t_new(1);
end

%% plot
if doPlot  
  figure
  subplot(numSubplots,1,1)
  plot(t_new,A_new); hold on
  plot(t_new_noDelay,A_filtered);
  xlim([t_new(1) t_new(end)])
  xlabel('t [s]')
  ylabel('A [m^2]')
  title(['Flame Surface (Case ',num2str(GFcase.solver.resume2case),')'])
  legend({'Original Data','Filtered Data (Moving Average)'});
  
%   subplot(numSubplots,1,2)
%   plot(t_new,dAdt_wavelet); hold on
%   plot(t_new,dAdt_SavGol);
%   plot(t_new_noDelay,dAdt_MACD);
%   xlim([t_new(1) t_new(end)])
%   xlabel('t [s]')
%   ylabel('dAdt [m^2/s]')
%   title(['1st Derivative Flame Surface (Smoothing Parameter: ',num2str(scaleSmoothing),')'])
%   legend({['Wavelet (',wavelet,')'] , 'Savitzky Golay Filter' , 'Moving Average and Central Differences' });
  
  if plotRef
    subplot(numSubplots,1,3)
    plot(data_V(:,1),data_V(:,2))
    xlabel('t [s]')
    ylabel('v_{ref}')
    title('Reference Velocity')
  end
end

end

