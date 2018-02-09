function [out] = identifyGFLAMEOutputData( GFcase , varargin )
%PLOTSURFACEDATA Identifies flame transfer function (FTF) from broad band excitation. Takes most settings from
% GFcase struct. 
%
% Function requires siEstimateFIR() and siIdentifyFlame() from TFD tools to work properly!
% 
% Optional Inputs:
%   case            : Case number (default from GFcase)
%   root            : Case root folder (path) (default from GFcase)
%   timeDelay       : Force time delay [s]  (default from GFcase.vel.tauBB)
%   timeDelayCoeff  : Force time delay [# of coefficients] (default from GFcase.vel.tauBB) 
%                     Default way for time delays!
%   maxFreq         : Maximum Frequency of identified data
%   t0              : First time step of time series (default from GFcase, t_transient0)
%   tend            : Last time step of time series (default is 'last')
%   Ts              : sampling time for iddata object generation (default is mean Ts of written out surface)
%   Name            : Name of the case for plotting
%   pSI             : Fraction of avaiable data which is used for identification. The rest is used as test
%                     data. (Default is 0.8)
%   compOrders      : Vector which contains the number of impulse coefficients
%   idModel         : Defines model for identification (Default is FIR from impulseest())
%   doPlot          : If specified, output is plotted (Default from GFcase.soler.doPlot)
%
% Outputs:
%   out       : Struct which contains...
%     - casename    : Cell array with name strings 
%     - models      : Cell array with identified transfer functions (idtf objects)
%     - maxFreq     : Maximum frequency of identified models
%     - Ts          : sampling rate of identified models
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% Parse varargin
% Should case be change?
ind = find(strcmp(varargin,'case'),1);
if ~isempty(ind)
  % Change run output folder
  GFcase.solver.resume2case = varargin{ind+1};
else
  % Load case specified by case number
  GFcase.solver.resume2case = GFcase.p.caseNumber;
end

% Should root folder be changed? be change?
ind = find(strcmp(varargin,'root'),1);
if ~isempty(ind)
  % Change root folder
  GFcase.p.caseRootDir = varargin{ind+1};
else
  % Leave everything as it is
end

% Set time delay in s
ind = find(strcmpi(varargin,'timeDelay'),1);
if ~isempty(ind)
  % Use user specified time delay
  tau = varargin{ind+1};
  timeDelayMode = 'timeDelay';
else
  % use default
  tau = GFcase.vel.tauBB;
  timeDelayMode = 'timeDelayCoeff';
end

% Set time delay in # of coefficients
ind = find(strcmpi(varargin,'timeDelayCoeff'),1);
if ~isempty(ind)
  % Use user specified number of time delay coefficients
  nk = varargin{ind+1};
  timeDelayMode = 'timeDelayCoeff';
else
  % Default value nk
  nk = -1;
end


% Maximum Frequency
ind = find(strcmp(varargin,'maxFreq'),1);
if ~isempty(ind)
  maxFreq = varargin{ind+1};
else
  % default value
  maxFreq = GFcase.vel.f_maxBB;
end


% First time step of time series ('first' for start at first time step)
ind = find(strcmp(varargin,'t0'),1);
if ~isempty(ind)
  t0 = varargin{ind+1};
else
  % default value
  t0 = GFcase.vel.t_transient0;
end

% First time step of time series
ind = find(strcmp(varargin,'tend'),1);
if ~isempty(ind)
  tend = varargin{ind+1};
else
  % default value
  tend = 'last';
end

% Sampling Time
ind = find(strcmp(varargin,'Ts'),1);
if ~isempty(ind)
  Ts = varargin{ind+1};
else
  % default value
  Ts = 'mean';
end

% Name
ind = find(strcmp(varargin,'Name'),1);
if ~isempty(ind)
  Name = varargin{ind+1};
else
  % default value
  Name = '2D-Gequation Broad';
end

% How Much of the time series should be used for identification? The rest
% is used for model checking (if set to 1 no check is performed!)
ind = find(strcmpi(varargin,'pSI'),1);
if ~isempty(ind)
  % Use user specified limits
  portionSIData = min( abs( varargin{ind+1} ) , 1 );
else
  % use default
  portionSIData = 0.8;
end

% How many models with different  orders of their FIR should be
% compared? Default is 1 of order 30;
ind = find(strcmpi(varargin,'compOrders'),1);
if ~isempty(ind)
  % Use user specified orders (as 1x2 vector!)
  compOrders = varargin{ind+1};
else
  % use default
  compOrders = GFcase.vel.FIR_ord;
end

% What kind of model should be used for identification? (Default: FIR)
%     FIR , ARX , ARXregul or BJ
% Second input argument should be model parameter:
%   ARX:  1st: Order of the polynom A (autoregressive part) [ nA ]
%   BJ:   1st: Order of the polynom C (noise model)
%         2nd: Order of Polynoms D
%         3rd: Order of Polynom F                           [ nB nD nF ]
ind = find(strcmpi(varargin,'idModel'),1);
if ~isempty(ind)
  % Use user specified model for identification
  idModel = varargin{ind+1};
  idModelParam = varargin{ind+2};
else
  % use default
  idModel = 'FIR';
  idModelParam = [];
end

% Should results be plotted? Default: no
ind = find(strcmpi(varargin,'doPlot'),1);
if ~isempty(ind)
  % use default: Do plot
  plotData = 1;
else
  % Take settings from solver
  plotData = GFcase.solver.doPlot;
  
end

% Do you want to plot reference data (e.g. sweep or digitalized from a paper)? Specify path to *.mat file
% as string which contains an idfrd object 'h' with the data ( it will be plotted using bodeplot(h) )
ind = find(strcmpi(varargin,'refBode'),1);
if ~isempty(ind)
  % Use user specified path to idfrd object
  bodeRef = 1;
  bodeRefData = varargin{ind+1};
else
  % No data (default) 
  bodeRef = 0;
  bodeRefData = [];
end


%% Load and prepare data
% Change dir to case root
change2caseRootDir( GFcase )

% define where to load data from
GFcase.p.run_output_folder = [GFcase.p.caseFolder,'/',num2str( GFcase.solver.resume2case )];

% Read surface data from file
f1 = fopen([GFcase.p.run_output_folder,'/Surface.out'],'r');
data_A = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
fclose(f1);

% Read reference velocity data from file
f1 = fopen([GFcase.p.run_output_folder,'/V_ref.out'],'r');
data_V = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
fclose(f1);

% Cut data (t0 to tend)
if isnumeric(tend) && isnumeric(t0)
  data_A = data_A( data_A(:,1)>=t0 & data_A(:,1)<=tend , : );
  data_V = data_V( data_V(:,1)>=t0 & data_V(:,1)<=tend , : );
elseif strcmpi(tend,'last') && isnumeric(t0)
  data_A = data_A( data_A(:,1)>=t0 , : );
  data_V = data_V( data_V(:,1)>=t0 , : );
elseif strcmpi(t0,'first') && isnumeric(tend)
  data_A = data_A( data_A(:,1)<=tend , : );
  data_V = data_V( data_V(:,1)<=tend , : );
end


% Calculate constant sampling time
if strcmpi(Ts,'mean')
  Ts = mean( diff(data_A(:,1)) );
end

% Interpolate data_A and data_V to same t_vec (of data_A)
data_Vtmp = interp1(data_V(:,1),data_V(:,2),data_A(:,1),'linear','extrap');
data_V = []; data_V(:,1) = data_A(:,1); data_V(:,2) = data_Vtmp;

% Uncomment only if too many time steps were written out:
% data_V = data_V(1:10:end,:);

% Generate iddata object
data = iddata( data_A(:,2) , data_V(:,2) , [] ,'SamplingInstants',data_A(:,1),...
  'Name','GFLAME','ExperimentName',Name,...
  'TimeUnit','s','InputUnit','m/s','OutputUnit','W',...
  'InputName','uref','OutputName','dQ');
% Resample data to constant Ts
data = siInterpolate(data,Ts);
t_vec = data.SamplingInstants;




%% Estimate IR
% Settings
p.filterType = idModel;
p.compOrders = compOrders;
p.modelParam = idModelParam;

% If time delay coefficients are set, tau will not be considered!
if strcmp(timeDelayMode,'timeDelay')
  % Time delay [s]
  if tau < 0
    if strcmp(GFcase.p.flameType,'slit')
      % if tau<0 take convective time delay
      p.tau = GFcase.p.tau_c * 0.9;
    else
      % No time delay for conical flames
      p.nk = 0;
    end
    
  elseif tau == 0
    % Switch to coefficient mode with 0 coefficients (retain confidence intervals)
    p.nk = 0;
  else
    % otherwise take user specified time delay
    p.tau = tau;
  end
  
elseif strcmp(timeDelayMode,'timeDelayCoeff')
  % Time delay coefficients
  if nk < 0
    if strcmp(GFcase.p.flameType,'slit')
      % if nk<0 take convective time delay and multiply by Nyquist frequency
      p.nk = floor( GFcase.p.tau_c * 0.85 * maxFreq*2 );
    else
      % No time delay for conical flames
      p.nk = 0;
    end
    
  else
    % Take user input
    p.nk = nk;
  end
  
end


% Other data
p.filePath = data;
p.maxFrequ = maxFreq;
p.t0 = t_vec(1);

p.pSI = portionSIData;
p.referenceVel = 'uref';
p.plotTitle = ['GFLAME ',GFcase.p.flameType,' (',idModel,')'];
p.outputFolder = GFcase.p.run_output_folder;

if plotData
  p.doPlot = 'doPlot';
  % Load reference data if desired
  if bodeRef
    load(bodeRefData)
    p.refBode = 'refBode';
    p.hExp = h;
  end
else
  
  p.doPlot = 'doNotPlot';
  
end

% p.refBode = 'noRefBode';
% p.hExp = [];
% p.refModel = 'noRefModel';
% p.refModelData = [];
% p.refModelName = [];

studyName = ['GFLAME_',GFcase.p.flameType,'_',idModel];

% Identify TF
[ out ] = siIdentifyFlame( studyName , p );


%% Save data to disk
save([GFcase.p.run_output_folder,filesep,'response_broad.mat'],'data')


end

