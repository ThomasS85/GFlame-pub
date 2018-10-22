%% This script is called by integrate_1DNLGeq( ) and parses varargin 

% How many time steps shall be written out at maximum?
ind = find(strcmpi(varargin,'nOut'),1);
if ~isempty(ind)
  % Use user specified
  nOut = varargin{ind+1};
else
  % use default: Write out 100 at maximum
  nOut = 100;
end

% Specify excitation signal (vector which varies between -1 and 1; multiplied with amplitude)
%    -> excSig = [ t_vec(:) , Signal_vec(:) ];
ind = find(strcmpi(varargin,'excSig'),1);
if ~isempty(ind)
  % Use user specified
  excSig = varargin{ind+1};
else
  % use default: Step input
  % SR
  excSig = [  0         , 0  ;...
    0.02      , 0  ;...
    0.02+eps  , 1  ;...
    t_max     , 1  ];
end

% Include flame-flow feedback  -> DILATATION (irrotational)
%  -> 1st input: Method: 'panel' or 'point'
%  -> 2st input: Set expansion ratio E (E>=1)
%  -> 3nd input: Set kernel Width of desingularized sources (in physical domain)
%  -> 4rd input: Sources/ panels per milimeter that should be equidistantly distributed
%  -> 5th input: mean-flame normal distance between G=0 and sources
ind = find(strcmpi(varargin,'FFF_I'),1);
if ~isempty(ind)
  % Yes, include feedback
  FFF.do_I = 1;
  FFF.method_I = varargin{ind+1};
  FFF.myE = varargin{ind+2};
  FFF.myDK_I = varargin{ind+3};
  FFF.SPpMM = varargin{ind+4};
  FFF.dispXi = varargin{ind+5};
  plotSources = 'sourceL';
else
  % use default: No feedback
  FFF.do_I = 0;
  FFF.myE = 1;
  plotSources = 'none';
end

% Include flame-flow feedback  -> Baroclinic tourque (solenoidal)
%  -> 1st input: Set expansion ratio E (E>=1)
%  -> 2nd input: Set kernel Width of desingularized vortexes (in physical domain)
ind = find(strcmpi(varargin,'FFF_S'),1);
if ~isempty(ind)
  % Yes, include feedback
  FFF.do_S = 1;
  myE_tmp = varargin{ind+1};
  FFF.myDK_S = varargin{ind+2};
else
  % use default: No feedback
  FFF.do_S = 0;
  myE_tmp = 1;
  FFF.myDK_S = 0;
end
% If expansion ratio has already been set by FFF_I -> check consistency!
if isfield(FFF,'myE')
  if (myE_tmp ~= FFF.myE) && FFF.do_S
    error('Inconsistent definition of expansion ratio (between dilatation and baroclinic torque!)!')
  end
else
  FFF.myE = myE_tmp;
end

% CFL number?
ind = find(strcmpi(varargin,'CFL'),1);
if ~isempty(ind)
  % Use user specified
  CFL = varargin{ind+1};
else
  % use default: 0.5
  CFL = 0.5;
end

% Should solution be visualized?
ind = find(strcmpi(varargin,'doPlot'),1);
if ~isempty(ind)
  % Yes
  doPlot = 1;
else
  % use default: No
  doPlot = 0;
end

% Should diffusion term of G-equation be considered?
%  1st input: markstein length
%  2nd input: Should curvature contribution on source/panel strengths be contributed?  1  or 0
ind = find(strcmpi(varargin,'doDiff'),1);
if ~isempty(ind)
  % Yes
  doDiff = 1;
  my_lM = varargin{ind+1}; % Markstein length
  curvContri_sL = varargin{ind+2};
else
  % use default: No
  doDiff = 0;
  curvContri_sL = 0;
  my_lM = [];
end

% Should all xi snapshots be returned by function?
ind = find(strcmpi(varargin,'returnAll'),1);
if ~isempty(ind)
  % Return all snapshots of xi
  doReturnAll = 1;
else
  % use default: No
  doReturnAll = 0;
end

% Forcing amplitude
ind = find(strcmpi(varargin,'amp'),1);
if ~isempty(ind)
  % Use user specified
  myAmp = varargin{ind+1};
else
  % use default: 10% bulk flow velocity
  myAmp = 0.1;
end

% Should gif be written (only if plotting is activated)
ind = find(strcmpi(varargin,'doGif'),1);
if ~isempty(ind)
  % Yes
  doGif = 1;
  filename4gif = varargin{ind+1};
else
  % use default: No
  doGif = 0;
end

% Solve around a certain mean value -> provide myMean struct!
%   myMean.flcoord      - coordinates of mean flame position in laboratory coordiantes L1 [x1 x2] (Nx2)
%   myMean.u_parallel   - flame-parallel mean flow component [u1_parallel] (Nx1)
%   myMean.vel          - alternatively specify vel struct with mean flow data (SC4F_toolbox)
%
%   Note: N should be sufficiently high! (similar to length(xi_0))
%
ind = find(strcmpi(varargin,'mean'),1);
if ~isempty(ind)
  % User specified
  myMean = varargin{ind+1};
else
  % use default: straight line with angle alpha
  if ~strcmp(p.flameType,'flat')
    x1L_tmp = linspace(0,p.H_flame,length(xi_0)).';
    x2L_tmp = p.R_flame - x1L_tmp*tan(p.alpha);
    % flame tangent vector (L1)
    myT = [  cos(p.alpha) ; -sin(p.alpha) ];
  else
    x1L_tmp = zeros(length(xi_0),1);
    x2L_tmp = linspace(p.R_i,0,length(xi_0)).';
    % flame tangent vector (L1)
    myT = [  0 ; 1 ];
  end
  myMean.flcoord = [ x1L_tmp , x2L_tmp ];
  % Compute mean flow component in flame parallel direction
  myMean.u_parallel = [ ones(length(xi_0),1)*p.u_1_centerIn , zeros(length(xi_0),1) ]*myT;
end


% Interpolate specified mean flame front coordinates + parallel velocity field to specified flame coordinates
[ myMean.flcoord , myMean.u_parallel ] = interp2Dcurve_2equidistantGrid_FPB( myMean.flcoord.' , length(xi_0) , 'fields' , {myMean.u_parallel.'} );
myMean.u_parallel = myMean.u_parallel{1};
% Make sure flcoord is line vector
[z,s] = size(myMean.flcoord);
if z>s; myMean.flcoord = myMean.flcoord.'; end;
% Make sure u_parallel is row vector
[z,s] = size(myMean.u_parallel);
if z<s; myMean.u_parallel = myMean.u_parallel.'; end;
% compute mean flame coordinates -> uniform spacing between 0 and mean flame length
%  Note: above flame coordinates were interpolated to equidistant grid, hence, x1F now is equidistant as well!
x1F = [ 0 , cumsum( sqrt( diff(myMean.flcoord(1,:)).^2 + diff(myMean.flcoord(2,:)).^2 ) ) ];
% Get spatial spacing
dx1F = x1F(2) - x1F(1);