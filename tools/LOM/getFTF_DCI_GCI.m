function [ slit_Lambda , slit_V , con_V , con_Lambda ] = getFTF_DCI_GCI( tau_r , tau_c , varargin )
%getFTF_DCI_GCI Computes analytic FR/ IR for three flame topologies for perfectly premixed flames
% for given convective (tau_c) and disturbance (tau_d) time scale.
%
% Used for Paper 
%   Steinbacher et al., "Consequences of Flame Geometry for the Acoustic Response of Premixed Flames", CNF 2018
%
% If you use Gaussian kernel, make shure the number of coefficients used for sampling the IR is set
%   sufficiently high!! 
%     -> line 211: numCoef = 70;
%
% by Thomas Steinbacher 05/2017

%% Parse varargin
% For V-flame: What is the non-dimensional rod radius r_tilde?
ind = find(strcmpi(varargin,'V'),1);
if ~isempty(ind)
  % User defined
  r_tilde = varargin{ind+1};
else
  % zero (no rod)
  r_tilde = 0;
end

% Use convective instead of incomp. convecvtive velocity model
ind = find(strcmpi(varargin,'convective'),1);
if ~isempty(ind)
  % Use convective velocity model
  useConv = 1;
else
  % Use convective incompressible velocity model
  useConv = 0;
end

% Use Gaussian kernel?
ind = find(strcmpi(varargin,'Gaussian'),1);
if ~isempty(ind)
  % Yes
  myKernel = 'Gaussian';
  myKernelSet = varargin{ind+1}; % Settings: [mySig0 D_f mydt]
else
  % No, use Dirac
  myKernel = 'Dirac';
  myKernelSet = zeros(1,3);
end

% Return only handle for one certain kind of topology (required for optimization)
%   slit slitV con conV
ind = find(strcmpi(varargin,'retOnly'),1);
if ~isempty(ind)
  % Return only one certain topo (as first output)
  doRetOnly1 = 1;
  retOnly = varargin{ind+1}; 
else
  % Return all topologies
  doRetOnly1 = 0;
  retOnly = 'all';
end

% Return only handle for one type of function (required for optimization)
%   IR FR
ind = find(strcmpi(varargin,'retFun'),1);
if ~isempty(ind)
  % Return only one certain function handle
  doRetFun = 1;
  retFun = varargin{ind+1}; 
else
  % Return all handles (default)
  doRetFun = 0;
  retFun = 'all';
end

% Directly evaluate requested function handle using provided vector (times or ang frequency)?
%   Only available if 'retFun' and 'retOnly' is chosen, as well!
ind = find(strcmpi(varargin,'call'),1);
if ~isempty(ind)
  % Yes
  doCall = 1;
  call_vec = varargin{ind+1}; 
else
  % No
  doCall = 0;
  call_vec = [];
end

%% Prepare kernel
% Choose kernel function
if strcmpi(myKernel,'Dirac')
  % Dirac kernel
  myg = @(y,tau) y==0;
  myG = @(y,tau) myHeaviside2(y);
  myGamma = @(y,tau) myHeaviside2(y) .* y;
  % myKernelFreq = @(w) 1;
  % myDiracFreq = @(w) 1;
elseif strcmpi(myKernel,'Gaussian')
  % Settings Gaussian Kernel
  mySig0 = myKernelSet(1);
  D_f = myKernelSet(2);
  mydt = myKernelSet(3);
  % Gaussian Kernel
  mySig = @(tau) mySig0 + 2*D_f*tau;
  myg = @(y,tau) exp(- (y-mydt).^2 ./ (2*mySig(tau)) ) ./ sqrt(2*pi*mySig(tau));
  myG = @(y,tau) erf( (y-mydt) ./ sqrt(2*mySig(tau)) ) / 2 ;
  myGamma = @(y,tau) (y-mydt) .* erf( (y-mydt) ./ sqrt(2*mySig(tau)) ) / 2 +...
    sqrt( mySig(tau)/(2*pi) ) .* exp( -(y-mydt).^2 ./ (2*mySig(tau))  );
  % myKernelFreq = @(w) exp( 1i*w*(mySig0/2-mydt) ) .* erfz( (mySig0*1i*w-2*mydt)/sqrt(2*mySig0) );
  % myDiracFreq = @(w,t0) exp( 1i*w*(mySig0/2-t0) ) .* erfz( (mySig0*1i*w-2*t0)/sqrt(2*mySig0) );
end

%% Compute IRs
delta_tau = tau_r - tau_c;

% General
fac1 = 1 / (tau_r-tau_c);
fac2 = tau_c / (tau_r-tau_c) * (1-useConv);

% For V flames
fac_r0 = (1 - r_tilde) / (1 + r_tilde);
fac_r1 = r_tilde / (1 - r_tilde);
fac_r2 = 1 / (1 - r_tilde);

% Slit Lambda
slit_Lambda.flanks = @(t) fac1 * ( 1 + fac2 ) * ( myG(t-tau_c,tau_c) - myG(t-tau_r,tau_r) );
slit_Lambda.IR1 = @(t) -fac1 * tau_c * myg(t-tau_r,tau_r) * (1-useConv); 
slit_Lambda.h_u = @(t) slit_Lambda.flanks(t) + slit_Lambda.IR1(t);
% Secondary deflection
slit_Lambda.flanks_2 = @(t) -fac1 * fac2 * ( myG(t-tau_c-delta_tau,tau_c-delta_tau) - myG(t-tau_r-delta_tau,tau_r-delta_tau) );
slit_Lambda.IR1_2 = @(t) fac1 * tau_c * myg(t-tau_r-delta_tau,tau_r-delta_tau) * (1-useConv); 
slit_Lambda.h_u_2 = @(t) slit_Lambda.h_u(t) + slit_Lambda.flanks_2(t) + slit_Lambda.IR1_2(t);

% Slit V
slit_V.flanks = @(t) fac1 * ( 1 - fac2 ) * ( myG(t-tau_c,tau_c) - myG(t-tau_r,tau_r) );
slit_V.IR1 = @(t) -fac1 * tau_c * fac_r1 * myg(t-tau_r,tau_r) * (1-useConv);
slit_V.IR2 = @(t) fac1 * tau_c * fac_r2 * myg(t-tau_c,tau_c) * (1-useConv);
slit_V.h_u = @(t) slit_V.flanks(t) + slit_V.IR1(t) + slit_V.IR2(t);

% Conical Lambda
con_Lambda.flanks = @(t) 2*fac1 * ( 1 + fac2/2 ) * ( ...
  ( myGamma(t-tau_r,tau_r) - myGamma(t,0) ) / tau_r - ( myGamma(t-tau_c,tau_c) - myGamma(t,0) ) / tau_c );
con_Lambda.IR2 = @(t) 2*fac1 * 0.5*( myG(t,0) + ( myGamma(t-tau_c,tau_c) - myGamma(t,0) ) / tau_c ) * (1-useConv);
con_Lambda.IR1 = @(t) 2*fac1 * 0.5*tau_c/tau_r * ( myG(t-tau_r,tau_r) - myG(t,0) ) * (1-useConv);
con_Lambda.h_u = @(t) con_Lambda.flanks(t) + con_Lambda.IR1(t) + con_Lambda.IR2(t);

% Conical V
con_V.flanks = @(t) 2*fac1 * fac_r0 * ( 1 - fac2/2 ) * ( ...
  fac_r2*( myG(t-tau_c,tau_c) - myG(t-tau_r,tau_r) ) ...
  -(( myGamma(t-tau_r,tau_r) - myGamma(t,0) ) / tau_r - ( myGamma(t-tau_c,tau_c) - myGamma(t,0) ) / tau_c) );
con_V.IR2 = @(t)  2*fac1 * fac_r0 * (...
  fac_r2 * tau_c/2 * fac_r2 * myg(t-tau_c,tau_c) ...
  - 0.5*( fac_r1*myG(t,0) - fac_r2*myG(t-tau_c,tau_c) - ( myGamma(t-tau_c,tau_c) - myGamma(t,0) ) / tau_c ) )...
   * (1-useConv);
con_V.IR1= @(t) 2*fac1 * fac_r0 * (...
  -fac_r2 * tau_c/2 * fac_r1 * myg(t-tau_r,tau_r) ...
  - fac_r1/2 * tau_c/tau_r * ( myG(t-tau_r,tau_r) - myG(t,0) ) ) * (1-useConv);
con_V.imp_d = @(t) 2*fac1 * fac_r0 * ( fac_r2 * tau_c/2 * fac_r2 * myg(t-tau_c,tau_c) ); % Impulse only
con_V.imp_c = @(t) 2*fac1 * fac_r0 * (-fac_r2 * tau_c/2 * fac_r1 * myg(t-tau_r,tau_r) ); % Impulse only
con_V.h_u = @(t) con_V.flanks(t) + con_V.IR1(t) + con_V.IR2(t);



%% Compute FRs
if strcmpi(myKernel,'Dirac')
  % Definitions
  myg_FRF = @(w,tau) exp(-1i*w*tau);
  myG_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w);
  myGamma_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w).^2;
  
  %%%%%%%%%%%%
  % Just copy+paste IR and replace 'kernels' with 'kernels_FRF' and 't-...' by 'w'!!
  %%%%%%%%%%%%
  
  % Slit Lambda
  slit_Lambda.flanks_FRF = @(w) fac1 * ( 1 + fac2 ) * ( myG_FRF(w,tau_c) - myG_FRF(w,tau_r) );
  slit_Lambda.IR1_FRF = @(w) -fac1 * tau_c * myg_FRF(w,tau_r) * (1-useConv);
  slit_Lambda.F_u = @(w) slit_Lambda.flanks_FRF(w) + slit_Lambda.IR1_FRF(w);
  slit_Lambda.F_uLP = @(w) slit_Lambda.flanks_FRF(w) + slit_Lambda.IR1_FRF(w) ./ (1+1i*w*tau_r*1e-0);  
  % Secondary deflection
  slit_Lambda.flanks_FRF_2 = @(w) -fac1 * fac2 * ( myG_FRF(w,tau_c-delta_tau) - myG_FRF(w,tau_r-delta_tau) );
  slit_Lambda.IR1_FRF_2 = @(w) fac1 * tau_c * myg_FRF(w,tau_r-delta_tau) * (1-useConv);
  slit_Lambda.F_u_2 = @(w) slit_Lambda.F_u(w) + slit_Lambda.flanks_FRF_2(w) + slit_Lambda.IR1_FRF_2(w);
  
  % Slit V
  slit_V.flanks_FRF = @(w) fac1 * ( 1 - fac2 ) * ( myG_FRF(w,tau_c) - myG_FRF(w,tau_r) );
  slit_V.IR1_FRF = @(w) -fac1 * tau_c * fac_r1 * myg_FRF(w,tau_r) * (1-useConv);
  slit_V.IR2_FRF = @(w) fac1 * tau_c * fac_r2 * myg_FRF(w,tau_c) * (1-useConv);
  slit_V.F_u = @(w) slit_V.flanks_FRF(w) + slit_V.IR1_FRF(w) + slit_V.IR2_FRF(w);
  
  % Conical Lambda
  con_Lambda.flanks_FRF = @(w) 2*fac1 * ( 1 + fac2/2 ) * ( ...
    ( myGamma_FRF(w,tau_r) - myGamma_FRF(w,0) ) / tau_r - ( myGamma_FRF(w,tau_c) - myGamma_FRF(w,0) ) / tau_c );
  con_Lambda.IR2_FRF = @(w) 2*fac1 * 0.5*( myG_FRF(w,0) + ( myGamma_FRF(w,tau_c) - myGamma_FRF(w,0) ) / tau_c )...
     * (1-useConv);
  con_Lambda.IR1_FRF = @(w) 2*fac1 * 0.5*tau_c/tau_r * ( myG_FRF(w,tau_r) - myG_FRF(w,0) ) * (1-useConv);
  con_Lambda.F_u = @(w) con_Lambda.flanks_FRF(w) + con_Lambda.IR1_FRF(w) + con_Lambda.IR2_FRF(w);
  
  % Conical V
  con_V.flanks_FRF = @(w) 2*fac1 * fac_r0 * ( 1 - fac2/2 ) * ( ...
    fac_r2*( myG_FRF(w,tau_c) - myG_FRF(w,tau_r) ) ...
    -(( myGamma_FRF(w,tau_r) - myGamma_FRF(w,0) ) / tau_r - ( myGamma_FRF(w,tau_c) - myGamma_FRF(w,0) ) / tau_c) );
  con_V.IR2_FRF = @(w)  2*fac1 * fac_r0 * (...
    fac_r2 * tau_c/2 * fac_r2 * myg_FRF(w,tau_c) ...
    - 0.5*( fac_r1*myG_FRF(w,0) - fac_r2*myG_FRF(w,tau_c) - ( myGamma_FRF(w,tau_c) - myGamma_FRF(w,0) ) / tau_c ) )...
     * (1-useConv);
  con_V.IR1_FRF = @(w) 2*fac1 * fac_r0 * (...
    -fac_r2 * tau_c/2 * fac_r1 * myg_FRF(w,tau_r) ...
    - fac_r1/2 * tau_c/tau_r * ( myG_FRF(w,tau_r) - myG_FRF(w,0) ) ) * (1-useConv);
  con_V.F_u = @(w) con_V.flanks_FRF(w) + con_V.IR1_FRF(w) + con_V.IR2_FRF(w);
  
else
  % Numerical FR
  numCoef = 70;
  % Sample IRs
  t_IR = linspace(0,tau_r*1.5,numCoef);
  IR_con_lambda = con_Lambda.h_u(t_IR);
  IR_con_V = con_V.h_u(t_IR);
  IR_slit_lambda = slit_Lambda.h_u(t_IR);
  IR_slit_V = slit_V.h_u(t_IR);
  
  % Evaluate FR numerically
  [ con_Lambda.m ] = FTFfromIR( IR_con_lambda , t_IR );
  [ con_V.m ] = FTFfromIR( IR_con_V , t_IR );
  [ slit_Lambda.m ] = FTFfromIR( IR_slit_lambda , t_IR );
  [ slit_V.m ] = FTFfromIR( IR_slit_V , t_IR );
  
  % Provide function handle by interpolation
  con_Lambda.F_u = @(w) interp1(con_Lambda.m.omega_vec,real(con_Lambda.m.comp),w) + ...
    1i*interp1(con_Lambda.m.omega_vec,imag(con_Lambda.m.comp),w);
  con_V.F_u = @(w) interp1(con_V.m.omega_vec,real(con_V.m.comp),w) + ...
    1i*interp1(con_V.m.omega_vec,imag(con_V.m.comp),w);
  slit_Lambda.F_u = @(w) interp1(slit_Lambda.m.omega_vec,real(slit_Lambda.m.comp),w) + ...
    1i*interp1(slit_Lambda.m.omega_vec,imag(slit_Lambda.m.comp),w);
  slit_V.F_u = @(w) interp1(slit_V.m.omega_vec,real(slit_V.m.comp),w) + ...
    1i*interp1(slit_V.m.omega_vec,imag(slit_V.m.comp),w);

end

% Should only one certain function handle be returned?
if doRetFun
  switch retFun
    case 'IR'
      % Return IR
      con_Lambda = con_Lambda.h_u;
      con_V = con_V.h_u;
      slit_Lambda = slit_Lambda.h_u;
      slit_V = slit_V.h_u;
    case 'FR'
      % Return FR
      con_Lambda = con_Lambda.F_u;
      con_V = con_V.F_u;
      slit_Lambda = slit_Lambda.F_u;
      slit_V = slit_V.F_u;
    otherwise
      error('A certain function handle was requested but the specified output name does not exist!')
  end
end

% Should only a certain topology be returned?
if doRetOnly1
  switch retOnly
    case 'slit'
      % Do nothing
    case 'slitV'
      slit_Lambda = slit_V;
    case 'con'
      slit_Lambda = con_Lambda;
    case 'conV'
      slit_Lambda = con_V;
    otherwise
      error('Only one output was requested but the specified output name does not exist!')
  end
end

% Evaluate function handle with provided vector
if doCall && doRetOnly1 && doRetFun
  slit_Lambda = slit_Lambda(call_vec);
end

end

