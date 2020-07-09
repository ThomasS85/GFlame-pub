function [ slit , con_V , con_Lambda  ] = flameBaseResponseAnalytic( tau_c , varargin )
%FLAMEBASERESPONSEANALYTIC Function computes response in terms of heat release to a normalized perturbation of
%the flame base flame position xi/(s_L*tau_c) for 3 flame topologies
%
% By Thomas Steinbacher Aug 2017


%% Parse varargin
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

% For V-flame: What is the non-dimensional rod radius r_tilde?
ind = find(strcmpi(varargin,'V'),1);
if ~isempty(ind)
  % User defined
  r_tilde = varargin{ind+1};
else
  % zero (no rod)
  r_tilde = 0.05;
end


%% Prepare kernel
% Choose kernel function
if strcmpi(myKernel,'Dirac')
  % Dirac kernel for IR
  myg = @(y,tau) (y-tau)==0;
  myG = @(y,tau) myHeaviside2(y-tau);
  myGamma = @(y,tau) myHeaviside2(y-tau) .* (y-tau);
  
  % Kernels for FR
  myg_FRF = @(w,tau) exp(-1i*w*tau);
  myG_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w);
  myGamma_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w).^2;
  
elseif strcmpi(myKernel,'Gaussian')
  % Settings Gaussian Kernel
  mySig0 = myKernelSet(1);
  D_f = myKernelSet(2);
  mydt = myKernelSet(3);
  % Gaussian Kernel
  mySig = @(tau) mySig0 + 2*D_f*tau;
  myg = @(y,tau) exp(- (y-mydt-tau).^2 ./ (2*mySig(tau)) ) ./ sqrt(2*pi*mySig(tau));
  myG = @(y,tau) erf( (y-mydt-tau) ./ sqrt(2*mySig(tau)) ) / 2 ;
  myGamma = @(y,tau) (y-mydt-tau) .* erf( (y-mydt-tau) ./ sqrt(2*mySig(tau)) ) / 2 +...
    sqrt( mySig(tau)/(2*pi) ) .* exp( -(y-mydt-tau).^2 ./ (2*mySig(tau))  );

end


%% Define FTFs with general Kernel
% (1) Slit
slit.gen = @(myVar,g,G,Gamma) g(myVar,tau_c);

% (2) Conical
% con_Lambda.gen = @(myVar,g,G,Gamma) 2/tau_c * ( G(myVar,0) - G(myVar,tau_c) )*1.05 - g(myVar,tau_c)*0.1;
con_Lambda.gen = @(myVar,g,G,Gamma) 1/tau_c * ( G(myVar,0) - G(myVar,tau_c) );

% (3) Conuical V
con_V.gen =  @(myVar,g,G,Gamma) 1/r_tilde * ( g(myVar,tau_c) - (1-r_tilde)/tau_c * ( G(myVar,0) - G(myVar,tau_c) ) );
% con_V.gen =  @(myVar,g,G,Gamma) preConV * (2*g(myVar,tau_cV) - 1/tau_cV * ( G(myVar,0) - G(myVar,tau_cV) ) );


%% Compute IRs
% (1) Slit
slit.h_xi0 = @(t) slit.gen(t,myg,myG,myGamma); 
% (2) Conical
con_Lambda.h_xi0 = @(t) con_Lambda.gen(t,myg,myG,myGamma); 
% (3) Conuical V
con_V.h_xi0 = @(t) con_V.gen(t,myg,myG,myGamma); 

% Return weight of impulses if kernel is Dirac
if strcmpi(myKernel,'Dirac')
  slit.imp_c = @(t) myg(t,tau_c);
  con_V.imp_c = @(t) 1/r_tilde*myg(t,tau_c);
end

%% Compute FRs
if strcmpi(myKernel,'Dirac')
  % Analytical evaluation FR from IR
  % (1) Slit
  slit.F_xi0 = @(w) slit.gen(w,myg_FRF,myG_FRF,myGamma_FRF); 
  % (2) Conical
  con_Lambda.F_xi0 = @(w) con_Lambda.gen(w,myg_FRF,myG_FRF,myGamma_FRF); 
  % (3) Conical V
  con_V.F_xi0 = @(w) con_V.gen(w,myg_FRF,myG_FRF,myGamma_FRF); 
  
else
  % Numerical evaluation FR
  numCoef = 70;
  % Sample IRs
  t_IR = linspace(0,tau_c*1.5,numCoef);
  IR_con_lambda = con_Lambda.h_xi0(t_IR);
  IR_con_V = con_V.h_xi0(t_IR);
  IR_slit = slit.h_xi0(t_IR);
  
  % Evaluate FR numerically
  [ con_Lambda.m ] = FTFfromIR( IR_con_lambda , t_IR , 'noPhase0' );
  [ con_V.m ] = FTFfromIR( IR_con_V , t_IR , 'noPhase0' );
  [ slit.m ] = FTFfromIR( IR_slit , t_IR , 'noPhase0' );
  
  % Provide function handle by interpolation
  con_Lambda.F_xi0 = @(w) interp1(con_Lambda.m.omega_vec,real(con_Lambda.m.comp),w) + ...
    1i*interp1(con_Lambda.m.omega_vec,imag(con_Lambda.m.comp),w);
  con_V.F_xi0 = @(w) interp1(con_V.m.omega_vec,real(con_V.m.comp),w) + ...
    1i*interp1(con_V.m.omega_vec,imag(con_V.m.comp),w);
  slit.F_xi0 = @(w) interp1(slit.m.omega_vec,real(slit.m.comp),w) + ...
    1i*interp1(slit.m.omega_vec,imag(slit.m.comp),w);

end


end

