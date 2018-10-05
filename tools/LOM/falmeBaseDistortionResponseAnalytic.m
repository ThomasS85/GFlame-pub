function [  slit , con_V , con_Lambda ] = falmeBaseDistortionResponseAnalytic( tau_c , T )
%FALMEBASEDISTORTIONRESPONSEANALYTIC Function computes response in terms of heat release to a normalized perturbation of
% velocity at the flame base assuming a rectangular profile


%% Kernels
% Dirac kernel for IR
myg = @(y,tau) (y-tau)==0;
myG = @(y,tau) myHeaviside2(y-tau);
myGamma = @(y,tau) myHeaviside2(y-tau) .* (y-tau);

% Kernels for FR
myg_FRF = @(w,tau) exp(-1i*w*tau);
myG_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w);
myGamma_FRF = @(w,tau) exp(-1i*w*tau) ./ (1i*w).^2;

%% Define FTFs with general Kernel
% General definitions
tau_d = tau_c * T;
% (1) Slit
slit.gen = @(myVar,g,G,Gamma) 1/tau_d * ( G(myVar,tau_c-tau_d) - G(myVar,tau_c) );

% (2) Conical
con_Lambda.gen = @(myVar,g,G,Gamma) 1/(tau_c-tau_d/2) * ( G(myVar,0) - Gamma(myVar,tau_c-tau_d)/tau_d...
    + Gamma(myVar,tau_c)/tau_d );

% (3) Conical V
con_V.gen =  @(myVar,g,G,Gamma) 2*slit.gen(myVar,g,G,Gamma) - con_Lambda.gen(myVar,g,G,Gamma);

%% Compute IRs
% (1) Slit
slit.h_u = @(t) slit.gen(t,myg,myG,myGamma); 
% (2) Conical
con_Lambda.h_u = @(t) con_Lambda.gen(t,myg,myG,myGamma); 
% (3) Conuical V
con_V.h_u = @(t) con_V.gen(t,myg,myG,myGamma); 

%% Analytical evaluation FR from IR
% (1) Slit
slit.F_u = @(w) slit.gen(w,myg_FRF,myG_FRF,myGamma_FRF);
% (2) Conical
con_Lambda.F_u = @(w) con_Lambda.gen(w,myg_FRF,myG_FRF,myGamma_FRF);
% (3) Conical V
con_V.F_u = @(w) con_V.gen(w,myg_FRF,myG_FRF,myGamma_FRF);

end

