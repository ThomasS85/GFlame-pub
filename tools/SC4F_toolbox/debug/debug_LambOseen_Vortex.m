clear
close all

addPath2Bib('GFlame-pub')

%% Setup
% Set default values for parameters:    Steinbacher_2D_Cr66    or     Steinbacher_2D
myFlame = 'Steinbacher_2D_Cr66';
[ p ] = setUpPredefinedFlame( myFlame ); 

% Position Vortex
xi_s = 0.05*1i;

%% Compute Velocity field
% Settings
myAmp = 0.1;
dt = 1e-4;
excSig = [ 0 , 2/dt ; dt , 0 ; 1 , 0 ];
Gamma = p.u_1_bulkFeed * dt/2*myAmp*p.u_1_bulkFeed* ( excSig(1,2) + excSig(2,2) );
% radial axis
xi = 1i*linspace(imag(xi_s)-0.2,imag(xi_s)+0.2,1000);

% Compute point Vortex
[ u_vort , dOmega_dxi ] = evalVel_vortex( xi , xi_s , Gamma , p , 'xi' );
u_vort_xi = abs(dOmega_dxi);

% Compute Lamb Oseen Vortex
vortParams = {0.01};
[ u_vort_LO , dOmega_dxi_LO ] = evalVel_vortex( xi , xi_s , Gamma , p , 'xi' , 'lamb' , vortParams );
u_vort_xi_LO = abs(dOmega_dxi_LO);


%% Plot
myXlim = [imag(xi_s)-0.1 imag(xi_s)+0.1];
myYlim = [0 10];

figure('Position',[0 50 700 700]); hold on;
plot(imag(xi),u_vort_xi,'b-')
plot(imag(xi),u_vort_xi_LO,'g-.')
plot(imag(xi_s),0,'bo')
plot(-imag(xi_s),0,'ro')
xlim(myXlim);ylim(myYlim)
legend('Point Vortex','Lamb-Oseen Vortex')

figure('Position',[520 50 500 500]); hold on;
plot(imag(xi),u_vort_xi_LO,'-')
plot(imag(xi_s),0,'o')
xlim(myXlim);ylim(myYlim)