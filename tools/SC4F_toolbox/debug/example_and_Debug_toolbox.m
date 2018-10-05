clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCflameTools  (SCFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Requires 1DLinROM setup functions
addPath2Bib('GFLAME-pub')


%% Set parameters 
% Set default values for parameters
[ p ] = setUpPredefinedFlame( 'Steinbacher_2D' );


%% Be sure to use only coordinates system L2 for mapping!
% Coordinate System L1:
%    ^
%    |  |
% x2 |-- o                     x_L(1)
%    |      o                    :
%    |         o                 :
%    |            o              :
%  0 |__ . __ . __ . __ . __   x_L(end)
%    |
%    ---------------------> 
%       0                x1
%
% % Coordinate System L2:
%    ^
% x2 |__ . __ . __ . __ . __  x_L(1)
%    |          o               :
%    |       o                  :
%    |--- o                     :
%  0 |   |_________________    x_L(end)
%    |
%    ---------------------> 
%        0               x1


%% Define points in L1 system and check mapping
% Point in L1 system
x_L1 = [ 2*p.H_flame + 1i*p.R_i , 1.5*p.H_flame + 1i*p.R_i , -p.H_flame + 1i*0.5*p.R_i ];
% Map to image domain
[ xi_SC ] = SCmapInv_SCFT( x_L1 , p , 'L1' );
% Map back to physical domain
x_L1_2 = SCmap_SCFT( xi_SC , p , 'L1' );

disp(' ')
disp('Check mapping:')
disp(['x original: ',num2str(x_L1)])
disp(['x mapped  : ',num2str(x_L1_2)])

if norm(x_L1-x_L1_2)>1e-10
  error('Mapping did not work properly!')
end


%% Check velocity
% Points far upstream of flame (L1)
x_up_L1 = [ -p.H_flame + 1i*p.R_i , -p.H_flame + 1i*0.5*p.R_i , -p.H_flame ];
% Points far downstream of flame (L1)
x_down_L1 = [ p.H_flame*1.5 , p.H_flame*1.5 + 1i*0.5*p.R_a , p.H_flame*1.5 + 1i*p.R_a ];

% Map to image domain
[ xi_up ] = SCmapInv_SCFT( x_up_L1 , p , 'L1' ); 
[ xi_down ] = SCmapInv_SCFT( x_down_L1 , p , 'L1' );
% Map back to physical domain
x_up_L1_2 = SCmap_SCFT(xi_up, p , 'L1' );
x_down_L1_2 = SCmap_SCFT(xi_down, p , 'L1' );

disp(' ')
disp('Checking mapping for velocity evluation:')
disp(['x_up original: ',num2str(x_up_L1)])
disp(['x_up mapped  : ',num2str(x_up_L1_2)])
disp(' ')
disp(['x_down original: ',num2str(x_down_L1)])
disp(['x_down mapped  : ',num2str(x_down_L1_2)])


% Get mapping infos
[ s ] = return_SCmap_SCFT( p );

% Compute velocities in L1 system
u_in = 1;
u_up = evalVel_acoustic(x_up_L1,u_in,p,'L1');
u_down = evalVel_acoustic(x_down_L1,u_in,p,'L1');
% Check
u_up2 = s.ccVelIrr(xi_up,u_in);
u_down2 = s.ccVelIrr(xi_down,u_in);

disp(' ')
disp(['Velocities upstream   (',num2str(u_in),' expected): ',num2str(u_up)])
disp(['Velocities downstream (',num2str(u_in*p.Cr),' expected): ',num2str(u_down)])


%% Plot irrotational velocity field
myF = visualize_velField( p , 'fac' , 0.001 );

%% Plot solenoidal velocity field
Gamma = 0.1 * p.u_1_bulkFeed^2;
xs = [ 0.2*p.H_flame + 1i*0.8*p.R_a , 0.6*p.H_flame + 1i*0.8*p.R_a ];
myF2 = visualize_velField( p , 'fac' , 0.01 , 'vort' , {xs,Gamma});



