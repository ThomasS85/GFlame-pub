clear
close all

addPath2Bib('GFlame-pub')


%% Set parameters and nondimensionalise
% Set default values for parameters
[ p ] = setUpPredefinedFlame( 'Steinbacher_2D' );
tau_r = p.L_flame / ( p.u_1_bulkFeed * cos(p.alpha) );



%% Check mapping
% Point in L2 system
x1 = [ 0.5*p.H_flame + 1i*(p.R_a-p.R_i) , 1.5*p.H_flame + 1i*(p.R_a*0.5) , -0.1*p.H_flame + 1i*(p.R_a-0.5*p.R_i) ];
% Map to image domain
for ii=1:length(x1)
  [ xi_SC(ii) ] = SCmapInv_SCFT( x1(ii) , p );
end
% Map back to physical domain
[ s ] = return_SCmap_SCFT( p );
x1_2 = s.x_xi(xi_SC);

disp(' ')
disp('Checking mapping:')
disp(['x original: ',num2str(x1)])
disp(['x mapped  : ',num2str(x1_2)])

% Plot points
[ myF ] = sketchCombustor_SCFT( p );
xlim([-0.3*p.H_flame 2*p.H_flame])
for ii=1:length(x1)
  myCoord = x1(ii);
  plot( myCoord ,'xr')
  text(real(myCoord),imag(myCoord),num2str(ii))
end

[ myF2 ] = sketchCombustor_SCFT( p , 'image' );
ax=gca; xlim([ax.XLim(1) 30]); ylim( [ax.YLim(1) , ax.YLim(2)*2] )
for ii=1:length(x1)
  myCoord = xi_SC(ii);
  plot( myCoord ,'xr')
  text(real(myCoord),imag(myCoord),num2str(ii))
end


%% Check velocity
% Points far upstream of flame (L2)
x_up = [ -0.1 + 1i*(p.R_a-0.9*p.R_i) , -0.1 + 1i*p.R_a ];
% Points far downstream of flame (L2)
x_down = [ p.H_flame*1.3 + 1i*0.1*p.R_i , p.H_flame*1.3 + 1i*p.R_a ];

% Map
[ xi_up ] = SCmapInv_SCFT( x_up , p );
[ xi_down ] = SCmapInv_SCFT( x_down , p );

disp(' ')
disp('Checking mapping for velocity evluation:')
disp(num2str(s.x_xi(xi_up)));
disp(num2str(s.x_xi(xi_down)))
% Compute velocities (complex conjugate of result since we are in system L2)
u_in = 1;
u_up = conj(s.ccVelIrr(xi_up,u_in));
u_down = conj(s.ccVelIrr(xi_down,u_in));

disp(' ')
disp(['Velocities set (inlet) : ',num2str(u_in)])
disp(['Velocities upstream    : ',num2str(u_up)])
disp(['Velocities downstream  : ',num2str(u_down)])


%% Check for branch cuts
% Define points in image domain
xi1 = linspace( -4 , 6 , 2e2);
xi2 = linspace( -1 , 2 , 2e2);
[Xi1 , Xi2] = meshgrid(xi1,xi2);
Xi_all = Xi1(:) + 1i*Xi2(:);

% Map to X domain (L2)
[ s ] = return_SCmap_SCFT( p );
X_all = s.x_xi(Xi_all);
X_all = reshape(X_all,size(Xi1));

% Check dxi_dx
dXi_dX_all = s.dxi_dx(Xi_all);
dXi_dX_all = reshape(dXi_dX_all,size(Xi1));

% Check d2xi_dx2
d2Xi_dX2_all = s.d2xi_dx2(Xi_all);
d2Xi_dX2_all = reshape(d2Xi_dX2_all,size(Xi1));

% Check dx_dxi
dX_dXi_all = s.dx_dxi(Xi_all);
dX_dXi_all = reshape(dX_dXi_all,size(Xi1));


% plot
% X_all
figure('Name','real X_all');
surf(Xi1,Xi2,real(X_all));

figure('Name','imag X_all');
surf(Xi1,Xi2,imag(X_all))

% dXi_dX_all
figure('Name','real dXi_dX_all');
surf(Xi1,Xi2,real(dXi_dX_all))

figure('Name','imag dXi_dX_all');
surf(Xi1,Xi2,imag(dXi_dX_all))

% d2Xi_dX2_all
figure('Name','real d2Xi_dX2_all');
surf(Xi1,Xi2,real(d2Xi_dX2_all))

figure('Name','imag d2Xi_dX2_all');
surf(Xi1,Xi2,imag(d2Xi_dX2_all))

% dX_dXi_all
figure('Name','real dX_dXi_all');
surf(Xi1,Xi2,real(dX_dXi_all))

figure('Name','imag dX_dXi_all');
surf(Xi1,Xi2,imag(dX_dXi_all))



%% Examples
% Square root
sqrtXi = sqrt(Xi_all);
sqrtXi = reshape(sqrtXi,size(Xi1));

% Logarithm
logXi = log(Xi_all);
logXi = reshape(logXi,size(Xi1));

% cosh^-1
acoshXi = acosh(Xi_all);
acoshXi = reshape(acoshXi,size(Xi1));

% sqrtXi
figure('Name','real sqrtXi');
surf(Xi1,Xi2,real(sqrtXi))

figure('Name','imag sqrtXi');
surf(Xi1,Xi2,imag(sqrtXi))

% sqrtXi
figure('Name','real logXi');
surf(Xi1,Xi2,real(logXi))

figure('Name','imag logXi');
surf(Xi1,Xi2,imag(logXi))

% acoshXi
figure('Name','real acoshXi');
surf(Xi1,Xi2,real(acoshXi))

figure('Name','imag acoshXi');
surf(Xi1,Xi2,imag(acoshXi))
