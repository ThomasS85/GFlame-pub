function [ u_1s , u_ref ] = applyKernel_convectiveMods( data , schemeData , tau_vec )
%APPLYKERNEL_CONVECTIVEMODS Function computes axial velocity field based on the history of the forcing signal.
% If used with a Gaussian kernel, the spatial velocity field is convoluted with a Gaussian and the reference
% velocity is in this case the signal at the reference position BEFORE the convolution!
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.de).        //
% // Created, 17.01.2017 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


if strcmpi(schemeData.kernel,'Diffusion')  
  % Diffusion model: Gaussian impulse at domain inlet is diffused along its way according to 1D
  %     advection-diffusion equation
  
  % (1) Define Gaussian kernel
  wc = schemeData.kernelWidth^2;  % Initial variance of the Gaussian (x1=0)
  Df = schemeData.mu;             % Diffusivity/ kinematic viscosity
  myKernel = @(x1,tau) exp( - x1.^2 / (4*Df*tau+2*wc)  ) ./ sqrt( 4*pi*Df*tau + 2*pi*wc );
  
  % (2)  Enlarge length of tau_vec so that no boundary artefacts from convolution arise
  myFac1 = 1.5;
  % 2.1 Store some quantities
  tau_now = tau_vec(1,1);tau_xMax = tau_vec(end,1); d_tau = abs(tau_vec(2,1) - tau_vec(1,1));
  % 2.2 Define length of vectors which shall be appended -> Depends on initial width of Gaussian and its
  %  propagation speed
  tau_fil_front = 3*sqrt(wc) /  schemeData.K;
  tau_fil_back = 3*sqrt(2*Df*abs(tau_xMax-tau_now)+wc) /  schemeData.K;
  % 2.3 Future times (append at front of tau_vec)
  tau_future = tau_now+d_tau:d_tau:tau_now+myFac1*tau_fil_front;
  % 2.4 Past times (append at back of tau_vec)
  tau_past = tau_xMax-d_tau:-d_tau:tau_xMax-myFac1*tau_fil_back;
  % 2.5 Create new tau_vec
  tau_vec = [ tau_future(end:-1:1)'; tau_vec(:,1) ; tau_past' ];
  
  % (3) Get velocity distortion for the domain: substract tau_fil in order to ensure causality
  %  3.1 Get signal and skip tau from future times since these might not be available (Tax coupling!)
  myFac2 = 0.2; % small myFac2: small time delay-> causality bad (information in domain too early)
  % Delay time signal by tau_shift, so that causality still holds 
  tau_shifted = tau_vec - tau_fil_front*myFac2;
  % Get index in shifted tau_vec where future signal ends: Only taus starting from here can be tretrieved from getTransientInput
  ind_now = find(tau_shifted<=tau_now,1,'first');
  % Now get time signal for all past tau
  [ u_s_rel,schemeData] = getTransientInput( schemeData , tau_shifted(ind_now:end) , data );
  % 3.2 Append values which avoid boundary artefacts: linear interpolation (for m = 0 value stays constant!)
  m = (u_s_rel(2) - u_s_rel(1))*0;
  u_s_rel = [ -m*( ind_now-1:-1:1 )'+u_s_rel(1) ;...
    u_s_rel ];
  %   figure;plot(u_s_rel,'o:')
  u_1s_tmp = schemeData.vAmp * u_s_rel * schemeData.meanFlowSpeed;
 
  % (4) Get reference velocity at reference position (1,1) for component 1
  %     Note that this is the velocity fluctuation NOW. At x1=0 we have a signal delayed by tau_fil*myFac2!
  u_ref = u_1s_tmp( ind_now  );
  
  % (5) Diffuse signal according to axial position/ time delay tau (convolution)
  u_1s = zeros(size(tau_vec));
  for ii=1:length(tau_vec)
    % Loop over all axial locations and apply kernel
    % Note, the kernel is growing for mu>0. The minimum kernel size is limited by a minimum tau of 
    %  (2*schemeData.grid.dx(1)^2-2*wc)/4/Df which corresponds to a kernel width close to zero (zero if 
    %  schemeData.grid.dx(1) is zero!). Min. kernel width is schemeData.grid.dx(1)!
    u_1s = u_1s + myKernel( ((1:length(tau_vec))'-ii)*schemeData.grid.dx(1) ,...
                            max(tau_vec(length(tau_future)+1,1)-tau_vec(ii,1) , (2*schemeData.grid.dx(1)^2-2*wc)/4/Df ) ) ...
                  *schemeData.grid.dx(1)*u_1s_tmp(ii);
                
    % Debug
%     if mod(ii,ceil(length(tau_vec)/5))==0 
%       x_tmp = ((ii-40:ii+40)'-ii)*schemeData.grid.dx(1);
%       figure(f1);plot(x_tmp,myKernel(x_tmp,tau_vec(1,1)-tau_vec(ii,1)));xlabel('x [mm]')
%       figure(f2);plot(u_1s)
%     elseif ii==1
%       f1 = figure;hold on;
%       f2 = figure;hold on;plot(u_1s)
%       line([length(tau_future) length(tau_future)],[0 schemeData.vAmp*schemeData.meanFlowSpeed]*1.2,'LineStyle',':')
%     end
    
  end
  
  % (6) Cut u_1s so that its length is same same as the x1 vector's length
  u_1s = u_1s( length(tau_future)+1 : length(tau_future)+schemeData.grid.N(1) );
  
  % (7) Assign signal to velocity field
  u_1s = repmat( u_1s , 1 , schemeData.grid.N(2) );
  
  % Debug
%   figure('Position',[0 50 1000 800]); hold on;
%   x1_max = schemeData.grid.vs{1}(end); x1_min = schemeData.grid.vs{1}(1);dx1=schemeData.grid.dx(1);
%   yMax = schemeData.vAmp*schemeData.meanFlowSpeed;
%   x1 = schemeData.grid.vs{1,1};
%   x1_ext = [linspace(x1_min-dx1*length(tau_future),x1_min-dx1,length(tau_future))';x1;...
%     linspace(x1_max+dx1,x1_max+dx1*length(tau_past),length(tau_past))'];
%   p1=plot(x1,u_1s(:,1),'b');
%   p2=plot(x1_ext,u_1s_tmp,'r-.');
%   x_tmp = ((ii-40:ii+40)'-ii)*schemeData.grid.dx(1);
%   kernelMin=myKernel(x_tmp,0);kernelMax=myKernel(x_tmp,abs(tau_xMax-tau_now));
%   p3 = plot(x_tmp,kernelMin/max(kernelMin)*yMax,'g');
%   p4 = plot(x_tmp+schemeData.grid.vs{1}(end),kernelMax/max(kernelMax)*yMax,'g:');
%   line([x1_max x1_max],[-yMax yMax]*1.4,'LineStyle',':');
%   line([0 0],[-yMax yMax]*1.4,'LineStyle',':')
%   ylim([-yMax yMax]*1.4);xlabel('x_1 [mm]');ylabel('u'' [m/s]')
%   legend([p1;p2;p3;p4],{'convoluted';'original';'min. kernel width';'max. kerne width'})
%   disp('Debug Plot')
  
  
else
  % (C) Tradiational model: Dirac impulse
  
  % (C.1) Get distributed distortion signal
  [ u_s_rel,schemeData] = getTransientInput( schemeData , tau_vec , data );
  u_1s = schemeData.vAmp * u_s_rel * schemeData.meanFlowSpeed;
  
  % Get reference velocity at reference position (1,1) for component 1
  u_ref = u_1s(1);
  
end


end

