function [ myMean , myN , x1F , dx1F , ind_sources ] = getMeanFlameFlow( p , Nx1F , do_I , method_I , sourceInfo)
%GETMEANFLAMEFLOW Summary of this function goes here
%   Detailed explanation goes here

%% Prepare gas expansion modeling
% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );

% Potential flow
vel_mean.sourceDat.G = p.u_1_bulkFeed*p.R_i;
vel_mean.sourceDat.x = -inf;
vel_mean.sourceDat.xi = 0;
% Initialize vel
vel_mean.vortDat = [];
vel_mean.u_p = 0;
vel_mean.doMirror = 'doMirror';

if strcmpi(p.CombType,'backwardFacingStep')
  % Enable Kutta condition
  vel_mean.Kutta.H = 1;
  vel_mean.Kutta.active = 1;
  doKutta = 'doKutta';
else
  % Disable Kutta condition
  doKutta = 'noKutta';
end

% use default: straight line with angle alpha
if ~strcmp(p.flameType,'flat')
  x1L_tmp = linspace(0,p.H_flame,Nx1F).';
  x2L_tmp = p.R_flame - x1L_tmp*tan(p.alpha);
  % flame tangent vector (L1)
  myT = [  cos(p.alpha) ; -sin(p.alpha) ];
else
  x1L_tmp = zeros(Nx1F,1);
  x2L_tmp = linspace(p.R_i,0,Nx1F).';
  % flame tangent vector (L1)
  myT = [  0 ; 1 ];
end
myMean.flcoord = [ x1L_tmp , x2L_tmp ];

% Compute mean flow component in flame parallel direction
myMean.u_parallel = [ ones(Nx1F,1)*p.u_1_centerIn , zeros(Nx1F,1) ]*myT;

% Interpolate specified mean flame front coordinates + parallel velocity field to specified flame coordinates
[ myMean.flcoord , myMean.u_parallel ] = interp2Dcurve_2equidistantGrid_FPB( myMean.flcoord.' , Nx1F , 'fields' , {myMean.u_parallel.'} );
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
dx1F = x1F(2) - x1F(1);

% Compute normal/ tangential vector to mean flame front
[ myN ] = comNorm2FlameFront( myMean.flcoord , x1F );

% [nDat] = compNormal2Curve_FPB( myMean.flcoord );
% nDat.x1F_nodes =  x1F(1:end-1) + dx1F/2;
% % Get normals at x1F nodes via interpolation -> exclude boundaries!
% myN = interp1( nDat.x1F_nodes , nDat.n_vec(1,:)+1i*nDat.n_vec(2,:) , x1F(2:end-1) );
% % Add values at boundaries manually
% myN = [ myN(1) , myN , myN(end) ];
% % Normalize (probably not required, but to be sure...)
% myN = myN ./ sqrt( real(myN).^2 + imag(myN).^2 );

ind_sources = [];
if do_I
  % Evaluate mean flame front dilatation effects
  dispXi = sourceInfo(1);
  myR0 = sourceInfo(2);
  SpMM = sourceInfo(3);
  myC_L1 = myMean.flcoord(1,:)+ 1i*myMean.flcoord(2,:) + dispXi .* myN;   % Mean source contour
  [ myC_L1 ] = checkNcorrectLine_geometry( myC_L1 , p );
  myM = ( p.E - 1 ) * p.s_l_u;     % Source strength
  if strcmpi(method_I,'source')
    [ myMean.vel , ind_sources ] = distributeSources( [real(myC_L1);imag(myC_L1)] , myM , myR0 , SpMM , vel_mean , p , 'source' );
  else
    [ myMean.vel ] = distributePanels( [real(myC_L1);imag(myC_L1)] , myM , SpMM , vel_mean , p , 'source' );
  end
  
else
  myMean.vel = vel_mean;
end

end

