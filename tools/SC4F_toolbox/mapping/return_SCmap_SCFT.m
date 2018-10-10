function [ s ] = return_SCmap_SCFT( p )
%RETURN_SCMAP Function returns everything which is required to perform specific SC mapping and computation of
% velocities
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, Jan 2018                                  //
% ////////////////////////////////////////////////////////


%% Get appropriate mapping
if strcmpi(p.CombType,'backwardFacingStep')
  % Backward facing step combustor
  % Mapping by Brown/ Churchill fig. 30
  % Define mapping of point B
  s.h = 1 / p.Cr^2;
  % reference length to calculate real lengths
  s.l_ref = p.R_a / pi;
  
  % SC-Mapping (normalized coordinates)
  z1 = @(xi) (2*xi-s.h-1) / (s.h-1);
  z2 = @(xi) ((s.h+1)*xi-2*s.h) ./ ((s.h-1)*xi);
  s.x_xi = @(xi) s.l_ref * ( acosh( z1(xi) ) - 1/sqrt(s.h) * acosh( z2(xi) ) ) ;
  
  % Derivative of x with respect to xi
  s.dx_dxi = @(xi) s.l_ref * sqrt( ( xi - 1 ) ./ ( xi - s.h ) ) ./ xi;
  
  % Derivative of xi with respect to x
  % root_z1 = @(xi) sqrt( z1(xi)+1 ) .* sqrt( z1(xi)-1 );
  % root_z2 = @(xi) sqrt( z2(xi)+1 ) .* sqrt( z2(xi)-1 );
  % s.dxi_dx_old =  @(xi) (s.h-1)/2 ./ ( 1./root_z1(xi) - sqrt(s.h)./(xi.^2.*root_z2(xi))  ) / s.l_ref ;
  s.dxi_dx =  @(xi) 1/s.l_ref * sqrt( ( xi - s.h ) ./ ( xi - 1 ) ) .* xi;
  
  % Second derivative of xi with respect to x
  s.d2xi_dx2 = @(xi) 1/s.l_ref * sqrt( ( xi - 1 ) ./ ( xi - s.h ) ) ...
    .* ( 2*xi.^2 - xi * (3+s.h) + 2*s.h ) ./ ( 2*( xi - 1 ).^2 );
  
  % Routh's correction (xi and Gamma can be vector of same length)
  s.RouthsCorr = @(Gamma,xi)  -1i*Gamma/(4*pi) .* s.d2xi_dx2( xi ) ./ s.dxi_dx( xi );
  
  % Inverse mapping
  s.xi_x = @(x) SCmapInv_SCFT( x , p );             % Input x in L2 system
  s.xi_x_L1 = @(x) SCmapInv_SCFT( x , p , 'L1' );   % Input x in L1 system
  
  % Flow field due to source at -infinity in physical domain (Irrotational flow)
  s.dOmega_dxi_S = @(xi,u_in) s.l_ref * u_in / sqrt(s.h) ./ ( xi );
  s.ccVelIrr = @(xi,u_in) s.dOmega_dxi_S(xi,u_in) .* s.dxi_dx(xi); % complex conjugate in L2 and velocity in L1
  
  % Vertex and prevertex positions
  s.vertexes = s.l_ref * [ 1i*pi*( 1 - p.Cr ) , 0 , inf+1i*pi , -inf+1i*pi ];
  s.prevertexes = [ 1 , 1/p.Cr^2 , inf , 0 ];
  
elseif strcmpi(p.CombType,'duct')
  % Duct flame
  % Mapping by Brown/ Churchill fig. 6
  % reference length to calculate real lengths
  s.l_ref = p.R_i / pi;                        
  % SC-Mapping (normalized coordinates) xi->x
  s.x_xi = @(xi) s.l_ref * log(xi);
  % SC-Mapping (normalized coordinates) x->xi
  s.xi_x = @(x) exp( x / s.l_ref );
  % Derivative of x with respect to xi
  s.dx_dxi = @(xi) 1 / ( xi / s.l_ref );
  % Derivative of xi with respect to x
  s.dxi_dx =@(xi) xi / s.l_ref;
  % Second derivative of xi with respect to x
  s.d2xi_dx2 = @(xi) xi / s.l_ref^2; %-xi^2;                           
  
  % Routh's correction (xi and Gamma can be vector of same length)
  s.RouthsCorr = @(Gamma,xi)  -1i*Gamma/(4*pi) .* s.d2xi_dx2( xi ) ./ s.dxi_dx( xi );
  
else
  error('Unknown combustor geometry!')
end

end

