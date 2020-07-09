function dOmega_dxi = myLambOseenVortex(xi,xi_s,Gamma,r_0,varargin)
% Returns velocity due to one single Lamb-Oseen Vortex + its mirror vortex with kerne radius r_0 in image domain.

%% Parse varargin
% Should mirror vortexes and sources be added, in order to maintain zero flux BC at real axis?
ind = find(strcmpi(varargin,'noMirror'),1);
if ~isempty(ind)
  % No
  addMirror = 0;
else
  % Default: Yes
  addMirror = 1;
end

%% Return derivative of potential
if addMirror
  % Add mirror vortex
  dOmega_dxi = 1i * Gamma/(2*pi) * ( ( xi_s - conj(xi_s) ) ./ ...
    ( ( xi - xi_s ) .* ( xi - conj(xi_s) ) ) ...
    + exp( -abs(xi-conj(xi_s)).^2 ./ r_0.^2 ) ./ ( xi - conj(xi_s) ) ...
    - exp( -abs(xi-xi_s).^2. / r_0.^2 ) ./ ( xi - xi_s )    );
  % Limit for xi->xi_s
  dOmega_dxi(xi==xi_s) = - 1i * Gamma/(2*pi) *...
    ( 1 - exp(-abs(xi(xi==xi_s)-conj(xi_s)).^2./r_0.^2) ) ./ ( xi(xi==xi_s) - conj(xi_s) );
else
  % Add no mirror vortex
  dOmega_dxi = 1i * Gamma/(2*pi) * ( 1 - exp( -abs(xi-xi_s).^2. / r_0.^2 ) )  ./ ( xi - xi_s );
end


end

