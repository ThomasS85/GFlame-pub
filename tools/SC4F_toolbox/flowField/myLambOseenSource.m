function dOmega_dxi = myLambOseenSource(xi,xi_s,Gamma,r_0,varargin)
% Returns velocity due to one single Lamb-Oseen-like source + its mirror source with kerne radius r_0 in image domain.

%% Parse varargin
% Should mirror sources be added, in order to maintain zero flux BC at real axis?
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
  if imag(xi_s)==0
    % If Source is at real axis, mirroring is equal to doubling the source strength!
    dOmega_dxi =  Gamma/pi * ( 1 - exp( -abs(xi-xi_s).^2. / (r_0/2).^2 ) )  ./ ( xi - xi_s );
    dOmega_dxi(xi==xi_s) = 0;
  else
    % Add mirror source
    dOmega_dxi =  Gamma/(2*pi) * (( 1 - exp( -abs(xi-xi_s).^2. / (r_0/2).^2 ) )  ./ ( xi - xi_s )...
      +( 1 - exp( -abs(xi-conj(xi_s)).^2. / (r_0/2).^2 ) )  ./ ( xi - conj(xi_s )));
    % Limit for xi->xi_s
    dOmega_dxi(xi==xi_s) = Gamma/(2*pi) *...
      ( 1 - exp(-abs(xi(xi==xi_s)-conj(xi_s)).^2./(r_0/2).^2) ) ./ ( xi(xi==xi_s) - conj(xi_s) );
  end
else
  % Add no mirror source
  dOmega_dxi =  Gamma/(2*pi) * ( 1 - exp( -abs(xi-xi_s).^2. / (r_0/2).^2 ) )  ./ ( xi - xi_s );
  dOmega_dxi(xi==xi_s) = 0;
end


end
