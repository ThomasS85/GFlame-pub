function [ u_res , dOmega_dxi ] = evalFlowField_PointSource( myXi , vel , varargin )
%EVALFLOWFIELD Evaluates the flow field, resulting from a number of vortexes at x_v, a number of sources
% x_s and constant part u_p in the xi plane.
%
% Inputs:
%   myXi          - Array of complex numbers which specify spatial positions where velocity should be evaluated
%   vel.vortDat   - Struct which contains information about all vortexes
%                   .xi  : Nx1 complex array with positions of all vortexes [ 1+1i*2 ; 2+1i*3 ]
%                   .G   : Nx1 real array with circulations of each vortex [ 2 ; 4 ]
%                   .r0  : Nx1 real array with kernel width of each vortex [ 3 ; 4]
%                           If set to zero, a point vortex is assumed, otherwise a Lamb vortex
%  vel.sourceDat  - Struct which contains information about all (point) sources
%                   .xi   : Nx1 complex array with positions of all sources [ 1+1i*2 ; 2+1i*3 ]
%                   .G   : Nx1 real array with strengths of each source [ 2 ; 4 ]
%  vel.u_p        - Complex constant velocity which is added to each point
%
% Optional inputs:
%  myMapping      - Struct with mapping information. Required to compute kernel radius for Lamb-Orseen Vortex
%                   and to apply Routh's rule for point vortex
%
% Outputs
%  u_res      - Array of complex numbers with complex velocities at all positions myX
%  dOmega_dxi - Array of complex numbers with derivative of complex velocity potential at all positions myX
%
% by Thomas Steinbacher (02.2018)

%% Parse varargin
% Should mapping be applied? Only relevant for r0!
ind = find(strcmpi(varargin,'myMapping'),1);
if ~isempty(ind)
  % Use user specified mapping
  myMap = varargin{ind+1};
else
  % Default: No mapping
  myMap.dxi_dx = @(xi) 1;
end

% Unwrap velocity field data
vortDat = vel.vortDat;
sourceDat = vel.sourceDat;
u_p = vel.u_p;
doMirror = vel.doMirror;


%% Evaluate vortexes
dOmega_dxi_v = zeros(size(myXi));
if ~isempty(vortDat)
  % If r0 does not exists, assume point vortexes!
  if ~isfield(vortDat,'r0')
    vortDat.r0 = zeros(size(vortDat.xi));
  end
  
  % Now loop over vortexes
  for ii=1:length(vortDat.xi)
    % Compute contribution of i-th vortex
    if vortDat.r0(ii) == 0
      % Point Vortex
      if ~strcmpi(doMirror,'noMirror') && ~(imag(vortDat.xi(ii))==0)
        % Add mirror
        dOmega_dxi_new = -1i * vortDat.G(ii)/(2*pi) * ( vortDat.xi(ii) - conj(vortDat.xi(ii)) ) ./ ...
          ( ( myXi - vortDat.xi(ii) ) .* ( myXi - conj(vortDat.xi(ii)) ) );
        % Special treatment at vortex center (xi=xi_s): Evaluate mirror vortex contribution at vortex center
        dOmega_dxi_new(myXi==vortDat.xi(ii)) = 1i * vortDat.G(ii)/(2*pi) ./ ...
          ( myXi(myXi==vortDat.xi(ii)) - conj(vortDat.xi(ii)) );
      else
        % Add no mirror
        dOmega_dxi_new = -1i * vortDat.G(ii)/(2*pi) ./ ( myXi - vortDat.xi(ii) );
        % Special treatment at vortex center (xi=xi_s): Set velocity to zero
        dOmega_dxi_new(myXi==vortDat.xi(ii)) = 0;
      end
    else
      % Lamb-Orseen Vortex
      r0 = abs( myMap.dxi_dx(vortDat.xi(ii))  ) * vortDat.r0(ii);
      dOmega_dxi_new = myLambOseenVortex( myXi , vortDat.xi(ii) , vortDat.G(ii) , r0 , doMirror );
    end
    
    % Add contribution of i-th vortex
    dOmega_dxi_v = dOmega_dxi_v + dOmega_dxi_new;
    
  end
end

clear dOmega_dxi_new;


%% Evaluate sources
dOmega_dxi_s = zeros(size(myXi));
if ~isempty(sourceDat)
  % If r0 does not exists, assume point sources!
  if ~isfield(sourceDat,'r0')
    sourceDat.r0 = zeros(size(sourceDat.xi));
  end
  
  % Now loop over sources
  for ii=1:length(sourceDat.xi)
    if  sourceDat.r0(ii) == 0
      % Compute contribution of i-th source
      if ~strcmpi(doMirror,'noMirror') %&& ~isinf(sourceDat.x(ii))
        % Add mirror
        if imag(sourceDat.xi(ii))==0
          % If Source is at real axis, mirroring is equal to doubling the source strength!
          dOmega_dxi_new = sourceDat.G(ii)/pi ./ ( myXi - sourceDat.xi(ii) );
          % Set velocity to zero at vortex center
          dOmega_dxi_new( myXi==sourceDat.xi(ii) ) = 0;
        else
          dOmega_dxi_new = sourceDat.G(ii)/(2*pi) * ( 2*myXi - sourceDat.xi(ii) - conj(sourceDat.xi(ii)) ) ./ ...
            ( ( myXi - sourceDat.xi(ii) ) .* ( myXi - conj(sourceDat.xi(ii)) ) );
          % Only evaluate mirror source contribution at source center
          dOmega_dxi_new(myXi==sourceDat.xi(ii)) = sourceDat.G(ii)/(2*pi) ./ ...
            ( myXi(myXi==sourceDat.xi(ii)) - conj(sourceDat.xi(ii)) );
        end
      else
        % Add no mirror
        dOmega_dxi_new = sourceDat.G(ii)/(2*pi) ./ ( myXi - sourceDat.xi(ii) );
        % Set velocity to zero at vortex center
        dOmega_dxi_new( myXi==sourceDat.xi(ii) ) = 0;
      end
    else
      % Distribute source strength analog to Lamp-Orseen Vortex
      r0 = abs( myMap.dxi_dx(sourceDat.xi(ii))  ) * sourceDat.r0(ii);
      dOmega_dxi_new = myLambOseenSource( myXi , sourceDat.xi(ii) , sourceDat.G(ii) , r0 , doMirror );
    end
    
    % Add contribution of i-th source
    dOmega_dxi_s = dOmega_dxi_s + dOmega_dxi_new;
    
  end
end


%% Add constant contribution and return resulting velocity
% Velocity is complex conjugate of dOmega_dxi
dOmega_dxi = dOmega_dxi_v + dOmega_dxi_s;
u_res = conj(dOmega_dxi) + u_p;


end

