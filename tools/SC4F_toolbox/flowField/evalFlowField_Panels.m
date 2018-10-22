function [ u_res , dOmega_dxi ] = evalFlowField_Panels( myXi , vel  )
%EVALFLOWFIELD Evaluates the flow field, resulting from a number of source or vortex panels in the xi plane.
%
% Inputs:
%   myXi          - Array of complex numbers which specify spatial positions where velocity should be evaluated
%   vel.panel   - Struct which contains information about all vortexes
%         -> Assume we have NS panels:
%           { sType_1    ,  sType_2     ,  ...  ;    string that defines panel type: 'source' or 'vort'
%             strength_1 , strength_2   ,  ...  ;     Nx1    vector including all panel strengthss
%             nodesX_1   ,  nodesX_2    ,  ...  ;    (N+1)x1 vector in x-coord. including all nodes as complex numbers
%             nodesXi_1  ,  nodesXi_2   ,  ...  ;    (N+1)x1 vector in xi-coord. including all nodes as complex numbers
%             nDat_1     ,  nDat_2      ,  ...   }   output of compNormal2Curve() that includes all normal/ tangential
%                                                     vectors
%
% Outputs
%  u_res      - Array of complex numbers with complex velocities at all positions myXi
%  dOmega_dxi - Array of complex numbers with derivative of complex velocity potential at all positions myXi
%
% by Thomas Steinbacher (10.2018)

%% Check if any panel is contained
if ~ isfield(vel,'panel')
  u_res = 0;
  dOmega_dxi = 0;
  return
end

%% Add mirror panels?
if ~strcmpi(vel.doMirror,'noMirror')
  % Yes
  fac1 = 1;
else
  % No
  fac1 = 0;
end

%% Evaluate panels
% Initialize derivative potential
dOmega_dxi = zeros(size(myXi));
% Loop over all panels
for pp=2:length(vel.panel(1,:))
  % Based on panel-type evaluate veloctiy field
  if strcmpi(vel.panel{1,pp},'vort')
    % (A) VORTEX PANEL
    dOmega_dxi_vortP_up = @(xi,xi1,xi2) 1i / (2*pi) * ...
      log( ( xi - xi2 ) ./ ( xi - xi1 ) ) * exp(-1i*angle(xi2 -xi1));
    dOmega_dxi_vortP_down = @(xi,xi1,xi2) -1i / (2*pi) * ...
      log( ( xi - conj(xi2) ) ./ ( xi - conj(xi1) ) ) * exp( 1i*angle(xi2 -xi1)) * fac1;
    dOmega_dxi_fun = @(G_p,xi,xi1,xi2) G_p*(dOmega_dxi_vortP_up(xi,xi1,xi2) + dOmega_dxi_vortP_down(xi,xi1,xi2));
  elseif strcmpi(vel.panel{1,pp},'source')
    % (B) SOURCE PANEL
    dOmega_dxi_sourceP_up = @(xi,xi1,xi2) -1 / (2*pi) * ...
      log( ( xi - xi2 ) ./ ( xi - xi1 ) ) * exp(-1i*angle(xi2 -xi1));
    dOmega_dxi_sourceP_down = @(xi,xi1,xi2) -1 / (2*pi) * ...
      log( ( xi - conj(xi2) ) ./ ( xi - conj(xi1) ) ) * exp( 1i*angle(xi2 -xi1)) * fac1;
    dOmega_dxi_fun = @(Q_p,xi,xi1,xi2) Q_p* (dOmega_dxi_sourceP_up(xi,xi1,xi2) + dOmega_dxi_sourceP_down(xi,xi1,xi2));
  else
    warning('Unknown panel type: skipped!')
  end
  
  % Loop over panels and evaluate their respective effect
  for ii=1:length(vel.panel{2,pp})
    dOmega_dxi = dOmega_dxi + dOmega_dxi_fun( vel.panel{2,pp}(ii) , myXi , vel.panel{4,pp}(ii) , vel.panel{4,pp}(ii+1) );
  end
  
end


%% Return resulting velocity
% Velocity is complex conjugate of dOmega_dxi
u_res = conj(dOmega_dxi);


end

