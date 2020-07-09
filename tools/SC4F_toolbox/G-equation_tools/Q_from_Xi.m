function [ Qs_Qm ] = Q_from_Xi( Xi ,x1F , p , varargin )
%Q_FROM_XI Evaluates the (linearized) surface integral based on flame front perturbation Xi (flame normal)
%  given over x1F (flame coordinate)
%
% Inputs:
%   - Xi      : Distribution of xi over x1F
%   - x1F     : Spatially discretized flame coordinates with constant spacing (dx1F)
%
%   Outputs:
%           - Qs_Qm         : Normalized heat release rate fluctuation
%
%
% TO DO: INCLUDE CORRECT CURVATURE DEPENDENT FLAME SPEED TREATMENT !!!!!!!!!!!!!!!!!!!!!!!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 10.11.2017                                //
% // Last modified: 10.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////

%% Parse varargin
% Should difusion be considered? (change in flame speed)
ind = find(strcmpi(varargin,'doDiffusion'),1);
if ~isempty(ind)
  % Yes
  curvContri = varargin{ind+1};
else
  % Default: No
  curvContri = 0;
end


%% Evaluuate heat release
if strcmp(p.Fdim,'2D')
  % Slit flames
  Qs_Qm = Xi(end) / (p.L_flame*tan(p.alpha));
  
elseif strcmp(p.Fdim,'3D')
  % Conical flames
  switch p.flameType
    case 'inverseV'
      Qs_Qm = 2 * cos(p.alpha) / (p.R_flame*p.L_flame) * trapz(x1F,Xi);
    case 'V'
      Qs_Qm = 2 * cos(p.alpha) / (p.R_flame*p.L_flame) * ( p.L_flame*Xi(end) - trapz(x1F,Xi) );
    otherwise
      error('Unknown flame type!')
  end
end

end

