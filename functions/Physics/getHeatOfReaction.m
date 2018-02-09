function [ delta_h_R ] = getHeatOfReaction( mixture , Phi )
%GETHEATOFREACTION Computes heat of reaction 
%
% Inputs:
%   mixture     : Name of the mixture
%   Phi     : Equivalence ratio (scalar or matrix)
%
% Outputs:
%   delta_h_R : Heat of reaction [J/kg]
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.04.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

if strcmpi(mixture,'MethaneAir')
  % For methan at 300K and atmospheric pressure
  %   Lieuwen, T., 2003. Modeling Premixed Combustion - Acoustic Wave Interactions:
  %     A Review. J. of Propul. and Power 19, 765–781.
  
  % Parameters
  a = 2.9125e6;
  b = 0.05825;
  
  delta_h_R = a * min(1,Phi) ./ ( 1 + b*Phi );
  
else
  error('Heat of reaction calculation for chosen mixture not implemented!')
end


end

