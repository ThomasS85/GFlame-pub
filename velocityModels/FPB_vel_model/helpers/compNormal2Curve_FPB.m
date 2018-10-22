function [nDat] = compNormal2Curve_FPB( myCurve )
% Function computes normal and tangential vectors to a curve myCurve
% 
% Inputs:
%   - myCurve : Points of my curve [ x1_1,x1_2,... ; x2_1,x2_2,... ]; 
%
% by Thomas Steinbacher

% Negative x1 values sign since elements are ordered starting from flame tip to base
nDat.d_x1 = diff(myCurve(1,:));
nDat.d_x2 = diff(myCurve(2,:));

% Compute normal and tangential vectors at positions between to nodes
nDat.n_vec = [ -nDat.d_x2 ; nDat.d_x1 ];
nDat.t_vec = [ nDat.d_x1 ; nDat.d_x2 ];
% Normalize to length 1
normFac = sqrt( nDat.d_x1.^2 + nDat.d_x2.^2 ); % Negative sign so that normal vectors point downstream direction
nDat.n_vec(1,:) = nDat.n_vec(1,:) ./ normFac; nDat.n_vec(2,:) = nDat.n_vec(2,:) ./ normFac;
nDat.t_vec(1,:) = nDat.t_vec(1,:) ./ normFac; nDat.t_vec(2,:) = nDat.t_vec(2,:) ./ normFac;

% Compute nodes at which normal vectors where evaluated
nDat.nodes = [ myCurve(1,1:end-1) + nDat.d_x1/2  ; myCurve(2,1:end-1) + nDat.d_x2/2 ];
% Number of nodes
nDat.Nnodes = length(nDat.nodes(1,:));

end