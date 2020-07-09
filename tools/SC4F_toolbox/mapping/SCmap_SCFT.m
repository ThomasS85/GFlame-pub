function [ x ] = SCmap_SCFT( xi , p , varargin )
%SC_MAPPING Maps a given complex vector with locations in the physical
%       domain to the image domain using inverse Schwarz-Christoffel mapping
%
%   Mapping xi -> x
%
%   Input:  - xi          : Complex vector with coordinates in image
%                           domain xi = xi1 + i*xi2
%           - p           :  Struct with flame settings as returned from setUpPredefinedFlame()
%
%   Output: - x           : Complex vector with coordinates in physical
%                           domain (L2 or L1 if desired)
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 28.05.2015 as part of GFLAME 0.1          //
% // Last modified: 28.05.2015 by steinbacher           //
% ////////////////////////////////////////////////////////

%% Check varargin
% Should result of mapping be mapped to L1 system?
ind = find(strcmpi(varargin,'L1'),1);
if ~isempty(ind)
  % Do map from L2 to L1
  mapL2L1 = 1;
else
  % Do not map to L2
  mapL2L1 = 0;
end


%% Perform SC Mapping
% Get mapping infos
[ s ] = return_SCmap_SCFT( p );
% Perform SC-Mapping (normalized coordinates)
x = s.x_xi(xi);
% Map to L1 if desired
if mapL2L1
  x = L2_to_L1(x,p);
end


end

