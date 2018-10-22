function [ vel , ind_sources , L_f ] = distributeSources( myC_L1 , myM , myR0 , SpMM , vel , p , sType )
%DISTRIBUTESOURCES Distributes Sources along a given flame contour line
%
% Inputs:
%   myC_L2 - Coordinates of flame front contour (L1):  [x1 ; x2]  (2xN)
%   myM    - Source strength distribution:        scalar or 1xN vector
%   myR0   - Source radius distribution:          scalar or 1xN vector
%   SpMM   - Sources per milimeter that should be equidistantly distributed
%   vel    - velocity field struct (SC4F_toolbox)
%   p      - Struct containing all required flame info
%   sType  - Distribute sources or vortexes:  'source'   'vort'   ?
%
% Outputs:
%   vel    - velocity field struct that includes ne sources (SC4F_toolbox)
%   ind_sources - Indices at which gas expansion sources have been inserted into vel.sourceDat.x (.G, .xi
%                 and .r0 accordingly)
%   L_f    - length of flame contour
%
% by Thomas Steinbacher Oct 2018
%
% To Dos:
%   - Use line sources
%


%% Compute Number of sources required
% Check input dimeniosn and correct if required
[z,s] = size(myC_L1);
if z>s; myC_L1 = myC_L1.'; end;
% Compute length of curve
L_f = sum( sqrt(sum(diff(myC_L1,[],2).^2,1)) );
% Copute number of sources required
NS = ceil( SpMM * L_f * 1e3 )+1;

% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );


%% Map contour line and distributed data to equidistant grid
doInterpFields = 'noFields';
myFields = cell(2,1);
% radius distribution
if length(myR0)~=1
  doInterpFields = 'fields';
  myFields{1} = myR0;
else
  myFields{1} = myR0 * ones(1,NS);
end

% source strength distribution
if length(myM)~=1
  doInterpFields = 'fields';
  myFields{2} = myM;
else
  myFields{2} = myM * ones(1,NS);
end

% Now interpolate to equidistant 1D grid
[ myC_L1 , myFields_int ] = interp2Dcurve_2equidistantGrid_FPB( myC_L1 , NS , doInterpFields , myFields );

% Length of each "source panel"
myDs = sqrt( diff(myC_L1(1,1:2))^2 + diff(myC_L1(2,1:2))^2 );

% Length flame contour
% L_f = myDs * (NS-1);

% Write interpolated fields to fields cell arrays if interpolation desired
if strcmpi(doInterpFields,'fields')
  myFields = myFields_int;
end

% Source positions should be right between two nodes
myC_L1_v(:,1) = myC_L1(1,1:end-1) + diff(myC_L1(1,:))/2;
myC_L1_v(:,2) = myC_L1(2,1:end-1) + diff(myC_L1(2,:))/2;

% Debug I
% figure;hold on;
% plot(myC_L2(:,1),myC_L2(:,2),'bo');plot(myC_L2_v(:,1),myC_L2_v(:,2),'rx')
% legend({'Flame points','Nodes for sources'})

if strcmpi(sType,'source')
  % Distribute sources
  % Add sourceDat entries if empty
  if isempty(vel.sourceDat)
    vel.sourceDat.x = [];
    vel.sourceDat.xi = [];
    vel.sourceDat.G = [];
  end
  
  % Initial number of sources
  NS0 = length(vel.sourceDat.G);
  if isempty(NS0); NS0=0; end;
  % Write coordinate as complex numbers (L2 system)
  vel.sourceDat.x = [ vel.sourceDat.x ; L1_to_L2(myC_L1_v(:,1) + 1i*myC_L1_v(:,2),p) ];
  % Map points to image domain
  vel.sourceDat.xi = [ vel.sourceDat.xi ; myMap.xi_x(vel.sourceDat.x(end-NS+2:end)) ];
  % Indices where sources have been added
  ind_sources = length(vel.sourceDat.x)-NS+2:length(vel.sourceDat.x);
  
  % Debug II
  % figure;plot(vel.sourceDat.x(2:end))
  % figure;plot(vel.sourceDat.xi(2:end))
  
  
  % Compute source strengths
  % Source radii
  if ~isfield(vel.sourceDat,'r0')
    vel.sourceDat.r0 = zeros(length(vel.sourceDat.G),1);
  end
  myR0_vec = ( myFields{1}(1:end-1) + myFields{1}(2:end) ) / 2;
  vel.sourceDat.r0 = [ vel.sourceDat.r0 ; myR0_vec.' ];
  
  % Source strength distribution
  %   average strength of two neighbouring nodes (linear interpolation)
  myM_vec = ( myFields{2}(1:end-1) + myFields{2}(2:end) ) / 2;
  vel.sourceDat.G = [ vel.sourceDat.G ; myDs*myM_vec.' ];
    
  
elseif strcmpi(sType,'vort')
  % Distribute vortexes     
  % Add vortDat entries if empty
  if isempty(vel.vortDat)
    vel.vortDat.x = [];
    vel.vortDat.xi = [];
    vel.vortDat.G = [];
  end
  
  % Write coordinate as complex numbers (L2 system)
  vel.vortDat.x = [ vel.vortDat.x ; L1_to_L2(myC_L1_v(:,1) + 1i*myC_L1_v(:,2),p) ];
  % Map points to image domain
  vel.vortDat.xi = [ vel.vortDat.xi ; myMap.xi_x(vel.vortDat.x(end-NS+2:end)) ];
  % Indices where sources have been added
  ind_sources = length(vel.vortDat.x)-NS+2:length(vel.vortDat.x);
  
  % Debug II
  % figure;plot(vel.vortDat.x(2:end))
  % figure;plot(vel.vortDat.xi(2:end))
  
  
  % Compute source strengths
  % Source radii
  if ~isfield(vel.vortDat,'r0')
    vel.vortDat.r0 = zeros(length(vel.vortDat.G),1);
  end
  myR0_vec = ( myFields{1}(1:end-1) + myFields{1}(2:end) ) / 2;
  vel.vortDat.r0 = [ vel.vortDat.r0 ; myR0_vec.' ];
  
  % Source strength distribution
  %   average strength of two neighbouring nodes (linear interpolation)
  myM_vec = ( myFields{2}(1:end-1) + myFields{2}(2:end) ) / 2;
  vel.vortDat.G = [ vel.vortDat.G ; myDs*myM_vec.' ];
  
else
  error('Unknown source type!')
end





end

