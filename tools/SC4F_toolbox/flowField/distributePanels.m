function [ vel , L_f ] = distributePanels( myC_L2 , myM , PpMM , vel , p , sType , varargin )
%DISTRIBUTESOURCES Distributes Panels along a given flame contour line
%
% Panels are defined within a cell array vel.panel:
% -> Assume we have NS panels:
%   { sType_1    ,  sType_2     ,  ...  ;    string that defines panel type: 'source' or 'vort' 
%     strength_1 , strength_2   ,  ...  ;     Nx1    vector including all panel strengthss
%     nodesX_1   ,  nodesX_2    ,  ...  ;    (N+1)x1 vector in x-coord. including all nodes as complex numbers
%     nodesXi_1  ,  nodesXi_2   ,  ...  ;    (N+1)x1 vector in xi-coord. including all nodes as complex numbers
%     nDat_1     ,  nDat_2      ,  ...   }   output of compNormal2Curve() that includes all normal/ tangential
%                                             vectors
% 
% 
% Inputs:
%   myC_L2 - Coordinates of flame front contour (L2):  [x1 ; x2]  (2xN)
%   myM    - Source strength distribution:        scalar or 1xN vector
%   PpMM   - Panels per milimeter that should be equidistantly distributed
%   vel    - velocity field struct (SC4F_toolbox)
%   p      - Struct containing all required flame info
%   sType  - Distribute sources or vortexes:  'source'   'vort'   ?
%
% Outputs:
%   vel    - velocity field struct that includes ne panels (SC4F_toolbox)
%   ind_sources - Indices at which gas expansion sources have been inserted into vel.sourceDat.x or 
%                  vel.vortDat (.G, .xi and .r0 accordingly)
%   L_f    - length of flame contour
%
% by Thomas Steinbacher Oct 2018
%
% To Dos:
%   - Use line sources
%

%% Parse varargin
% At which position ii of vel.panel{:,ii+1} should new panel be written? LEAVEs NO EMPTY COLUMNS!
ind = find(strcmpi(varargin,'write2'),1);
if ~isempty(ind)
  % User specified mean flame position
  write2 = varargin{ind+1};
else
  % No mean flame position
  write2 = [];
end

% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );


%% Compute Number of sources required
% Check input dimeniosn and correct if required
[z,s] = size(myC_L2);
if z>s; myC_L2 = myC_L2.'; end;
% Compute length of curve
L_f = sum( sqrt(sum(diff(myC_L2,[],2).^2,1)) );
% Copute number of panels required
NS = ceil( PpMM * L_f * 1e3 );


%% Map contour line and distributed data to equidistant grid
doInterpFields = 'noFields';
% source strength distribution
if length(myM)~=1
  doInterpFields = 'fields';
  myFields{1} = myM;
else
  myFields{1} = myM * ones(1,NS+1);
end

% Now interpolate to equidistant 1D grid
[ myC_L2_intp , myFields_int ] = interp2Dcurve_2equidistantGrid_FPB( myC_L2 , NS+1 , doInterpFields , myFields );

% Length of each "source panel"
% myDs = sqrt( diff(myC_L2_intp(1:2,1))^2 + diff(myC_L2_intp(1:2,2))^2 );

% Write interpolated fields to fields cell arrays if interpolation desired
if strcmpi(doInterpFields,'fields')
  myFields = myFields_int;
end

% Debug I
% figure;hold on;
% plot(myC_L2(1,:),myC_L2(2,:),'bx:');
% plot(myC_L2_intp(:,1),myC_L2_intp(:,2),'ro--');
% legend({'Original','Interpolated'})

% Add panel if none exists so far
if ~isfield(vel,'panel')
  vel.panel = {'type';'strengths';'x';'xi';'normals/tangential'};
end
% Get number of panels that have already been defined
numbP = length(vel.panel(1,:))+1;

% Write panel to desired position
if ~isempty(write2)
  numbP = min( numbP , write2+1 );
end

% Initialize new panel -> Write panel type
vel.panel{1,numbP} = sType;

% Write panel nodes in physical domain
vel.panel{3,numbP} = L1_to_L2( (myC_L2_intp(1,:) + 1i*myC_L2_intp(2,:)).' , p );

% Write panel nodes in image domain
vel.panel{4,numbP} = myMap.xi_x(vel.panel{3,numbP});

% Compute and write panel normal/tangential information
[nDat] = compNormal2Curve( myC_L2_intp );
vel.panel{5,numbP} = nDat;

% Write panel strengths
%   average strength of two neighbouring nodes (linear interpolation)
myM_vec = ( myFields{1}(1:end-1) + myFields{1}(2:end) ) / 2;
myFacMap = abs( myMap.dx_dxi( vel.panel{4,numbP} ) );
myFacMap = ( myFacMap(1:end-1) + myFacMap(2:end) ) / 2;
vel.panel{2,numbP} = myFacMap .* myM_vec.';


end

