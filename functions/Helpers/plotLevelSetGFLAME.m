function [ h ] = plotLevelSetGFLAME( GFcase ,varargin )
%PLOTLEVELSETGFLAME Plots data in a GFcase struct
%   
%
%% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 02.03.2015 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% unwrap data
g = GFcase.grid;
data = GFcase.data;

%% Parse varargin
% Plot G-field or only zero contour?
ind = find(strcmpi(varargin,'plotG'),1);
if ~isempty(ind)
  % Plot G-field
  plotG = 1;
else
  % Default: No
  plotG = 0;
end

% Show zero plane?
ind = find(strcmpi(varargin,'show0'),1);
if ~isempty(ind)
  % Show zero plane
  show0 = 1;
else
  % Default: No
  show0 = 0;
end

% Also plot grid?
ind = find(strcmpi(varargin,'showGrid'),1);
if ~isempty(ind)
  % Show grid
  showGrid = 1;
else
  % Default: No
  showGrid = 0;
end

% Also visualize grid?
ind = find(strcmpi(varargin,'showGrid'),1);
if ~isempty(ind)
  % Plot G-field
  showGrid = 1;
else
  % Default: No
  showGrid = 0;
end



%% Start plotting
hold on
if plotG
  % Plot G-field
  visualizeLevelSet(g, data, 'surf', 0, 'GFLAME');
  view(-37.5,30);
  if show0
    % Plot zero plane
    surf( g.xs{1}, g.xs{2}, zeros(size(g.xs{1})) );
  end 
end

% Plot zero contour (default)
h = visualizeLevelSet(g, data, 'contour', 0, 'GFLAME');


% Plot grid if desired
if showGrid
  visualizeGrid(g,[],'g');
end

end

