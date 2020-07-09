function [ grid ] = map_grid(X,Y,p,varargin )
% this function maps the entire grid from the physical in the image domain
% required:
% X  x-coordinates
% Y  y-coordinates
% p contains information about the flame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  by Axel Zimmermann 05.2018 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse varargin
ind = find(strcmpi(varargin,'L1'),1);
if ~isempty(ind)
  % Do map from L1 to L2
  mapL1L2 = 1;
else
  % Do not map to L2; x already in L2!
  mapL1L2 = 0;
end

%% Map grid
myGridFile = [p.run_output_folder,filesep,'mappedGrid_',num2str(length(X(:,1))*length(X(1,:))),'.mat'];
if exist(myGridFile,'file') == 2
  % If grid could be loaded from file -> do it!
  load(myGridFile)
  
else
  % Writing grid as a complex-coordinates
  grid.x=X+1i*Y;
  
  if mapL1L2 && strcmpi(p.CombType,'backwardFacingStep')
    grid.x = L1_to_L2(grid.x,p);
  end
 
  % preparing output
  [size_y,size_x] = size(grid.x);
  grid.xi = zeros(size_y,size_x);
  grid_neu = zeros(size_y,size_x);
  
  % Now map
  disp('Mapping grid... (this may take some time)')
  for ii=1:1:size_x
    grid.xi(:,ii) = SCmapInv_SCFT( grid.x(:,ii) , p , 'L1' );
    grid_neu(:,ii) = SCmap_SCFT( grid.xi(:,ii) , p , 'L1' );
  end
  
  if find((grid.x-grid_neu)>1e-7,1)
    error('>Mapping was not successful')
  else
    disp('  >mapping grid was successful! Grid written to disc!');
    save(myGridFile,'grid')
  end
  
end

end

