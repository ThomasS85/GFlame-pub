function [ band ] = identifyInterfaceBand( data , bandWidth )
%EXTRACTBAND Identifies a band around the zero interface of width 
%
% This function identifies a band of width 2*bandWidth around the zero
% level iso line an returns a logic matrix with information where this band
% is located
%
% Sparse should be used here since a lot of indexing has to be performed
% (performance tests have been carried out)
%
%   Only for dimension 2!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 30.03.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

[rowData,colData] = size(data);

%% Find nodes at interface

% Look in x1 direction
cutData = data( 1:end-1 , : );
nodesTrueX1 = (cutData .* data(2:end,:) < 0);
% Copy the last line
nodesTrueX1(end+1,:) = nodesTrueX1(end,:);

% Look in x2 direction
cutData = data( : , 1:end-1 );
nodesTrueX2 = (cutData .* data(:,2:end) < 0) ;
% Copy the last row
nodesTrueX2(:,end+1) = nodesTrueX1(:,end);

% Nodes with zero are on the interface
nodesTrueZero = ( data == 0 );

% Write band (faster with full since there is much indexing taking
% place in the following. This has been tested!)
band = ( nodesTrueX1 | nodesTrueX2 | nodesTrueZero );

 

%% Add nodes in a band around interface
% Find indices of non zero entries
[rows,cols]=find(band);

% Add true nodes to upper end x1 direction
% remove those which are closer than bandWidth to the end
nodesTrue = rows<rowData-bandWidth & rows>bandWidth;
rowsU = rows(nodesTrue);
colsU = cols(nodesTrue);

for ii=1:length(rowsU)
  band( rowsU(ii)-bandWidth:rowsU(ii)+bandWidth , colsU(ii) ) = 1;
end

% Treat all nodes close to boundaries
boundaryRows = rows( ~nodesTrue );
boundaryCols = cols( ~nodesTrue );

for ii=1:length(boundaryRows)
  % close to upper or lower boundary?
  if boundaryRows(ii) <= bandWidth
    % Closer to lower boundary
    band( 1:2*bandWidth , boundaryCols(ii) ) = 1;
  else
    % Closer to upper boundary
    band( boundaryRows(ii)-bandWidth:end , boundaryCols(ii) ) = 1;
  end
end



% Add true nodes to upper end x2 direction
% remove those which are closer than bandWidth to the end
nodesTrue = cols<colData-bandWidth & cols>bandWidth;
rowsU = rows(nodesTrue);
colsU = cols(nodesTrue);

for ii=1:length(colsU)
  band( rowsU(ii) , colsU(ii)-bandWidth:colsU(ii)+bandWidth ) = 1;
end

% Treat all nodes close to boundaries
boundaryRows = rows( ~nodesTrue );
boundaryCols = cols( ~nodesTrue );

for ii=1:length(boundaryCols)
  % close to upper or lower boundary?
  if boundaryCols(ii) <= bandWidth
    % Closer to lower boundary
    band( boundaryRows(ii) , 1:2*bandWidth ) = 1;
  else
    % Closer to upper boundary
    band( boundaryRows(ii) , boundaryCols(ii)-bandWidth:end ) = 1;
  end
end



end

