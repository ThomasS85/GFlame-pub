function [ C ] = extract_IsoLine( grid, data , level)
%EXTRACTISOLINE Function extracts iso line from given 2D G-field
%
% Input:
% grid    -struct with grid information
%         .x grid coordinates in physical domain
% data    - G-Field 
%
% 
%  Output:
%  C      - isoLine coordinates in physical domain
%
%
%
%
% ////////////////////////////////////////////////////////
% // Main program created by:                           //
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    // 
% // Created, 01.12.2015 as part of GFLAME 0.1          //
% // Last modified: 05.06.2018 by Axel Zimmermann       //
% ////////////////////////////////////////////////////////


 %% Extract iso line
C = contourc(grid.vs{1}, grid.vs{2}, data', [level,level]);
  %% Check if multiple iso lines were returned and return only the longest one
  doIterate = 1; ii = 1;nE_tmp = 0;
  while doIterate
    % Get number of elements and start index of iso line
    nE(ii) = C(2,ii+nE_tmp);
    indI(ii) = ii+nE_tmp;
    if sum(nE)< length(C(2,:))-ii
      % Another iso line is stored in C
      nE_tmp = nE_tmp + nE(ii);
      ii = ii + 1;
    else
      % No other iso line is stored in C
      doIterate = 0;
    end
  end
  % Extract iso line with maximum index
  indMax = indI(nE==max(nE));
  C = C(:,indMax:indMax+max(nE));
  
  % Check order how x1 values are sorted -> Sort ascendingly!
  if C(1,2)>C(1,end)
    C = [ C(:,1) , C(:,end:-1:2) ];
  end
end

