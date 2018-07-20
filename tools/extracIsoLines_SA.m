function [ cFuns ] = extracIsoLines_SA( grid , data , level )
%EXTRACISOLINES_SA Extracts level-iso line of a surface data and returns all lines separately stored in a
% cell-array cFuns 
%
%   -> code taken from lineInegral_contour() (replace there by call of extracIsoLines_SA??!)
%
% Inputs:
%   grid     - grid struct
%   data     - surface
%   level    - level of iso-line to extract
%
% Outputs:
%   cFuns    -  cell array with all iso-lines
%
% by Thomas Steinbacher July 2018

%% Extract iso line
C = contourc(grid.vs{1}, grid.vs{2}, data', [level,level]);


%% Check if multiple iso lines were returned and split them
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
% write out number of contour lines
numC = ii;


%% Convert contour lines and save all to cell array cFuns
cFuns = cell(numC,1);
for ii = 1:numC
  % Convert jjth contour line to unique functions in x1 (assumption: Contour lines are continuous!)
  
  % (1) Get contour line to temporarily vector (x1 and x2 components separated)
  tempC_x1 = C( 1 , indI(ii)+1:indI(ii)+nE(ii) );
  tempC_x2 = C( 2 , indI(ii)+1:indI(ii)+nE(ii) );
  
  % debug
%   figure;plot(tempC_x1,tempC_x2,'b');hold on;
%   for jj = 1:length(tempC_x1)
%     plot(tempC_x1(jj),tempC_x2(jj),'rx')
%     waitforbuttonpress 
%   end
  
  % (2) Write all contour lines to cell array of unique functions
  cFuns{ii} = [ tempC_x1 ; tempC_x2 ];
  
end

end

