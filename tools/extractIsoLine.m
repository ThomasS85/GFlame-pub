function [ C ] = extractIsoLine( GFcase , caseNum , t_load , level,varargin)
%EXTRACTISOLINE Function extracts iso line from given 2D G-field
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2015 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Which fields should be decomposed?
ind = find(strcmpi(varargin,'changeOriginX1'),1);
if ~isempty(ind)
  % Use user specified fields
  deltaOriginX1 = varargin{ind+1};
else
  % use default decompositions: Helmholtz Hodge
  deltaOriginX1 = 0;
end

% Load only 1 time step or all?
if t_load<=0
  % Load all time teps
  % define where to load data from
  caseNumber = GFcase.p.caseNumber;
  run_output_folder = ['.',filesep,num2str( caseNumber )];
  % Find available timesteps
  d = dir(run_output_folder);
  isub = ~[d(:).isdir]; %# returns logical vector indicating all not directories
  t_load = {d(isub).name}'; % write out all files to cell aray of cells
  t_load = regexp(t_load,'^\d+\.\d+','match'); % delete all entries but those with numeric data
  t_load(cellfun('isempty',t_load)) = []; % delete empty entries
  t_load = cellfun(@cell2mat,t_load,'UniformOutput',false); % convert to array of strings
  t_load = cellfun(@str2num,t_load); % Convert to vector of numerics
  
  % Initialize cell arry for output
  C_all = cell(length(t_load),1);
end


%% load data
for tt = 1:length(t_load)
  % get path to matlab file which is to be loaded
  fileLoad = [GFcase.p.caseRootDir,filesep,num2str(caseNum),filesep,num2str(t_load(tt),12),'.mat'];
  if exist(fileLoad,'file') == 2
    load(fileLoad,'data')
  else
    error('Unable to load chosen time step data!')
  end
  
  %% Extract iso line
  [C,~] = contour(GFcase.grid.xs{1}, GFcase.grid.xs{2}, data, [level,level]);
  hold on;
  
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
  
  
  %% Change origin if desired
  C(1,2:end) = C(1,2:end) + deltaOriginX1;
  
  % Write to cell array for output
  if length(t_load)>1
    C_all{tt,1} = C; 
  end
  
end
%% Save iso line(s) to disk
if length(t_load)>1
  fileSave = [GFcase.p.caseRootDir,filesep,num2str(caseNum),filesep,'all','_contour.mat'];
   save(fileSave,'C_all')
else
  fileSave = [GFcase.p.caseRootDir,filesep,num2str(caseNum),filesep,num2str(t_load(tt)),'_contour.mat'];
   save(fileSave,'C')
end
  
  
end

