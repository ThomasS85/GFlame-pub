function copyCase( source , destination )
%COPYCASE Copies a case i to a case j updating all the information in
%GFcase
%
% Use this function e.g. if you want to copy a steady state solution and
% use this for sevaral excitations or settings
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 12.06.2015 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Check if source and target are integer
if mod(source,1)==0 && mod(destination,1)==0
  sourceFold = num2str(source);
  destinationFold = num2str(destination);
else
  error('Either source or target are not integer values!')
end

%% Check if source exits
if ~ (exist(['.',filesep,sourceFold],'file') == 7)
  error('Source folder does not exist!')
end

%% Check if target already exists
if exist(['.',filesep,destinationFold],'file')==7
  error('Target folder already exists!')
end

%% Copy files
copyfile(sourceFold,destinationFold)

%% update data in destination folder
% find *.mat files
d = dir(destinationFold);
isub = ~[d(:).isdir]; %# returns logical vector indicating all not directories
nameFolds = {d(isub).name}'; % write out all files to cell aray of cells
nameFolds = regexp(nameFolds,'^\d+\.\d+','match'); % delete all entries but those with numeric data
nameFolds(cellfun('isempty',nameFolds)) = []; % delete empty entries
nameFolds = cellfun(@cell2mat,nameFolds,'UniformOutput',false); % convert to array of strings

for ii = 1:length(nameFolds)
  % Load GFcase from *.mat file
  load([destinationFold,filesep,nameFolds{ii},'.mat']);
  
  % update information
  GFcase.p.caseNumber = destination;
  GFcase.p.run_output_folder = [GFcase.p.caseFolder,filesep,destinationFold];
  % get source logfile name and update it
  logFileName  = strsplit(GFcase.p.logFileName,filesep);
  logFileName = logFileName{end};
  GFcase.p.logFileName = [GFcase.p.caseFolder,filesep,destinationFold,filesep,logFileName];
  
  % Save GFcase again
  save([destinationFold,filesep,nameFolds{ii},'.mat'],'GFcase','data');
  
end


end

