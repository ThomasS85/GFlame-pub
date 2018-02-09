function showLog( GFcase, caseNum )
%CASEINFO Shows (latest) logfile of a specific run 
%   
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 07.04.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% Change dir to case root
change2caseRootDir( GFcase )
% get case folder
caseFolder = [GFcase.p.caseFolder,'/',num2str( caseNum )]; 
% Load case 
load([caseFolder,'/init.mat'],'GFcase');
% get name of logfile
logFile = GFcase.p.logFileName;
type(logFile)

end

