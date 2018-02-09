function change2caseRootDir( GFcase )
%CHANGE2CASEROOTDIR Changes working directory to the root directory of the
%case
%   Many helper functions require the case root directory as working
%   directory

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


cd(GFcase.p.caseRootDir)

end

