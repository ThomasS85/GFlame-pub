function [ GFcase ] = prepare4output( GFcase , varargin )
%PREPARE4OUTPUT prepares everything for output of GFLAME. 
%   Function creates an output folder if there is none and initiates log
%   file. This function is always called when a GFLAME routine with output
%   is started (like integrateLevelSet or loadFOAMCase). This way during
%   set up (setUpCase) nothing has to be written to hard disc).
%
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 15.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% unwrap data
p = GFcase.p;

% look if message is provided in varargin
if nargin > 1
  msg = varargin{1};
else
  msg = '<not provided>';
end

%% Create folder for output if not yet existent
if ~(exist(p.run_output_folder,'dir') == 7)
  mkdir(p.run_output_folder);
end

%% Logging
if ~p.logFileInit
  % Get filename for logfile
  p.logFileName = [p.run_output_folder,'/log.out'];
  ii = 1;
  while exist(p.logFileName,'file')
    p.logFileName = [p.run_output_folder,'/log',num2str(ii),'.out'];
    ii = ii + 1;
  end
  % Wrap data
  GFcase.p = p;
  % initialize logfile
  logging( GFcase , 0 , msg )
  % Set indicator if logfile has already been initialized to 1 (true)
  p.logFileInit = 1;
  
end

%% wrap data
GFcase.p = p;

end

