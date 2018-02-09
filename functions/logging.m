function logging( GFcase , code , varargin )
%LOGGING creates log file for GFLAME and exits programm in case of error
%
% Codes for logging:
%     0 : Initialize log file
%     1 : Append message to log file
%     2 : Append warning to log file
%     3 : Append error to log file and terminate praogramme with error (not
%         used at the moment!)
%     4 : Programme finished successfull
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Get filename
filename = GFcase.p.logFileName;


% look if message is provided in varargin
if nargin > 2
  msg = varargin{1};
else
  msg = '<not provided>';
end

% opens text file
if code == 0
  % Initialise logging facility
  fileID = fopen(filename,'w');
  
  fprintf(fileID,['Programmstart: ', datestr(clock), '\n']);
  fprintf(fileID,['\tInitiated by ',upper(msg),'\n']);
  fprintf(fileID,'---------------------------------------');
  
  % Print genereal settings p
  fprintf(fileID,'\n-General settings:');
  fieldNames = fieldnames(GFcase.p);
  for ii = 1:length(fieldNames)
    if sum( strcmp( fieldNames{ii} , GFcase.p.logfile ) )
    fprintf(fileID,'\n\t%25s :\t%s',fieldNames{ii},num2str( GFcase.p.(fieldNames{ii}) ) );
    end
  end
  % Print solver settings solver
  fprintf(fileID,'\n-Solver settings:');
  fieldNames = fieldnames(GFcase.solver);
  for ii = 1:length(fieldNames)
    if sum( strcmp( fieldNames{ii} , GFcase.solver.logfile ) )
      fprintf(fileID,'\n\t%25s :\t%s',fieldNames{ii},num2str( GFcase.solver.(fieldNames{ii}) ) );
    end
  end
  % Print settings for reinitialisation
  fprintf(fileID,'\n-Reinitialization settings:');
  fieldNames = fieldnames(GFcase.solver.reinit);
  for ii = 1:length(fieldNames)
    if sum( strcmp( fieldNames{ii} , GFcase.solver.reinit.logfile ) )
      fprintf(fileID,'\n\t%25s :\t%s',fieldNames{ii},num2str( GFcase.solver.reinit.(fieldNames{ii}) ) );
    end
  end
  % Print vel settings vel
  fprintf(fileID,'\n-Velocity model settings:');
  fieldNames = fieldnames(GFcase.vel);
  for ii = 1:length(fieldNames)
    if sum( strcmp( fieldNames{ii} , GFcase.vel.logfile ) )
      fprintf(fileID,'\n\t%25s :\t%s',fieldNames{ii},num2str( GFcase.vel.(fieldNames{ii}) ) );
    end
  end
  % print grid information
  fprintf(fileID,'\n-Grid information:');
  fprintf(fileID,'\n\t%25s :\t%sx%s','grid size',num2str( GFcase.grid.N(1) ) , num2str( GFcase.grid.N(2) ) );
  fprintf(fileID,'\n\t%25s :\t%s','Boundary Conditions x1',func2str(GFcase.grid.bdry{1,1}) );
  fprintf(fileID,'\n\t%25s :\t%s','Boundary Conditions x2',func2str(GFcase.grid.bdry{2,1}) );
  fprintf(fileID,'\n---------------------------------------');
  
else
  % Append file
  fileID = fopen(filename,'a');
end


% log
switch code
  case 0
    fprintf(fileID,'\n');
  case 1
    % Write Log
    fprintf(fileID,['\n\t', msg]);
  case 2
    % Warning
    fprintf(fileID,['\nWarning: ', msg]);
  case 3
    % Error
    fprintf(fileID,['\nERROR: ', msg]);
    fprintf(fileID,'\n\nProgram terminated!\n');
    fprintf(fileID,'---------------------------------------');
    fprintf(fileID,['\n', datestr(clock)]);
    fclose(fileID);
    error(msg)
  case 4
    % Programme finishd successfull
    fprintf(fileID,'\n\nProgram finishd successfully!\n');
    fprintf(fileID,'---------------------------------------');
    fprintf(fileID,['\n', datestr(clock)]);
  otherwise
    error('unknown code for logging!');
end

fclose(fileID);

end

