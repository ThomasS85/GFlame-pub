function file = openAndCheckFile( solver, fileName  )
%CHECKOUTFILE Opens file for appending and checks if file for output has additional entries and deletes
%the if true
%   Can be either file for flame surface or reference velocity

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 08.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% open file for appending since data will be appended by GFLAME)
file = fopen(fileName,'r');

% read header
header = fgetl(file);

% Read data from file
data = cell2mat( textscan(file,'%f , %f','Headerlines',0) );

% Find position of start time 
dt_max = max(diff(data(:,1)));
ind_t0 = find (data(:,1) - solver.t0 <= dt_max ,1,'last');
% ind_t0 = find (data(:,1) <= solver.t0 ,1,'last');

% close file
fclose(file);

% Check in case ind_t0 is last entry if time gap to solver.t0 is not too
% big (maximum 10 times the last saved time step!)
if ind_t0 == size(data,1) && (solver.t0-data(ind_t0,1))>10*(data(ind_t0,1)-data(ind_t0-1,1))
  error(['Starting time of simulation is far bigger than last saved time step data in <',fileName,'>!'])
end

% if this is not at the end, deleta entries after this and update file
if ind_t0 < size(data,1)
  % Delete old file
  file = fopen(fileName,'w');
  % delte entries
  data(ind_t0+1:end,:) = [];
  % write data 2 file
  fprintf(file,header);
  fprintf(file,'\n%.8f , %.8f',data');
  
else
  % Open the old file for appending data
  file = fopen(fileName,'a');
  
end
  
end

