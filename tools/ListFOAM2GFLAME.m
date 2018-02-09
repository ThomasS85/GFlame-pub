function ListFOAM2GFLAME( fileSource , fileDestination, varargin )
%LISTFOAM2FFLAME Convets a file containing an openFOAM list to a list
%readable by GFLAME
%   
% Using default settings file is just converted to different format
% The following options can be chosen:
%   1) Resampe data to deltax spacing : ('resample',dx)
%   2) Shift signal to specified start time: ('shiftT0',t0)

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Do resample time series?
ind = find(strcmp(varargin,'resample'),1);
if ~isempty(ind)
  dx = varargin{ind+1};
else
  % Don't resample
  dx = 0;
end

% Shift start of transient signal?
ind = find(strcmp(varargin,'shiftT0'),1);
if ~isempty(ind)
  % Set desired t0
  t0 = varargin{ind+1};
else
  % Don't shift
  t0 = -1;
end

% Postpone signal due to to convective delay
ind = find(strcmp(varargin,'postponeSig'),1);
if ~isempty(ind)
  dt_postpone = varargin{ind+1};
else
  % Don't shift
  dt_postpone = -1;
end

% Minimum time which should be included in time series. if < than t0, fill
% with zeros
ind = find(strcmp(varargin,'solverT0'),1);
if ~isempty(ind)
  solverT0 = varargin{ind+1};
else
  % Don't shift
  solverT0 = -1;
end


% Open and read FOAM list
if exist(fileSource,'file')==2
  f1 = fopen(fileSource,'r');
else
  error(['Please make sure a file <',fileSource,'> with time series exists in the case root directory.'])
end
data_mat = cell2mat( textscan(f1,'( %f %f )','Headerlines',2) );
fclose(f1);

% Resample if desired
if dx>0
  x_vec_new = ( data_mat(1,1):dx:data_mat(end,1) )';
  data_mat_new(:,2) = interp1(data_mat(:,1),data_mat(:,2),x_vec_new);
  data_mat_new(:,1) = x_vec_new;
  data_mat = data_mat_new;
end

% Shift if desired so that transient time series starts at t0. If chosen t0
% of solver is smaller than shifted t0, fill everything between solver t0
% and transient t0 with zeros (later)
if t0>=0
  dt_shift1 = data_mat(1,1) - t0;
  data_mat(:,1) = data_mat(:,1) - dt_shift1;
end

% Shift and add zeros for shifted time span if desired; doesn't change t0
if dt_postpone > 0
  % Original start time
  t_start = data_mat(1,1);
  % Time discretisation
  dt = data_mat(2,1) - data_mat(1,1);
  % Shift
  data_mat(:,1) = data_mat(:,1) + dt_postpone;
  % add new times
  newTimes = (t_start : dt : data_mat(1,1)-dt )' ;
  data_mat = [ zeros( length(newTimes) , 2 ) ; data_mat ];
  data_mat(1:length(newTimes),1) = newTimes;
end

% Add zeros from t0 solver to t0 transient start
if solverT0<t0
  % Find time step just before transient start
  t_constLast = data_mat(1,1)-(data_mat(2,1)-data_mat(1,1));
  if t_constLast > solverT0
    % Add time step before transient start with zero excitation
    data_mat = [ solverT0 , 0 ; t_constLast , 0 ; data_mat ];
  else
    % Transient start was too close to solverT0
    data_mat = [ solverT0 , 0 ; data_mat ];
  end
end
  
% Write 2 file
[pathstr] = fileparts(fileDestination);
if exist(pathstr,'dir')==7
  f1 = fopen(fileDestination,'w');
  fprintf(f1,'xData , yData');
  fprintf(f1,'\n%.8f , %.8f',[data_mat(:,1)';data_mat(:,2)']);
  % fprintf(f1,'\n( %.9f %.9f )',[t_vec';u_vec']);
  fclose(f1);
else
  error('Please initialize simulation with v.type=constant in the velocity settings!')
end

end

