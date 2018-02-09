function [ p ] = setUserParameters( p,userParameters )
%SETUSERPARAMETERS Replaces all desired fields of p with values provided
% by the user

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% Write field names of p to array
fieldNames = fieldnames(p);

% Define black list (settings which cannot be modified
blacklist = {'dim'};

% Remove all non string elements in userParameters (which are parameters)
n = 1;
for jj = 1:length(userParameters)
  if ~ischar(userParameters{jj})
    nonStringElement{n} = userParameters{jj};
    userParameters{jj} = ['nonStringVal_',num2str(n)];
    n = n+1;
  end
end

% iterate over all fields of p
for ii = 1:length(fieldNames)
  % Iterate ofer all user supplied fields
  for jj = 1:length(userParameters)
    % Check if one supplied field matches a field in p and is not in
    % blacklist!
    if  strcmpi(fieldNames{ii},userParameters{jj}) && ~sum(strcmpi(blacklist,userParameters{jj}))
      % If so overwrite the matching field with user supplied data
      if strfind(userParameters{jj+1},'nonStringVal_')==1
        tmpVal = strsplit(userParameters{jj+1},'_');
        tmpVal = str2num(tmpVal{end});
        p.(fieldNames{ii}) = nonStringElement{tmpVal};
      else
        p.(fieldNames{ii}) = userParameters{jj+1};
      end
      break;
    end
    
  end
  
end


end

