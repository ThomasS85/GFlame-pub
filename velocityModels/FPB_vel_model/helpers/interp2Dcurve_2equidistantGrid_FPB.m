function [ myCurve_new , fieldData ] = interp2Dcurve_2equidistantGrid_FPB( myCurve , nPoints , varargin )
%INTERP2DCURVE_2EQUIDISTANTGRID Interpolates a given curve to an equistant grid
% 
% Inputs:
%   - myCurve : Points of my curve [ x1_1,x1_2,... ; x2_1,x2_2,... ]; 
%   - nPoints : Number of points
%
% by Thomas Steinbacher

%% Parse varargin
% Should fields be interpolated? -> Specify as cell array where each element contains a field (formate
%  according to myCurve)
ind = find(strcmpi(varargin,'fields'),1);
if ~isempty(ind)
  % Yes
  myFields = varargin{ind+1};

else
  % Default: No
  myFields = [];
end

% Specify interpolated grid?  
%   SepCurve.ds0   : Shift of starting point
%   SepCurve.dsEnd : Shift of end point
ind = find(strcmpi(varargin,'specgrid'),1);
if ~isempty(ind)
  % Yes
  SpeCurve = varargin{ind+1};

else
  % Default: No
  SpeCurve.ds0 = 0;
  SpeCurve.dsEnd = 0;
end

% Make sure flcoord is line vector
[z,s] = size(myCurve);
if z>s; myCurve = myCurve.';doTransp=1;else; doTransp=0; end;


%% Interpolate curve to equidistant grid
% Interpolate to grid specified by user
myDs = sqrt(sum(diff(myCurve,[],2).^2,1));
myDs = [0, myDs]; % add starting point
myS = cumsum(myDs);
myS_new = linspace(SpeCurve.ds0,myS(end)-SpeCurve.dsEnd, nPoints);
myCurve_new = interp1(myS, myCurve.', myS_new);


% Interpolate fields to grid
if ~isempty(myFields)
  fieldData = cell(length(myFields),1);
  for ii=1:length(myFields)
    if length(myFields{ii})==length(myS)
      fieldData{ii} = interp1(myS, myFields{ii}, myS_new);
    else
      % Do not interpolate if length is wrong!
      fieldData{ii} = myFields{ii};
    end
  end
  
else
  fieldData = [];
end

%% If input was row/column vector, also return same format
if ~doTransp
  % Input was column vecor -> tranpose!
  myCurve_new = myCurve_new.';
end

end

