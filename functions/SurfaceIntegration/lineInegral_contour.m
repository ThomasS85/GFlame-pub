function [ integralLine , data ] = lineInegral_contour( grid, data , myFun , ~ )
%LINEINEGRAL_DIRECT Function computes the surface integral of a function myFun along the zero isoline of the
% data field.
%
% Only 2D!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Extract iso line
C = contourc(grid.vs{1}, grid.vs{2}, data', [0,0]);


%% Check if multiple iso lines were returned and split them
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
% write out number of contour lines
numC = ii;


%% Convert contour lines and save all to cell array cFuns
cFuns = cell(numC,1);
for ii = 1:numC
  % Convert jjth contour line to unique functions in x1 (assumption: Contour lines are continuous!)
  
  % (1) Get contour line to temporarily vector (x1 and x2 components separated)
  tempC_x1 = C( 1 , indI(ii)+1:indI(ii)+nE(ii) );
  tempC_x2 = C( 2 , indI(ii)+1:indI(ii)+nE(ii) );
  
  % debug
%   figure;plot(tempC_x1,tempC_x2,'b');hold on;
%   for jj = 1:length(tempC_x1)
%     plot(tempC_x1(jj),tempC_x2(jj),'rx')
%     waitforbuttonpress 
%   end
  
  % (2) Write all contour lines to cell array of unique functions
  cFuns{ii} = [ tempC_x1 ; tempC_x2 ];
  
end


%% Define function which has to be integrated
if isequal(size(myFun), size(data))
  % myFun is an array of the same size as data
  
  % Only perform calculation in a small band around zero line
  bandWindow = 2;
  band = ( identifyInterfaceBand( data , bandWindow ) );
  x1_band = grid.xs{1}(band);
  x2_band = grid.xs{2}(band);
  myFun_band = myFun(band);

  %  -> interpolate myFun to center points on contour lines if myFun is not a scalar
  F = scatteredInterpolant( x1_band , x2_band , myFun_band );
  
elseif isscalar(myFun)
  % myFun is a scalar
  F = @(x1,x2) ones(size(x1)) * myFun;
  
else
  % myFun is NEITHER an array of the same size as data NOR a scalar -> error
  error('Specified function which should be integrated must either be a scalar or have the same size as data!')
  
end


%% Iterate over all contour lines and integrate myFun along line
integralLine = 0;
for ii = 1:numC 
    % Get dx1, dx2 and compute length of line segments
    dx1Vals = diff( cFuns{ii}(1,:) );
    dx2Vals = diff( cFuns{ii}(2,:) );
    ds = sqrt( dx1Vals.^2 + dx2Vals.^2 );
    
    % Evaluate myFun at all points of contour line
    funC = F( cFuns{ii}(1,:) , cFuns{ii}(2,:) );
    
    % Use trapezoidal method to integrate myFun along contour line
    integralLine = integralLine + sum( ( funC(1:end-1) + funC(2:end) ) / 2 .* ds );
    
    % Debug
%     if 1
%       f = figure('Position',[0 50 1200 1000]);
%       subplot(3,1,1); plot(cFuns{ii}(1,:),cFuns{ii}(2,:));hold on;
%       plot(cFuns{ii}(1,1),cFuns{ii}(2,1),'ro');
%       plot(cFuns{ii}(1,10),cFuns{ii}(2,10),'rx');
%       xlabel('x_1')
%       subplot(3,1,2); plot(cFuns{ii}(1,:),funC/2/pi);hold on;
%       plot(cFuns{ii}(1,1),funC(1)/2/pi,'ro');
%       plot(cFuns{ii}(1,10),funC(10)/2/pi,'rx');
%       plot( cFuns{ii}(1,2:end)-0.5*dx1Vals , ( funC(1:end-1) + funC(2:end) ) / 2 /2/pi, 'g:');
%       title('funC/2/\pi');xlabel('x_1');
%       subplot(3,1,3); plot( ds);  title('ds');
%       close(f)
%     end
    
end
  

end

