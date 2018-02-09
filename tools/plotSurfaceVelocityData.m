function plotSurfaceVelocityData( GFcase , varargin )
%PLOTSURFACEDATA Plots temporal evolution of the flame surface of an case
%   

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

%% Parse varargin
% Should case be change?
ind = find(strcmp(varargin,'case'),1);
if ~isempty(ind)
  % Change run output folder
  GFcase.solver.resume2case = varargin{ind+1};
else
  % Load case specified by case number
  GFcase.solver.resume2case = GFcase.p.caseNumber;
end

% Should root folder be changed? be change?
ind = find(strcmp(varargin,'root'),1);
if ~isempty(ind)
  % Change root folder
  GFcase.p.caseRootDir = varargin{ind+1};
else
  % Leave everything as it is
end

% Smoothing parameter
ind = find(strcmp(varargin,'smoothing'),1);
if ~isempty(ind)
  scaleSmoothing = varargin{ind+1};
  doSmooth = 1;
else
  % default value
  scaleSmoothing = 20;
  doSmooth = 0;
end



% Change dir to case root
change2caseRootDir( GFcase )

% define where to load data from
GFcase.p.run_output_folder = [GFcase.p.caseFolder,'/',num2str( GFcase.solver.resume2case )]; 

% Read surface data from file
f1 = fopen([GFcase.p.run_output_folder,'/Surface.out'],'r');
data_A = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
fclose(f1);

% Read reference velocity data from file
f1 = fopen([GFcase.p.run_output_folder,'/V_ref.out'],'r');
data_V = cell2mat( textscan(f1,'%f , %f','Headerlines',1) );
fclose(f1);


% Smooth data using the Savitzky Golay Filter
if doSmooth
  degreePoly = 3;
  derivative = 0;
  frameSize = 2*scaleSmoothing + 1; % Must be odd!
  A_SavGol = savitzkyGolayFilt(data_A(:,2)',degreePoly,derivative,frameSize,[],2);
end

% plot data
figure
subplot(2,1,1)
plot( data_A(:,1) , data_A(:,2) ); hold on
if doSmooth
  plot( data_A(:,1) , A_SavGol );
  legend({ 'Original data' , ['Savitzky Golay Filter (Smoothing: ',num2str(scaleSmoothing),')'] },'Location','southeast')
else
  legend('Original data','Location','southeast')
end
xlabel('t [s]')
ylabel('A_{fl}')
title(['Flame Surface (Case ',num2str(GFcase.solver.resume2case),')'])


subplot(2,1,2)
plot(data_V(:,1),data_V(:,2))
xlabel('t [s]')
ylabel('v_{ref}')
title('Reference Velocity')

end

