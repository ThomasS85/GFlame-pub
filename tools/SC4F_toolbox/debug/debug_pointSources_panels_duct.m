clear
close all

addPath2Bib('GFlame-pub')

%% Settings
% Include gas expansion?
do_I = 0;

% panel or source method (gas expansion)?   panel   or source
method_I = 'source';
% Set kernel radius desingularized source
myR0 = 3e-3;
% Markstein length
my_lM = 1e-4;
% Number of sources per milimeter
SpMM = 10;
PpMM = 2;
% Distance between sources and G=0
dispXi = 0*5e-4;

% Integration settings
t_max = 0.03;

% Consider impact of curvature on sL (for heat relase and source strengths)?
curvContri_sL = 1;
% Spatial resolution Xi
Nx1F = 200;


%% Setup
% Set default values for parameters
[ p ] = setUpPredefinedFlame( 'Axel_DuctFlame' );

% Get mapping infos
[ myMap ] = return_SCmap_SCFT( p );

% Wrap source info to array
sourceInfo(1) = dispXi;
sourceInfo(2) = myR0;
if strcmpi(method_I,'panel')
  SpMM = PpMM;
end
sourceInfo(3) = SpMM;
 
% Set up mean
[ myMean , myN , x1F , dx1F , ind_sources ] = getMeanFlameFlow( p , Nx1F , do_I , method_I , sourceInfo);

% Grid for test plot
[ Xi_all , X_up_mesh , X_down_mesh ] = make_grid4plot( p );

% Remove source at -infintiy
myMean.vel.sourceDat.G = [];
myMean.vel.sourceDat.x = [];
myMean.vel.sourceDat.xi = [];


%% Add a vortex
vStrength = 0.03;
myR0 = 1e-3;
% center
myMean.vel.vortDat.G = vStrength; 
myMean.vel.vortDat.x = 1i*p.R_i/2; 
myMean.vel.vortDat.xi = myMap.xi_x(myMean.vel.vortDat.x);
myMean.vel.vortDat.r0 = myR0; 

% upstream
Nvort0 = length(myMean.vel.vortDat.G);
myMean.vel.vortDat.G = [ myMean.vel.vortDat.G ; vStrength ]; 
myMean.vel.vortDat.x = [ myMean.vel.vortDat.x ; -p.R_i*0.7+1i*p.R_i/2 ]; 
myMean.vel.vortDat.xi = [ myMean.vel.vortDat.xi ; myMap.xi_x(myMean.vel.vortDat.x(Nvort0+1:end)) ]; 
myMean.vel.vortDat.r0 = [ myMean.vel.vortDat.r0 ; myR0 ]; 

% downstream
Nvort0 = length(myMean.vel.vortDat.G);
myMean.vel.vortDat.G = [ myMean.vel.vortDat.G ; vStrength ]; 
myMean.vel.vortDat.x = [ myMean.vel.vortDat.x ; p.R_i*0.7+1i*p.R_i/2 ]; 
myMean.vel.vortDat.xi = [ myMean.vel.vortDat.xi ; myMap.xi_x(myMean.vel.vortDat.x(Nvort0+1:end)) ]; 
myMean.vel.vortDat.r0 = [ myMean.vel.vortDat.r0 ; myR0 ];


%% Add a source
sStrength = 0.03;
% % center
% Nvort0 = length(myMean.vel.sourceDat.G);
% myMean.vel.sourceDat.G = [ myMean.vel.sourceDat.G ; sStrength ]; 
% myMean.vel.sourceDat.x = [ myMean.vel.sourceDat.x ; 1i*p.R_i/2 ]; 
% myMean.vel.sourceDat.xi = [ myMean.vel.sourceDat.xi ; myMap.xi_x(myMean.vel.sourceDat.x(Nvort0+1:end)) ]; 
% 
% % upstream
% Nvort0 = length(myMean.vel.sourceDat.G);
% myMean.vel.sourceDat.G = [ myMean.vel.sourceDat.G ; sStrength ]; 
% myMean.vel.sourceDat.x = [ myMean.vel.sourceDat.x ; -p.R_i*0.7+1i*p.R_i/2 ]; 
% myMean.vel.sourceDat.xi = [ myMean.vel.sourceDat.xi ; myMap.xi_x(myMean.vel.sourceDat.x(Nvort0+1:end)) ]; 
% 
% % downstream
% Nvort0 = length(myMean.vel.sourceDat.G);
% myMean.vel.sourceDat.G = [ myMean.vel.sourceDat.G ; sStrength ]; 
% myMean.vel.sourceDat.x = [ myMean.vel.sourceDat.x ; p.R_i*0.7+1i*p.R_i/2 ]; 
% myMean.vel.sourceDat.xi = [ myMean.vel.sourceDat.xi ; myMap.xi_x(myMean.vel.sourceDat.x(Nvort0+1:end)) ]; 


%% Plot
% test Plot mean
x1L_flame = myMean.flcoord(1,:) + 1i*myMean.flcoord(2,:);
% Plot sources
plotSources = 'sourceL';
% Evaluate and visualize velocity field
% u1
myF1 = figure('Position',[0 50 800 550],'Color','w');hold on;
draw_FlowField( X_up_mesh , X_down_mesh , myMean.vel , x1L_flame , p , 'scaleU' , 0.001 ,...
  plotSources , ind_sources,'Xi',Xi_all,'plotField','u1','doQuiver',0.005,'fig',myF1,'cLim',[-1 1]*0.4);

% u2
myF1 = figure('Position',[820 50 800 550],'Color','w');hold on;
draw_FlowField( X_up_mesh , X_down_mesh , myMean.vel , x1L_flame , p , 'scaleU' , 0.001 ,...
  plotSources , ind_sources,'Xi',Xi_all,'plotField','u2','doQuiver',0.005,'fig',myF1,'cLim',[-1 1]*0.4);