function [  schemeData   ] = initialize_FPB(time,data, schemeData )
%INITIALIZE_FPB initializes the velocity field parameters for velocityFieldFirstPrincipleBased
%
% Inputs:
%   - data        - G-Field matrix
%
%   - schemeData  - Struct which contains information about all nescessary
%                   parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  by Axel Zimmermann (07.2018)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID MAPPING     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAP the grid
grid = map_grid(schemeData.grid.xs{1}+0.0000001,schemeData.grid.xs{2},schemeData.p );
%WRITE the grid in schemeData
schemeData.FPB.grid = grid;

%map grid for the plotting if requested
if strcmp(schemeData.FPB.plot.plot_physical,'y')
  %Define new grid, in the entire domain
  solver.x1Lim = [schemeData.grid.min(1) schemeData.grid.max(1)];
  solver.x2Lim = [schemeData.grid.min(2) schemeData.p.R_a ];
  solver.dim = schemeData.grid.dim;
  solver.dxi = schemeData.grid.dx ;
  [ schemeData.FPB.plot.myGrid ] = generateGrid( solver );
  
  %Map new grid to the image domain
  schemeData.FPB.plot.grid = map_grid(schemeData.FPB.plot.myGrid.xs{1}+0.0000001,...
    schemeData.FPB.plot.myGrid.xs{2},schemeData.p );
  
  %save additional data, necessary for the velocity calculation in the
  %plotting programm
  schemeData.FPB.plot.R_i = schemeData.p.R_i;
  [schemeData.FPB.plot.myMap] = return_SCmap_SCFT( schemeData.p );
  schemeData.FPB.plot.meanFlowSpeed = schemeData.meanFlowSpeed ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION ABOUT THE SOURCES    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%source RADIUS  %
%%%%%%%%%%%%%%%%%
%reading defined source radius
velPar.sourceDat.kw = schemeData.FPB.source_kw;
velPar.sourceDat.myR = schemeData.FPB.source_myR;

switch velPar.sourceDat.kw
  case 'no_radius'
    velPar.sourceDat.radius_s = 0;
  case 'small'
    velPar.sourceDat.radius_s = 0.00025;
  case 'normal'
    velPar.sourceDat.radius_s = 0.0005;
  case 'high'
    velPar.sourceDat.radius_s = 0.0008;
  case 'my_radius'
    velPar.sourceDat.radius_s = schemeData.FPB.source_myR;
  otherwise
    error('source radius is not defined!')
    
end

%%%%%%%%%%%%%%%%%%%%
%AMOUNT of sources %
%%%%%%%%%%%%%%%%%%%%
switch schemeData.FPB.source_amount
  case 'high'
    velPar.sourceDat.amount = 200;
  case 'normal'
    velPar.sourceDat.amount = 150;
  case 'low'
    velPar.sourceDat.amount = 100;
  case 'my_amount'
    velPar.sourceDat.amount = schemeData.FPB.source_myA;
  otherwise
    error('number of sources is not defined!')
end

%%%%%%%%%%%%%%%%%%%
%source POSITION  %
%%%%%%%%%%%%%%%%%%%
% extract SOURCE position
[ C ] = extract_IsoLine_FPB( schemeData.grid, data ,velPar.sourceDat.radius_s);

% interpolate position of sources to an equidistant distance
[velPar.sourceDat.x,~]= interp2Dcurve_2equidistantGrid_FPB([C(1,2:end);C(2,2:end)] ,velPar.sourceDat.amount);
velPar.sourceDat.x = transpose(velPar.sourceDat.x(:,1)+ 1i*velPar.sourceDat.x(:,2));

%map source coordinates to image domain
velPar.sourceDat.xi = SCmapInv_SCFT( velPar.sourceDat.x , schemeData.p , 'L1' );

%check if mapping was successfull
velPar.sourceDat.x_check = SCmap_SCFT( velPar.sourceDat.xi , schemeData.p , 'L1' );
if find((velPar.sourceDat.x-velPar.sourceDat.x_check)>1e-7,1)
  error(' source mapping was not successful')
end

velPar.sourceDat.G = zeros(size(velPar.sourceDat.x));

%MIRROR for image domain
velPar.doMirror ='doMirror';

% Source Radius
velPar.sourceDat.r0 = zeros(size(velPar.sourceDat.x))+velPar.sourceDat.radius_s;


%% mean flow source
velPar.sourceDat.x(end)=-inf;
velPar.sourceDat.xi(end)=0;
velPar.sourceDat.r0(end)=0;
velPar.sourceDat.G(end) = schemeData.p.R_i*(schemeData.meanFlowSpeed);
velPar.u_p=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLAME POSITION EXTRACTION                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(schemeData.FPB.do_proceed,'y')
  % extract flame coordinates from G-field
  [ C2 ] = extract_IsoLine_FPB( schemeData.grid, data ,0);
  [C2,~]= interp2Dcurve_2equidistantGrid_FPB([C2(1,2:end);C2(2,2:end)] ,150);
  C2 = transpose(C2);
  Flame.coor{1, 1} =C2(:,1:end);
  
  % complex flame coordinatesd in physical domain
  Flame.x = (Flame.coor{1, 1}(1,:)+1i*Flame.coor{1, 1}(2,:));
  
  %complex flame coordinates in image domain
  Flame.xi = SCmapInv_SCFT( Flame.x , schemeData.p , 'L1');
  schemeData.FPB.Flame = Flame;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KUTTA-CONDITION         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define position of the Kutta panel in both domains
velPar.Kutta.position_Kutta_x = 0+1i*(schemeData.p.R_i);
velPar.Kutta.position_Kutta_xi = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION ABOUT THE VORTEXES       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% should shear layer be computed ?
if strcmp(schemeData.FPB.shear_layer,'do_sl') && strcmpi(schemeData.p.CombType,'backwardFacingStep')
  
  % Evaluate velocity and potential at Kutta point xiKutta for the
  % vortex strength
  velPar.vortDat={};
  myMap=return_SCmap_SCFT( schemeData.p );
  [ ~ , dOmega_dxi_xiKutta_sum ] = evalFlowField_PointSource( velPar.Kutta.position_Kutta_xi , velPar , 'myMapping' , myMap );
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INITIALIZE SHEAR LAYER   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Define shear layer angle and length in image domain
  MyBeta = pi/9;
  H = abs( myMap.dxi_dx( myMap.xi_x( myMap.x_xi(velPar.Kutta.position_Kutta_xi) + myMap.l_ref * 0.6 ) ) ) * myMap.l_ref*0.5;
  
  %get a first approximation of the shear layer
  velPar.vortDat.xi=linspace (1+H*exp(1i*MyBeta),1+H*exp(1i*MyBeta)*1700,3000);
  velPar.vortDat.x = SCmap_SCFT( velPar.vortDat.xi , schemeData.p , 'L1' );
  
  % now interpolate the vortex location in the physical domain
  [velPar.vortDat.x,~]= interp2Dcurve_2equidistantGrid_FPB([real(velPar.vortDat.x);imag(velPar.vortDat.x)] ,100);
  velPar.vortDat.x = transpose(velPar.vortDat.x(:,1)+ 1i*velPar.vortDat.x(:,2));
  
  % now map vortex position in the image domain once again
  velPar.vortDat.xi = SCmapInv_SCFT( velPar.vortDat.x , schemeData.p , 'L1' );
  
  % compute strength of Kutta panel
  Kutta_G = dOmega_dxi_xiKutta_sum * pi /...
    ( 1i * sqrt(H ) * ( exp(-1i*MyBeta) -  exp(1i*MyBeta) ) );
  
  %distance between the vortexes for the strength
  distance = abs(velPar.vortDat.x(1:end-1)-velPar.vortDat.x(2:end));
  distance = [0,distance,0];
  
  % Strength of each vortex
  % (ds*Kutta_panel_strength/Kutta_panel_length)
  velPar.vortDat.G = (distance(1:end-1)+distance(2:end))/2*Kutta_G/abs(velPar.vortDat.x(1)-schemeData.p.R_i*1i);
  
  
  %vortex RADIUS
  %initializing vortex radius
  velPar.vortDat.r0  = zeros(length(velPar.vortDat.x),1);
  %reading defined radius for vortex
  velPar.vortDat.kw = schemeData.FPB.vortex_kw;
  velPar.vortDat.myR = schemeData.FPB.vortex_myR;
  %radius for vortex center
  switch velPar.vortDat.kw
    case 'no_radius'
      radius_v = 0;
    case 'small'
      radius_v = 0.0001;
    case 'normal'
      radius_v = 0.0005;
    case 'high'
      radius_v = 0.001;
    case 'my_radius'
      radius_v = velPar.vortDat.myR;
    otherwise
      error('kernel width of vortex not defined')
  end
  velPar.vortDat.r0 = velPar.vortDat.r0 + radius_v;
  
  
else
  %%%%%%%%%%%%%%%%%%%%
  % NO SHEAR LAYER   %
  %%%%%%%%%%%%%%%%%%%%
  velPar.vortDat = {};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARING FOR OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVING source and vortex information
schemeData.FPB.velPar = velPar;


end

