function [ vel , schemeData , u_ref ] = velocityFieldFirstPrincipleBased(time, data, schemeData )
% velocityFieldFirstPrincipleBased evaluates the velocity field due to bulk stream, flame influence and vortices
%
% Inputs:
%
%   - t           - actual time
%   - data        - G-Field matrix
%
%   - schemeData  - Struct which contains information about all nescessary
%                   parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%  by Axel Zimmermann (06.2018)  %
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE VELOCITY FIELD IF NECESSARY                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isfield(schemeData.FPB,'velPar'))
  schemeData  = initialize_FPB(time, data, schemeData );
end
velPar=schemeData.FPB.velPar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION ABOUT THE VORTEXES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% should shear layer be computed ?
if strcmp(schemeData.FPB.shear_layer,'do_sl')
  % UNWRAP existing vortex parameters into velPar
  velPar.vortDat = schemeData.FPB.velPar.vortDat;
else
  %do not compute shear layer
  velPar.vortDat = {};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION ABOUT THE SOURCES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract amount of sources
velPar.sourceDat.amount= schemeData.FPB.velPar.sourceDat.amount;


%%%%%%%%%%%%%%%%%%
%source POSITION 
%%%%%%%%%%%%%%%%%%
% extract source position
[ C ] = extract_IsoLine_FPB( schemeData.grid, data ,schemeData.FPB.velPar.sourceDat.radius_s);

% interpolate position of flame front to equidistant distance
[velPar.sourceDat.x,~]= interp2Dcurve_2equidistantGrid_FPB([C(1,2:end);C(2,2:end)] ,velPar.sourceDat.amount);
velPar.sourceDat.x = transpose(velPar.sourceDat.x(:,1)+ 1i*velPar.sourceDat.x(:,2));

%complex source coordinates in image domain
velPar.sourceDat.xi = SCmapInv_SCFT( velPar.sourceDat.x , schemeData.p , 'L1','xi0',transpose(schemeData.FPB.velPar.sourceDat.xi(1:end-1)) );

%%%%%%%%%%%%%%%%%%
%source STRENGTH 
%%%%%%%%%%%%%%%%%%
% distance between sources
velPar.sourceDat.source_distance = abs(velPar.sourceDat.x(1:end-1)-velPar.sourceDat.x(2:end));

if strcmpi(schemeData.p.CombType,'backwardFacingStep')
  velPar.sourceDat.source_distance = [real(velPar.sourceDat.x(1)),velPar.sourceDat.source_distance,imag(velPar.sourceDat.x(end))];%  axial
elseif strcmpi(schemeData.p.CombType,'duct')
  velPar.sourceDat.source_distance = [schemeData.p.R_i-imag(velPar.sourceDat.x(1)),velPar.sourceDat.source_distance,imag(velPar.sourceDat.x(end))];%  axial
end

% curvature calculation for a varying s_L
if strcmp(schemeData.FPB.do_source_curvature,'y')
  kappa = LineCurvature2D([real(transpose(velPar.sourceDat.x)) transpose(imag(velPar.sourceDat.x))]);
  kappa = smoothdata (kappa ,'lowess',30);
  velPar.s_l_u = schemeData.p.s_l_u*(1-kappa*schemeData.p.marksteinLength);
else
  velPar.s_l_u = zeros(size(velPar.sourceDat.x))'+schemeData.p.s_l_u ;
end

%strength of sources
velPar.sourceDat.G = (velPar.sourceDat.source_distance(1:end-1)+velPar.sourceDat.source_distance(2:end))...
  /2.*schemeData.p.s_l_u'*(schemeData.p.E-1);

%strength of first sources equals 0, compare CFD data
if strcmpi(schemeData.p.CombType,'backwardFacingStep')
  velPar.sourceDat.G(1:2)=0;
end


%MIRROR for image domain
velPar.doMirror ='doMirror';

%%%%%%%%%%%%%%%%%%
%source RADIUS  
%%%%%%%%%%%%%%%%%%
velPar.sourceDat.radius_s = schemeData.FPB.velPar.sourceDat.radius_s;
velPar.sourceDat.r0 = zeros(size(velPar.sourceDat.x))+velPar.sourceDat.radius_s;


%% mean flow source
%Get distributed distortion SIGNAL
if strcmp(schemeData.type,'transient')
  [ ~ , u1s_ref ] = applyKernel_convectiveMods( data , schemeData , time );
else
  u1s_ref = 0;
end

%adds mean flow source
velPar.sourceDat.x(end+1)=-inf;
velPar.sourceDat.xi(end+1)=0;
velPar.sourceDat.r0(end+1)=0;
velPar.sourceDat.G(end+1) = schemeData.p.R_i*(u1s_ref+schemeData.meanFlowSpeed);
velPar.u_p = schemeData.meanFlowSpeed+u1s_ref;

% save acoustic signal for plotting
if    strcmp(schemeData.FPB.plot.plot_physical,'y')
  schemeData.FPB.plot.ac= u1s_ref;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLAME POSITION       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(schemeData.FPB.do_proceed,'y')
  %%%%%%%%%%%%%%%%%%
  % Proceed flame velocity normal to flame front in the domain  
  %%%%%%%%%%%%%%%%%%
  % extract flame coordinates from G-field
  [ C2 ] = extract_IsoLine_FPB( schemeData.grid, data ,0);
  [C2,~]= interp2Dcurve_2equidistantGrid_FPB([C2(1,2:end);C2(2,2:end)] ,150);
  C2 = transpose(C2);
  Flame.coor{1, 1} =C2(:,1:end);
  
  % complex flame coordinatesd in physical domain
  Flame.x = (Flame.coor{1, 1}(1,:)+1i*Flame.coor{1, 1}(2,:));
  
  %complex flame coordinates in image domain
  Flame.xi = SCmapInv_SCFT( Flame.x , schemeData.p , 'L1','xi0',transpose(schemeData.FPB.Flame.xi));
  
  %velocity should be evaluated only at the flame front
  eval_vel.x = Flame.x;
  eval_vel.xi = Flame.xi;
  
else
  
  %%%%%%%%%%%%%%%%%%
  % evaluate velocity in the whole domain  
  %%%%%%%%%%%%%%%%%%
  eval_vel.x = schemeData.FPB.grid.x;
  eval_vel.xi = schemeData.FPB.grid.xi;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE VELOCITY  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Struct with mapping information
[myMap] = return_SCmap_SCFT( schemeData.p );

if strcmpi(schemeData.p.CombType,'backwardFacingStep')
  [ u_ges , velPar , ~ ] = evalFlowField_physicalDomain(...
    eval_vel.xi , velPar , myMap,'xi','L1','doKutta',velPar.Kutta.position_Kutta_xi,'noPanel2Point');
  
elseif strcmpi(schemeData.p.CombType,'duct')
  [ u_ges , velPar , ~ ] = evalFlowField_physicalDomain(...
    eval_vel.xi , velPar , myMap,'xi','L1','noKutta','noPanel2Point');
  u_ges= conj(u_ges);
  
end

%velocity FIELD
vel{1}=real(u_ges);
vel{2}=imag(u_ges);

%%%%%%%%%%%%%%%%%%
%PROCEED VELOCITY normal to flame front in whole domain
%%%%%%%%%%%%%%%%%%
if strcmp(schemeData.FPB.do_proceed,'y')
  clear vel
  
  valC{1, 1} = transpose(real(u_ges));
  valC{2, 1} = transpose(imag(u_ges));
  
  if strcmpi(schemeData.p.CombType,'backwardFacingStep')
    %change boundary conditions
    grid = schemeData.grid;
    grid.bdryData{2, 1}=grid.bdryData{1, 1};
    grid.bdry = {@addGhostExtrapolate ; @addGhostExtrapolate};
    %proceed the velocity
    [ vel ] = expandVelocityFromZero( Flame.coor , valC , grid , data );
    
  elseif strcmpi(schemeData.p.CombType,'duct')
    %proceed the velocity
    [ vel ] = expandVelocityFromZero( Flame.coor , valC , schemeData.grid , data );
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARING FOR OUTPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INLET velocity
u_ref = velPar.u_p;

%SAVING source and vortex information
schemeData.FPB.velPar = velPar;

if strcmp(schemeData.FPB.do_proceed,'y')
  schemeData.FPB.Flame = Flame;
end
%pcolor(real(schemeData.FPB.grid.x),imag(schemeData.FPB.grid.x),vel{1, 1})
%shading interp

end

