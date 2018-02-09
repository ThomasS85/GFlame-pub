function [ A_surf , data ] = calcSurfaceArea( data , grid , schemeData , p )
%CALCSURFACEAREA Integrates function along the flame surface from 2D data and geometry
% information (slit or conical)
%
% Inputs:
%   - data          : G-field data (matrix)
%   - grid          : Grid struct
%   - schemeData    : Struct with all scheme data including a struct...
%       .solver                 : Settings for delta function
%   - p             : Struct with general settings (flame and domain type)
%
% Outputs:
%   - A_surf        : Flame surface area [m^2] or heat release []
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


%% Compute flame surface area or heat release?
if strcmp(schemeData.solver.outputY,'area')
  % Compute surface area
  funI = 1;
  
else
  % Compute heat release
  
  % Equivalence ratio (for ER perturbation call function which evaluates Phi-field here!)
  Phi = p.Phi_m;
  
  % Get factors for integration
  %   (1) Density [kg/m^3]
  [ ~ , rho_in ] = getMoleMassFromPhi( p.mixture , Phi );
  
  %   (2) specific heat of reaction [J/kg]
  [ delta_h_R ] = getHeatOfReaction( p.mixture , Phi );
  
  %   (3) Flame speed: From equivalance ratio, Temperature and pressure [m/s]
  [ s_l_u ] = getLaminarFlameSpeed( p.mixture , Phi , p.T_in , p.p_in );
  
  %   (4) Flame speed: From curvature [m/s]
  if schemeData.solver.curvature == 1 
    % Find index of convective schemeData
    tmp = cellfun(@(x) strcmp(func2str(x),'termCurvature') , schemeData.innerFunc );
    indSD = find(tmp,1,'first');
    % Compute flame curvature 
    [curF,~] = feval(schemeData.innerData{indSD,1}.curvatureFunc,grid,data);
    % Apply Gauss filter to smooth data
    if schemeData.solver.sigmaG > 0
      curF = imgaussfilt(curF,schemeData.solver.sigmaG);
    end
    % Compute flame speed from curvature
    s_l_u = s_l_u .* ( 1 + p.marksteinLength * curF );
    % Debug plot
    % visualizeLevelSet(grid, s_l_u, 'surf')
    
  end
   
  % Now compute function which should be integrated along the flame surface
  %  (scalar or matrix of size(data)
  funI = rho_in .* delta_h_R .* s_l_u;
  
end


%% Perform line integral

if strcmp(p.flameType,'slit')
  % SLIT burner (extension in x3 direction is assumed to be 1)
  
  % Surface is length of line
  if strcmpi(schemeData.solver.SurfInt_method,'deltaFun')
    % Use delta function (1st or 2nd order)
    [integralLine , data] = lineIntegral( grid, data , 1*funI , schemeData.solver );
    
  elseif strcmpi(schemeData.solver.SurfInt_method,'contour')
    % Use integration of contour line
    [integralLine , data] = lineInegral_contour( grid, data , 1*funI , schemeData.solver );
    
  end
  
  if strcmp(p.domainType,'sym')
    % Symmetric geometry
    A_surf  = 2 * integralLine;
  else
    % Full geometry
    A_surf  = integralLine;
  end
  
elseif strcmp(p.flameType,'conical')
  % CONICAL burner
  
  % Surface is weighted by radius (x2 coordinate)
  if strcmpi(schemeData.solver.SurfInt_method,'deltaFun')
    % Use delta function (1st or 2nd order)
    [integralLine , data] = lineIntegral( grid, data , funI.*abs( grid.xs{2,1} ) , schemeData.solver );
    
  elseif strcmpi(schemeData.solver.SurfInt_method,'contour')
    % Use integration of contour line
    [integralLine , data] = lineInegral_contour( grid, data , funI.*abs( grid.xs{2,1} ) , schemeData.solver );
    
  end
  
  if strcmp(p.domainType,'sym')
    % Symmetric geometry (radius r = x2 )
    A_surf = 2 * pi * integralLine;
  else
    % Full geometry (radius r = x2 )
    A_surf = pi * integralLine;
  end
  
else
  error('Unkonwn flame geoemtry!')
  
end


end

