function [ myGrid ] = setBoundaryConditions( myGrid , p , solver , data )
%SETBOUNDARYCONDITIONS Sets boundary conditions for G-equation

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


myGrid.bdry = { solver.boundaryCon{1}; solver.boundaryCon{2} };
myGrid.bdryData = { [];  [] };

if strcmp(solver.holdFlame,'radial') && ~(strcmpi(p.geom,'V')||strcmpi(p.geom,'M'))
  if strcmp(p.domainType,'sym') 
    % Symmetric domain
    myGrid.bdryData{2}.symmetryPlane = 'lower';
    myGrid.bdryData{2}.upperVector = data(:,end);
    myGrid.bdryData{2}.towardZero = 0;

  else
    % Full domain
    myGrid.bdryData{2}.symmetryPlane = 'lower';
    myGrid.bdryData{2}.lowerVector = data(:,1);
    myGrid.bdryData{2}.upperVector = data(:,end);
    myGrid.bdryData{2}.towardZero = 0;

  end

elseif strcmp(solver.holdFlame,'axial') || (strcmpi(p.geom,'V')||strcmpi(p.geom,'M'))
  % same for symmetric or full
  myGrid.bdryData{1}.lowerVector = data(1,:);
  myGrid.bdryData{1}.towardZero = 0;
end

end

