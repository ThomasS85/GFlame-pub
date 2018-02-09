function [ data ] = setInitialConditions( myGrid , p , solver , makeSignedDist )
%SETINITIALCONDITIONS Sets default initial conditions for flame

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% In order to respect flame lift-off, take the reference position at flame base
myRef = solver.x1Lim(1);

if strcmp(solver.initial,'default')
  if strcmp(solver.holdFlame,'radial') && ~(strcmpi(p.geom,'V')||strcmpi(p.geom,'M'))
    % creats a line for initial conditions with slope of 1 in x1-direction which is flame normal direction
    % (signed distance function)
    %     m = 1;   % G<0 : fresh gas
    %     offSet = p.H_flame * 0.01;
    %     t = - ( p.liftOff + offSet) * m;  % offset needed to stabilise flame for liftOff=0!
    %     data =  m * grid.xs{1} + t;
    
    %     % Creates skewd line as initial condition
    %     m = -tan(p.alpha);
    %     x2_0 = solver.x2Lim(2);
    %     offSet = p.H_flame * 0.01;
    %     data = grid.xs{2} - m * ( grid.xs{1} - offSet ) - x2_0;
    
    % Creates parabula for initial condition. Slope normal to flame front is not signed distance!
    const = solver.x1Lim(2)*0.5 ;
    prefactor = -const / solver.x2Lim(2)^2;
    offSet = p.H_flame * 0.01;
    data =  myGrid.xs{1}-myRef - prefactor*myGrid.xs{2}.^2 - const - offSet;
    
    
  elseif strcmp(solver.holdFlame,'axial') || (strcmpi(p.geom,'V')||strcmpi(p.geom,'M'))
    if strcmpi(p.geom,'Vinv')
      % Creates parabula for initial condition. Slope normal to flame front is not signed distance!
      const = solver.x1Lim(2)*0.5 ;
      prefactor = -const / p.R_flame^2;
      data =  myGrid.xs{1}-myRef - prefactor*myGrid.xs{2}.^2 - const;
      
      %     % Creates skewd line as initial condition
      %     m = -tan(p.alpha);
      %     x2_0 = p.R_flame;
      %
      %     if strcmp(p.domainType,'sym')
      %       % Symmetric domain
      %       data = grid.xs{2}(grid.xs{2}>=0) - m * grid.xs{1}(grid.xs{2}>=0) - x2_0;
      %
      %     elseif strcmp(p.domainType,'full')
      %       % Full domain
      %
      %       % Initialize data
      %       data = zeros(grid.shape);
      %       % x2>=0: upper part
      %       data(grid.xs{2}>=0) = grid.xs{2}(grid.xs{2}>=0) - m * grid.xs{1}(grid.xs{2}>=0) - x2_0;
      %       % x2<0:  lower part
      %       data(grid.xs{2}<0) =  -grid.xs{2}(grid.xs{2}<0) - m * grid.xs{1}(grid.xs{2}<0) - x2_0 ;
      %     end
      
      
    elseif strcmpi(p.geom,'V')
      % Creates skewd line as initial condition
      m = tan(p.alpha);
      if strcmp(p.domainType,'sym')
        % Symmetric domain
        data = -(myGrid.xs{2}(myGrid.xs{2}>=0) - m * (myGrid.xs{1}(myGrid.xs{2}>=0)-myRef) - p.R_rod);
        
      elseif strcmp(p.domainType,'full')
        % Full domain
        
        % Initialize data
        data = zeros(myGrid.shape);
        % x2>=0: upper part
        data(myGrid.xs{2}>=0) = -(myGrid.xs{2}(myGrid.xs{2}>=0) - m * (myGrid.xs{1}(myGrid.xs{2}>=0)-myRef) - p.R_rod);
        % x2<0:  lower part
        data(myGrid.xs{2}<0) = -(-myGrid.xs{2}(myGrid.xs{2}<0) - m * (myGrid.xs{1}(myGrid.xs{2}<0)-myRef) - p.R_rod);
      end
      % reshape vector to matrix
      data = reshape(data,myGrid.N');
      
      
    elseif strcmpi(p.geom,'M')
      % Create parabula as initial condition
      myConst = (p.R_i-p.R_rod);
      myFac = solver.x1Lim(2)*0.5 ;
      myLower = p.R_rod - 5*myGrid.dx(1);
      myupper = p.R_i + 5*myGrid.dx(1);
      if strcmp(p.domainType,'sym')
        % Symmetric domain
        myCond = myGrid.xs{2}>myLower & myGrid.xs{2}<myupper;
        data = zeros(myGrid.shape);
        data(myCond) = myGrid.xs{1}(myCond)-myRef - myFac*sin( pi/myConst * (myGrid.xs{2}(myCond) - p.R_rod) );
      elseif strcmp(p.domainType,'full')
        % Full domain
        myCondU = myGrid.xs{2}>myLower & myGrid.xs{2}<myupper;
        myCondL = -myGrid.xs{2}>myLower & -myGrid.xs{2}<myupper;
        data = zeros(myGrid.shape);
        data(myCondU) = myGrid.xs{1}(myCondU)-myRef - myFac*sin( pi/myConst * (myGrid.xs{2}(myCondU) - p.R_rod) );
        data(myCondL) = myGrid.xs{1}(myCondL)-myRef - myFac*sin( pi/myConst * (-myGrid.xs{2}(myCondL) - p.R_rod) );
      end
      
    end
    
  end
  
  % Initialite G-field as signed distance function
  if makeSignedDist
    % Initialise data: Make sure G is a signed distance function
    disp('Initialising G-field as signed distance function...')
    % Set temporalily extrapolation boundary conditions (otherwise the default periodic BC might lead to bad
    % results of reinitialisation close to boundaries!)
    myGrid.bdry{1,1} = @addGhostExtrapolate;
    myGrid.bdry{2,1} = @addGhostExtrapolate;
    
    % Perform reinitialisation (improves quality of boundary condition)
    data = signedDistanceIterativeSubCellFix(myGrid, data, 'veryHigh', max(myGrid.dx(1) , myGrid.dx(2))*30, 1e-4);
    disp('G-field initialised!')
  else
    warning('Initial G-field was not initialised to a signed distance function!')
  end
  
elseif strcmp(solver.initial,'useInit')
  % uses a init.mat file in the root directory for initial conditions
  load('init.mat','data')
  
  % Check if data has the same size as the grid
  [z,s] = size(data);
  
  if (myGrid.N(1) ~= z) || (myGrid.N(2) ~=s)
    error('Provided data for initial conditions is not compatible with grid!')
  end
  
else
  % Sets the whole -field to -1
  data =  -1 * myGrid.xs{1};
  
end


end

