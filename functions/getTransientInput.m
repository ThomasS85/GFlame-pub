function [ u_vec , schemeData ] = getTransientInput( schemeData , t_vec , ~ )
%TRANSIENTINPUT Function returns excitation signal u_vec (between -1 and 1) for a given time vector t_vec from
%   different sources: file, harmonic, step or TaX
%
% All Functions must be able to return results for all time (if needed extrapolate!)
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 09.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


if strcmp(schemeData.transientDef,'file')
  % From file (needs path2File). Amplitude in file should be 1  (from -1 to 1)
  % Everything before the beginning of the signal is zero
  u_vec = readNInterpolateDatafromFile( t_vec , schemeData.path2File , 'extrap' , 0 );
  
  
elseif strcmp(schemeData.transientDef , 'harmonic' )
  % transient v: sin function (needs omega and t_transient0)
  % transient part starts after this time
  t_on = schemeData.t_transient0;
  % Calculate output. Everything before t_on is 0
  u_vec = sin( schemeData.omega * (t_vec-t_on) ) .* myHeaviside2(t_vec-t_on);
  

elseif strcmp(schemeData.transientDef , 'sweep' ) || strcmp(schemeData.transientDef , 'broad' )
  % From file (needs path2File). Amplitude in file should be 1  (from -1 to 1)
  % Everything before the beginning of the signal is zero
  u_vec = readNInterpolateDatafromFile( t_vec , schemeData.path2File , 'extrap' , 0 );
  
  
elseif strcmp(schemeData.transientDef , 'step' )
  % Transient behaviour: Step function (needs sharpness and t_transient0)
  % transient part starts after this time
  t_on = schemeData.t_transient0;
  u_vec = ( atan( ( t_vec - t_on )*schemeData.sharpness ) / pi + 0.5 );
  
  
elseif strcmpi(schemeData.transientDef,'TaX')
  % from Tax.
  
  % To do
  %   - use same time integration scheme as GFLAME: For full time step take linear combination of intermediate
  %     results (save those + check if odeCFL1,2 or three is beeing used)
  %   - Clean up: Put some commands to external functions (if everything works)
  %   - Verification!
  %   - Check if flame surface evaluation is not too noisy!
    
  % If this is the first time step, return initial uref data  
  t_now = max(t_vec); % Get actual time
  if schemeData.t_vec(1) < t_now
    
    % (1) Create time vector sampled by schemeData.Ts from t_vec(1) to t_now    
    N = ceil( ( t_now - schemeData.t_vec(1) ) / schemeData.Ts_ac ) + 1;
    t_sim_vec = linspace(schemeData.t_vec(1),t_now,N);
    Ts = t_sim_vec(2) - t_sim_vec(1); % new Ts
    
    % (2) Fit of polynomal to past flame surface data and use this to extrapolate A'/A_m
    polyOrder = 2;                    % 1=linear , 2 = quadratic , 3 = cubic , ...
    dataLength = (polyOrder+1) * 3;   % Time series 3x as much as necessary (smoothing)
    if length(unique(schemeData.t_vec)) < dataLength
      % Switch back to zero order polynomial fit if not enough data has been computed yet
      polyOrder = 0;                    
      dataLength = (polyOrder+1) * 2;
    end
    % Build matrix for least square
    A_poly_mat = zeros(dataLength,polyOrder+1);
    for i=1:polyOrder+1
      A_poly_mat(:,polyOrder-i+2) = schemeData.t_vec(1:dataLength).^(i-1);
    end
    % Least squaree fit: Get coefficients of polynomial
    polyCoeff = A_poly_mat \ schemeData.A_surf(1:dataLength);
    
    % Evaluate fit to surface area/ heat release for acoustic simulation time steps
    A_surf_sim = polyval(polyCoeff,t_sim_vec)';
    % Get normalized area fluctuation
    A_relFluc_sim = (A_surf_sim - schemeData.A_mean) / schemeData.A_mean;
    
    % debugging: 
%     f=figure; plot(t_sim_vec,A_surf_sim);hold on; plot(schemeData.t_vec(1:dataLength),schemeData.A_surf(1:dataLength),'o');...
%     legend(['extrapolation fit (order ',num2str(polyOrder),')'],'fit data','Location','NorthWest'); xlabel('t');...
%     ylabel('A_{surf}')
%     close(f);
    
    
    % (3) Simulation of TaX model from last time step to now.
    [ u_rel_sim , x ] = schemeData.sys_acoustic.lsim( A_relFluc_sim , Ts , 'RK4' );
    
    % Update system with this solution
    schemeData.sys_acoustic.x0 = x;
    
    % Append solution to uref and t_vec
    uref_sim = [ u_rel_sim(end) ; (schemeData.uref-schemeData.meanFlowSpeed) / schemeData.meanFlowSpeed ];
    t_vec_sim = [ t_now ; schemeData.t_vec ];
    
    
    % (4) Construct output u_vec by interpolation
    % Get unique elements in schemeData.t_vec (in the beginng everything might be equal!)
    [ t_vec_sim_unique , ia ] = unique(t_vec_sim);
    % Interpolate
    u_vec = interp1( t_vec_sim_unique ,uref_sim(ia) , t_vec , 'linear' , 0 );
    
  else
    % First time step and initialisation just interpolate given data
    % Get unique elements in schemeData.t_vec (in the beginng everything might be equal!)
    [ ~ , ia ] = unique(schemeData.t_vec);
    if length(ia) > 1
      % Interpolate
%       uref = interp1( schemeData.t_vec(ia) , schemeData.uref(ia) , t_vec , 'linear' , 0 );
      uref = interp1( schemeData.t_vec(ia) , schemeData.uref(ia) , t_vec , 'linear' , 'extrap' );
      u_vec = ( uref - schemeData.meanFlowSpeed ) / schemeData.meanFlowSpeed;
    else
      % Just return reference velocity at first position
      u_vec = ones(size(t_vec)) * ( schemeData.uref(1) - schemeData.meanFlowSpeed ) / schemeData.meanFlowSpeed;
    end
    
  end
  
else
  error('Chosen transientDef not defined!')
end


end

