function [ y , schemeData] = writeInOut2File( t , y , schemeData )
%WRITEINOUT2FILE Function evaluates writes input and output data of GLFLAME after an integrated time step.
% The data consists of u_ref(t) and A(t). Used as a post-time step function! 
% 
% Constructed as a post time step routine which is called by odeCFLn and
% therefore has the following in- and outputs:
%
% Inputs: 
%   - t             : Actual time
%   - y             : G-field data (vector)
%   - schemeData    : Struct with all scheme data including a struct...
%       .writeInOut       : Which contains...
%         .fileA                : File handle for flame surface output
%         .fileV                : File handle for reference velocity
%       .solver           : Solver settings
%         .intervalWriteOI      : Defines output intervall of u_ref(t)/ A(t)
%       .innerData{ii}    : Contains convectiveData (ii will be computed) which contains (among others)...
%         .velocity             : Function handle to velocity model function
%         .uref                 : Actual reference velocity at t (line vector or scalar)
%         .A_surf               : Actual flame surface area at t (line vector or scalar)
%         .t_vec                : Actual time t (line vector or scalar)
%         .grid                 : Grid for surface evaluation
%         .p                    : Struct with general settings
%       .counter          : Counts time steps performed
%       
%
% Outputs:
%   - y             : G-field data (vector)
%   - schemeData    : Modified version of input
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.de).        //
% // Created, 23.03.2016 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////


% Function internal settings
% Should reinitialised field be passed back to callinf function? 
%   Disadvantage: Since extrapolation BC are used this could lead to an internal field which is not compatible
%                 with BC. This might lead to noise (solution: change BC in delta function; This, however, 
%                 might effect the quality of the signed distance function close to the boundaries)
%   Advantage:    Since reinitialisation is computational expensice, the reinitialized field could also be
%                 passd to the calling function and is thus be used for the solution. However, the control of
%                 the reinitialisation would become more complex.
doReinit = 0;  % Recommended: 0

% Should data be evaluated and written to file?
if ~mod(schemeData.counter-1 , schemeData.solver.intervalWriteOI )
  % Yes, evaluate data and write to file
  
  
  % (1) Find index of convectiveData in innerData struct
  tmp = cellfun(@(x) strcmp(func2str(x),'termConvection') , schemeData.innerFunc );
  indSD = find(tmp,1,'first');
  
  if isempty(indSD)
    error('No convective term specified or name of convective term unknown!')
  end
  
  
  % (2) Get back the correctly shaped data array
  data = reshape(y, schemeData.innerData{indSD}.grid.shape);

  
  % (3) Compute velocity field and write reference velocity
  [ vel , convectiveData , uref ] =...
    feval( schemeData.innerData{indSD}.velocity , t , data, schemeData.innerData{indSD} );
  schemeData.innerData{indSD} = convectiveData;
  schemeData.innerData{indSD}.velocityField = vel;
  
  % Compute reference velocity and save to first position in uref
%   uref = schemeData.innerData{indSD}.velocityField{ref.component}(ref.pos(1),ref.pos(2));
   
  % (4) Compute and write flame surface area (data is also reinitialized around 0-iso line!)
  if doReinit
    [ A_surf , data] = calcSurfaceArea( data ,schemeData.innerData{indSD}.grid , schemeData ,...
      schemeData.innerData{indSD}.p );
    % Transform data back to vector (which has been reinitialized)
      y = data(:);
      
  else
    [ A_surf ] = calcSurfaceArea( data ,schemeData.innerData{indSD}.grid , schemeData ,...
      schemeData.innerData{indSD}.p );
    
  end

  % (5) Write time to time vector
  schemeData.innerData{indSD}.uref = [ uref ; schemeData.innerData{indSD}.uref(1:end-1) ];
  schemeData.innerData{indSD}.A_surf = [ A_surf ; schemeData.innerData{indSD}.A_surf(1:end-1) ];
  schemeData.innerData{indSD}.t_vec = [ t ; schemeData.innerData{indSD}.t_vec(1:end-1) ];
  
  % (6) Write to file
  fprintf( schemeData.writeInOut.fileV , '\n%.8f , %.8f' , t , full(schemeData.innerData{indSD}.uref(1)) );
  fprintf( schemeData.writeInOut.fileA , '\n%.8f , %.8f' , t , full(schemeData.innerData{indSD}.A_surf(1)) );
  
end

% Raise counter by 1
schemeData.counter = schemeData.counter + 1;


end

