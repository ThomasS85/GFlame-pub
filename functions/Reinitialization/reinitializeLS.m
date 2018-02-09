function [ y, schemeData ] = reinitializeLS( ~ , y, schemeData  )
% function data = reinitializeLS(grid, data0, accuracy, tMax, errorMax)
% reinitDemo: reinitializes the G-field.
%
% Inputs: 
%   - y             : G-field data (vector)
%   - schemeData    : Struct with all scheme data including a struct...
%       .solver           : Solver settings
%       .innerData{ii}    : Contains convectiveData (ii will be computed) which contains (among others)...
%         .velocity             : Function handle to velocity model function
%         .uref                 : Actual reference velocity at t (line vector or scalar)
%         .A_surf               : Actual flame surface area at t (line vector or scalar)
%         .t_vec                : Actual time t (line vector or scalar)
%         .grid                 : Grid for surface evaluation
%         .p                    : Struct with general settings
%       
%
% Outputs:
%   - y             : G-field data (vector)
%   - schemeData    : Modified version of input
%
%
% Converts an implicit surface function into a signed distance function
%   by iterative solution of the reinitialization equation.
% Contrary to signedDistanceIterative the subcell fix is used in order to
%   minimize the movement of the zero level iso line
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 08.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% Should data be reinitialized?
if ~mod(schemeData.counter-1 , schemeData.solver.reinit.intervalStep )
  
  % (1) Get user defined values 
  %   Accuracy
  accuracy = schemeData.solver.reinit.accuracy;
  %   Time at which to halt the reinitialization iteration
  tMax = schemeData.solver.reinit.tMax;
  %   When to  assume that reinitialization has converged?
  errorMax = schemeData.solver.reinit.errorMax;
  
  % (2) Find index of convectiveData in innerData struct
  tmp = cellfun(@(x) strcmp(func2str(x),'termConvection') , schemeData.innerFunc );
  indSD = find(tmp,1,'first');
  
  if isempty(indSD)
    error('No convective term specified or name of convective term unknown!')
  end
  
  
  % (3) Perform reinitialisation
  grid = schemeData.innerData{indSD}.grid;
  y = signedDistanceIterativeSubCellFix(grid, y, accuracy, tMax, errorMax);
  
  
end

end