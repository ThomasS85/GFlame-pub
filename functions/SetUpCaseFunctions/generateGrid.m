function [ myGrid ] = generateGrid( solver )
%GENERATEGRID generates the Grid for the G-equation

% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

myGrid.dim = solver.dim;
myGrid.min = [ solver.x1Lim(1) ; solver.x2Lim(1) ];
myGrid.dx = solver.dxi;
% grid.max = [ solver.x1Lim(2) ; solver.x2Lim(2) ] - grid.dx;
myGrid.max = [ solver.x1Lim(2) ; solver.x2Lim(2) ] ;

myGrid = processGrid(myGrid);

end

