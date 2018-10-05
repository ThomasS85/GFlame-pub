function [ res ] = SC_inverse_mapping_Residuum( xi , x , s )
%SC_INVMAPRES returns the residuum of the implicit inverse SC mapping 
%     function for a given xi and x
%
%  Inputs:  - xi       : Guess for solution as vector [real , imag]
%           - x        : Right hand side of equation as complex number
%
%  Outputs: - res      : Residuum of complex inverse mapping equation 
%                         SCmapping(xi)-x
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////

% Convert xi vector to xi as complex number:
xi_compl = xi(1) + 1i*xi(2);
% real part
res(1) =  real(s.SC_mapping(xi_compl) - x);
% Imaginary part
res(2) =  imag(s.SC_mapping(xi_compl) - x);

end

