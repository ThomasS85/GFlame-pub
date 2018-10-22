function [ myLine ] = checkNcorrectLine_geometry( myLine , p )
%CHECKNCORRECTLINE_GEOMETRY Checks if a given line is conform with the combustor geoemtry, i.e. has no points
%outside the flow domain
%
% Inputs:
%   myLine  - line in L1 coordinates (complex vector)
%   p       - flame definitions (struct)
%

if strcmpi(p.CombType,'backwardFacingStep')
 % Find points that are outside of the domain
  pointsOut1 = imag(myLine)>p.R_i & real(myLine)<0;
  pointsOut2 = imag(myLine)<0;
  pointsOut3 = imag(myLine)>p.R_a & real(myLine)>0;
  
  % remove
  myLine = myLine(~pointsOut1 | ~pointsOut2 | ~pointsOut3);
  
elseif strcmpi(p.CombType,'duct')
  % Find points that are outside of the domain
  pointsOut1 = imag(myLine)>p.R_i;
  pointsOut2 = imag(myLine)<0;
  
  % remove
  myLine = myLine(~pointsOut1 | ~pointsOut2);
  
else
  error('Unknown combustor geometry!')
end

end

