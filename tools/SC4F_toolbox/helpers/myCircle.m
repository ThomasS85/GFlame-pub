function [ xOut , yOut ] = myCircle( x0 , y0 , myRad , numEls )
%MYCIRCLE Returns x and y coordinates of a circle with radius myRad around the points x0 and y0. numEls points
%  are created

myTheta = linspace(0,2*pi,numEls);
xOut = myRad * cos(myTheta) + x0;
yOut = myRad * sin(myTheta) + y0;

end

