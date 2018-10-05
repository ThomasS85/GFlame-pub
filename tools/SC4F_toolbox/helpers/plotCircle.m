function plotCircle( x,y,r )
%PLOTCIRCLE Plots a circle around (x,y) with radius r

% Angular step
angStep = 1e-2;
ang=0:angStep:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'b'); hold on;
plot(x,y,'bx'); 

end

