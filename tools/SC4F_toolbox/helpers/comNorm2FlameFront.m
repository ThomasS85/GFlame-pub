function [ myN ] = comNorm2FlameFront( flameFront_curve , x1F )
%COMNORM2FLAMEFRONT Computes normal sto flame front provided by 'flameFront_curve'

% Make sure xi_0 and x1F are line vectors
[z,s] = size(flameFront_curve);
if z>s; flameFront_curve = flameFront_curve.'; end;

% Get spatial spacing
dx1F = x1F(2) - x1F(1);

% Evaluate normal to flame front
[nDat] = compNormal2Curve_FPB( flameFront_curve );
nDat.x1F_nodes =  x1F(1:end-1) + dx1F/2;

% Get normals at x1F nodes via interpolation -> exclude boundaries!
myN = interp1( nDat.x1F_nodes , nDat.n_vec(1,:)+1i*nDat.n_vec(2,:) , x1F(2:end-1) );
% Add values at boundaries manually
myN = [ myN(1) , myN , myN(end) ];
% Normalize (probably not required, but to be sure...)
myN = myN ./ sqrt( real(myN).^2 + imag(myN).^2 );


%% Include BC
myN(1) = 1;
myN(end) = 1;

% Debug plot
% figure;hold on;axis equal;
% plot(flameFront_curve(1,:),flameFront_curve(2,:),'r')
% myScal = 2e-3;
% quiver(nDat.nodes(1,:),nDat.nodes(2,:),myScal*nDat.n_vec(1,:),myScal*nDat.n_vec(2,:),0)
% xlim([-10 10]*5e-3);ylim([0 1]*50e-3);

end

