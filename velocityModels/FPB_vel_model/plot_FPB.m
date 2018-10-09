function [vel] = plot_FPB(FPB,p)
% plot_FPB plots the requested velocity field due to the mean flow, the
% acoustics and the gas expansion
%
% Inputs:
%
%   - data        - G-Field matrix
%
%   - schemeData  - Struct which contains information about all nescessary
%                   parameters
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%  by Axel Zimmermann (09.2018)  %
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% velocity field due to ACOUSTICS
%should acoustic velocity field be plotted?
if strcmp(FPB.plot.plot_acoustics,'y')
  u_ac= FPB.plot.ac;
else
  u_ac= 0;
end

%% should mean flow be plotted?
if strcmp(FPB.plot.plot_mean_flow,'y')
  u_mean=FPB.plot.meanFlowSpeed;
else
  u_mean=0;
end

FPB.velPar.sourceDat.G(end) = (u_mean+u_ac)*FPB.plot.R_i;

%% velocity field due to Flame Front
%should source velocity field be plotted ?

if ~strcmp(FPB.plot.plot_source_field,'y')
  FPB.velPar.sourceDat.G(1:end-1)=0;
end

if strcmp(p.CombType,'backwardFacingStep')
  [ u_ges ,~ , ~ ] = evalFlowField_physicalDomain(...
    FPB.plot.grid.xi , FPB.velPar , FPB.plot.myMap,'xi','L1','doKutta',FPB.velPar.Kutta.position_Kutta_xi,'noPanel2Point');
elseif strcmp(p.CombType,'duct')
  [ u_ges ,~ , ~ ] = evalFlowField_physicalDomain(...
    FPB.plot.grid.xi , FPB.velPar , FPB.plot.myMap,'xi','L1','noKutta','noPanel2Point');
  u_ges=conj(u_ges);
end

%prepare velocity field for output
vel{1}=real(u_ges);
vel{2}=imag(u_ges);


end

