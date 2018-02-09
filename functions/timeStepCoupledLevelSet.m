function [ t , y , yTax ] = timeStepCoupledLevelSet( taxMod , tSpan , y0 , solver , vel )
%TIMESTEPCOUPLEDLEVELSET Integrates coupled Tax-GFLAME model from starting
%from t0 to t0+tspan
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 01.12.2014 as part of GFLAME 0.1          //
% ////////////////////////////////////////////////////////

% Loop
while t<tend

% TaX: Acoustic
[yTax,taxMod.x,~] = sss.sim_RK4(taxMod.A,taxMod.B,taxMod.C,taxMod.D,taxMod.E,u,taxMod.x(:,end-1),...
  vel.Ts_Tax,vel.Ts_Tax,false);

% Update reference velocty according to TaX output
for ii=1:length(solver.schemeData.innerFunc)
  if strcmp( func2str(solver.schemeData.innerFunc{ii}), '')
    solver.schemeData.innerData{ii}.uref = urefTAX;
  end
end

% Level-Set solver
[ t , y , solver.schemeData ] = feval(solver.integratorFunc, solver.schemeFunc, tSpan, y0,...
  solver.integratorOptions, solver.schemeData);

% Evaluate surface area

end

end

