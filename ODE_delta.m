function dXdt = ODE_delta(X,param)
Ca=X(1); Cb=X(2); Ta=X(3); Tb=X(4); N=X(5);
fa = X(6); fb = X(7);
%Cell Density A
% (N/(N+param.KN)) = N
dCadt =(1-fa)*param.mu*(N/(N+param.KN))*Ca - param.kay*Tb*Ca;
%Cell Density B
dCbdt =(1-fb)*param.mu*(N/(N+param.KN))*Cb - param.kay*Ta*Cb;
%Toxin A Concentration
dTadt = fa*(N/(N+param.KN))*Ca - (param.D)*Ta;
%Toxin B Concentration
dTbdt = fb*(N/(N+param.KN))*Cb - (param.D)*Tb;
% explicit Nutrients
dNdt = -1*(N/(N+param.KN))*(Ca + Cb) ;
%%%% put them all together
dXdt=[dCadt;dCbdt;dTadt;dTbdt;dNdt];
end
 
