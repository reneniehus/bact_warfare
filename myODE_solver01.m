function [t, X] = myODE_solver01(a_strat,b_strat,param,InitC,TIME,dt)
%my_own_ode45 Solving my ODEs with fixed timestep
StartTime = TIME(1); Endtime = TIME(2);
t = StartTime:dt:round(Endtime,1);     % time discretization
N = length(t);  % number of time points
X = zeros(7,N+10); % initialise X
% Initialise X
X(1:5,1) = InitC; % Ca, Cb, Ta, Tb, N
[fa,fb] = return_f(a_strat,b_strat,X(1:5,1));
X(6:7,1) = [fa;fb];
i = 1; % initiate step counter
dSV = inf; % initiate change condition
while t(i) < round(Endtime,1)
    i = i + 1;
    dSV = ODE_delta(X(:,i-1),param)*dt; % delta state variables
    X(1:5,i) = X(1:5,i-1) + dSV;
    [fa,fb] = return_f(a_strat,b_strat,X(1:5,i));
    %
    X_cells = X(1:2,i) ;
    X_cells( X_cells(:,1) < 1/1000000 ) = 0 ; % correct cells to that they go extinct
    X(1:2,i) = X_cells ;
    %
    X((X(:,i) < 1/100000000),i) = 0; % so that state variables are never < 0
    X(6:7,i) = [fa;fb];
end
t = t(:,1:i); % cut away empty time points
X = X(:,1:i); % cut away empty time points
end % end of function


