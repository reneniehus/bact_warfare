% Invasion analysis
clc; clear all
% close all
tic;
%%%%%%%% get filename %%%%%%%%
p = mfilename('fullpath');
[pathstr, name, ext] = fileparts(p); % CHANGE IF YOU RUN IT IN ZOOLOGY
plotting_0 = 0; % no plotting for Zoology server
rpath = [pathstr '/Results/'];
fpath = [pathstr '/graphs/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
% Run description
% set initial condition ODE
Endtime = 48; % given in hours
TIME = [0,Endtime];
dt = 0.1;       % time step
% -----------
param.Ca0 = 0.1; param.Cb0 = param.Ca0; Ta0 = 0; Tb0 = 0;  % define initial conditions
param.N0 = 1;
param.KN = 5; % half-saturation constant for nutrient-dependent growth
param.mu = 10; % max growth rate
param.kay = 30; %0.7; % how many cells are killed per unit toxin
param.D = 0.20; % loss of toxin
% -------------
InitC = [param.Ca0;param.Cb0;Ta0;Tb0;param.N0];


%% Algorithm to get as fast as possible to exact f*
% starts at highest and lowest f
for fresini = [0 1]
    singular = 0; % boolean saying if singular strategy is found
    newres = 1; % is there a new resident strategy
    fres = fresini; % resident strategy
    forward = 1; % indicates is new mutant will be higher or lower than resident
    step = 0.1; % initial step for the mutant. Initially coarse
    nbflip = 0; % record the number of consecutive flips in direction
    minstep = 0.00001; % the precision we want for the strategy value
    %
    while singular ~= 1 % do this while you haven't found the singular strategy
        if newres == 1 % new resident strategy; calculate new residents avarage fitness
            %%%%%%%% local one-on-one competition (ODE solving)%%%%%%%%%%%%
            [t,X] = myODE_solver01([fres 0 0 0 0 0 0 0 0],[fres 0 0 0 0 0 0 0 0],param,InitC,[0,Endtime],dt); % solve ODE b vs b
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            wresav = X(1,end); % residents w as final biomass
        end
        fmut = fres + forward*step; % pick an invader that differs from resident according to direction
        if fmut < 0 % if mutant is smaller than 0, set to 0
            fmut = 0;
        end
        if fmut > 1
            fmut = 1; % if mutant is bigger than 1, set to 1
        end
        %%%%%%%% local one-on-one competition (ODE solving)%%%%%%%%%%%%
        [t,X] = myODE_solver01([fmut 0 0 0 0 0 0 0 0],[fres 0 0 0 0 0 0 0 0],param,InitC,[0,Endtime],dt); % solve ODE b vs b
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        wmut = X(1,end); % record mutant biomass at end of competition
        previousforward = forward; % saver: update previous direction
        previousfres = fres; % saver: update previous resident
        if wmut <= wresav % resident stays the resident
            forward = (-1) * forward; % change the dircetion for mutant
            newres = 0;
        end
        if wmut > wresav % the mutant invades and replaces the resident
            fres = fmut; % the new resident will take the mutant value
            disp(['Resident:' num2str(fres)])
            if fres == 0 || fres == 1
                singular = 1; % you have found a singular strategy
            end
            newres = 1; % we have a new resident
        end
        %
        if forward ~= previousforward % compare with previous direction
            nbflip = nbflip + 1;
        end
        if forward == previousforward % compare with previous direction
            nbflip = 0;
        end
        if previousfres == fres && nbflip > 1 % singular strategy localised
            if step <= minstep % precision high enough
                singular = 1; % we have found an optimal strategy
            end
            if step > minstep % will redo the whole procedure with tinier step
                step = step/10; % increase precision 10 times
                nbflip = 0; % reset the number of consecutive flips
            end
        end
    end % while you search for singular strategy
    fres % record the final fres (ESS)
end
%%
disp(['Done after ' secs2hms(toc)]);
