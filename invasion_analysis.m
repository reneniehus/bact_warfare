% from Battle_bychone17_PIP_testpars


% Optimising an army of strategies
% Invasion analysis of the ODE (Niehus) Toxin ODE
clc; clear all
% close all
tic;
%%%%%%%% get filename %%%%%%%%
p = mfilename('fullpath');
[pathstr, name, ext] = fileparts(p); % CHANGE IF YOU RUN IT IN ZOOLOGY
%pathstr = '/home/rene/Dropbox/ZoologyProject/ToxinRegulation2015';%carful, this is for SSH running
%name = 'Battle_bychone10_ODE_eitherSignal_A_02space';
plotting_0 = 0; % no plotting for Zoology server
rpath = [pathstr '/Results/'];
fpath = [pathstr '/graphs/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endtimewarning = 0;
% Settings
max_end_time = 48;
% Run description
%% set initial condition ODE and PDE
Endtime = max_end_time; % given in hours
TIME = [0,Endtime];
dt = 0.1;       % time step
% -----------
param.Ca0 = 0.005; param.Cb0 = param.Ca0; Ta0 = 0; Tb0 = 0;  % define initial conditions
    param.gamma = 1; % the consumption of nutrients
    param.N0 = 1;
    param.KN = 5; % half-saturation constant for nutrient-dependent growth
    param.d = 0.2*1000; % distance between colonies and outside nutrient source (micrometer)
    param.cd = 0.2*1000; % width of colonies (micrometer)
    param.bd = 0*1000; % distance between the colonies (micromter)
    param.fraction_with_cells = (2*param.cd)/(2*param.cd+2*param.d+param.bd);
    %% set paramters
    param.mu = 10; % max growth rate
    param.q = 1; % toxin production parameter
    param.kay = 40; %0.7; % how many cells are killed per unit toxin
    param.beta = 0; % food value of killing a cell
    param.l = 0; % loss of toxin
    param.D = 0.20; % 0.04; % diffusion loss of toxin in ODE
    param.Diff_T = 15*1000; % Toxin diffusion (time unit: hours, space unit: micrometers)
    param.Diff_N = 10000000;  % Nutrient diffusion (time unit: hours, space unit: micrometers)
% -------------
InitC = [param.Ca0;param.Cb0;Ta0;Tb0;param.N0];
% nutrient amount available per all cells at any time; instant diffusion model
%%%%%%%%%%%%% Start different functions here


if 1 == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run f* algorithm
    %%
    %% Algorithm to get as fast as possible to exact f*
    % Settings
    space_on = 0;
    plotting = 0;
%     res_v = 0.2:0.2:0.8;
%     % get to initial resident strategy
%     resws = zeros(1,numel(res_v));
%     for res_str_c = 1 : numel(res_v) % find the highest average resident biomass
%         fb = res_v(res_str_c); % fb is resident strategy
%         % get to resident on-its-own fitness = final biomass
%         if space_on == 1
%             [t,X] = my_own_pde41K_sensALL([fb 0 0 0 0 0 0],[fb 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
%         else [t,X] = my_own_ode46K_sensALL([fb 0 0 0 0 0 0],[fb 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
%         end
%         resws(res_str_c) = X(1,end);
%     end
%     fresini = res_v(resws == max(resws)); % initial resident strategy
%     fresini = fresini(1); % in case there a multiple highest biomass residents
for fresini = [0 1]
    singular = 0; % boolean saying if singular strategy is found
    newres = 1; % is there a new resident strategy
    fres = fresini; % resident strategy
    forward = 1; % indicates is new mutant will be higher or lower than resident
    step = 0.1; % initial step for the mutant. Initially coarse
    nbflip = 0; % record the number of consecutive flips in direction
    minstep = 0.00001; % the precision we want for the strategy value
    while singular ~= 1 % do this while you haven't found the singular strategy
        if newres == 1 % new resident strategy; calculate new residents avarage fitness
            if space_on == 1
                [t,X] = my_own_pde41K_sensALL([fres 0 0 0 0 0 0],[fres 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
            else [t,X] = my_own_ode46K_sensALL([fres 0 0 0 0 0 0],[fres 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
            end
            wresav = X(1,end); % residents w as final biomass
        end
        fmut = fres + forward*step; % pick an invader that differs from resident according to direction
        if fmut < 0 % if mutant is smaller than 0, set to 0
            fmut = 0;
        end
        if fmut > 1
            fmut = 1; % if mutant is bigger than 1, set to 1
        end
        if space_on == 1
            [t,X] = my_own_pde41K_sensALL([fmut 0 0 0 0 0 0],[fres 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
        else [t,X] = my_own_ode46K_sensALL([fmut 0 0 0 0 0 0],[fres 0 0 0 0 0 0],param,InitC,[0,Endtime],dt,plotting); % solve ODE b vs b
        end
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
        % wmut = wres only when: both die or no one produces toxins
        % when wmut = wres
        %         if wmut == wresav %
        %             disp('wmut = wres ! Do something!')
        %         end
        
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
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run f* algorithm



disp(['Done after ' secs2hms(toc)]);
