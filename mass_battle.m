% Evolutionary tournament competing all 3 sensing strategies
clc; clear all
tic; % for timer
 
%%%%%%%% get filename %%%%%%%%
p = mfilename('fullpath');
[pathstr, name, ext] = fileparts(p); 
name = 'allevolvePNAS';
rpath = [pathstr '/Results/'];
fpath = [pathstr '/graphs/'];
name_sens = 'THREE' ; % as addition to file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% which sensor types will be in the tournament
random_guy_sensors = [1 1 0 1]; % [Nu TB TA QS]
% Nu: Nutrient sensing
% TB: Toxin sensing
% TA: own toxin sensing (not included)
% QS: Quorum sensing
 
% independent repeats of turnament (allows running more)
alg_num1 = 1;
alg_num2 = 5;
% number of evolution rounds (allows extending runs)
startround = 1;
endround = 50;
% tournament rules (control parameters)
n_rounds_migration = 20 ;
%
n = 60; %60 number of strategies in the pool
n_fights = 60/n; %16/n; % fraction of population that cells have to meet
top_f = 4/n; % 4/n fraction of best strategies that go straight into next gen
% mutation = 6/n
frac_mig = 10/n; % 10/n fraction of entirly new guesses
frac_mut = 46/n ; % 60-(4+10)
% ecological parameters
param.Ca0 = 0.1; param.Cb0 = param.Ca0; Ta0 = 0; Tb0 = 0;  % initial conditions
Endtime = 48; % duration of single pair-wise competition (in hours)
param.N0 = 1;
param.KN = 5; % half-saturation constant for nutrient-dependent growth
param.mu = 10; % max growth rate
param.kay = 20; %0.7; % how many cells are killed per unit toxin
param.D = 0.10; % loss of toxin
 
% for ODE solver
dt = 0.1;       % time stepping
InitC = [param.Ca0;param.Cb0;Ta0;Tb0;param.N0]; % initial conditions
 
% ranges [fnull | fN  UN | fTB  UTB | fTA  UTA | fQS QS]
% initiate
param_range = zeros(9,3); % give ranges and mutation steps
% Baseline investment range
my_max_tox_guess = 0.01; % my guess for highets tox concentration
% changed limits
param_range(1,:) = [0   1/50 1]; % fNull_range
% Sensing nutrients
param_range(2,:) = [-1   1/50 1]; % fNut_range
param_range(3,:) = [0.85   1/50    1];
% Sensing opposing toxins
param_range(4,:) = [-1   1/100 1]; % fTopp_range
param_range(5,:) = [0.00   my_max_tox_guess/50    my_max_tox_guess];
% Sensing own toxins
param_range(6,:) = [-1   1/50 1]; % fown_range
param_range(7,:) = [0.00   my_max_tox_guess/50    my_max_tox_guess];
% Sensing QS
param_range(8,:) = [-1   1/50 1]; % fOQ_range
param_range(9,:) = [0.00   1/50 1.2]; % fOQ_range
 
% useful things
numberofalgs = alg_num2 - alg_num1 + 1;
numberofrounds = endround - startround + 1; % number of rounds
for alg_runs = alg_num1 : alg_num2
    % load previous results
    loadname = [rpath name name_sens num2str(alg_runs) 't' num2str(startround - 1) '.mat'];
    if exist(loadname, 'file') % File exists.
        load(loadname); % loads: save(savename,'WARR_save','end_run','n','frac_mig','WARR');
        WARR_save = zeros(n,9+2+3,endround);
        WARR_save(:,:,1:(startround - 1)) = WARR_saveS;
        WARR = WARRS;
        disp('Loaded previous file!')
    else % File does not exist.
        if startround ~= 1
            disp('Previous file does not exist')
        end
        % 1) Generate initial random population of N strategies
        WARR = zeros(n,9 + 2); % initiate warrior matrix with tagPool tagTop
        WARR_save = zeros(n,9 + 2 + 3,numberofrounds); % save the population+fitness+toxs each round
        for seed1_WARR = 1 : n % seed the initial n warriors
            % first select sensing for this guy
            random_guy_sensors_i = [0 0 0 0];
            random_guy_sensors_i(randsample(find(random_guy_sensors ==1),1)) = 1 ;
            % then create the random guy until it works
            mut_strat_test = create_random_strategy(param_range,random_guy_sensors_i) ;
            light = test_strategy_range_qs(mut_strat_test,param_range);
            while light == 0 % while strategy gets red light try again
                mut_strat_test = create_random_strategy(param_range,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
                light = test_strategy_range_qs(mut_strat_test,param_range);
            end
            WARR(seed1_WARR,1:9) = mut_strat_test; WARR(seed1_WARR,10:11) = [0 0]; % new tags TagPool TagTop
        end % seed the initial n warriors
    end
    % MAIN LOOP
    FITS = zeros(n,1);
    TOXS = zeros(n,2); % for each strain context dep tox average + tox peak
    WARR_next = zeros(n,9 + 2); % initiate next-gen warrior matrix + 2 tags
    loopStart = tic; % TIC, pair loopStart
    for round_c = startround : endround % through rounds of algorithm
        disp(['Run:' num2str(alg_runs) '/' num2str(alg_num2) ' round ' num2str(round_c) '/' num2str(endround)])
        % 2) EVALUATE ALL STRATEGIES in k random battles (might fight yourself)
        save_biom = NaN(n,n); % matrix to save biomasses
        for through_warriors = 1 : n % give each warrior its fitness
            FITS(through_warriors) = 0; % set warrior's fitness to zero
            TOXS(through_warriors,:) = [0,0]; % can help to troubleshoot
            % get the opponents
            compets = randsample(n,(n*n_fights)); % vector of competitors, may include focal guy itself
            for battles = 1 : numel(compets) % through battles
                if ~isnan(save_biom(through_warriors,compets(battles))) % check if there is already a result saved
                    wmut =  save_biom(through_warriors,compets(battles));
                else
                    %%%run the local competition (solving ODE numerically)
                    [t,X] = myODE_solver01(WARR(through_warriors,:),WARR(compets(battles),:),param,InitC,[0,Endtime],dt);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    wmut = X(1,end); % record a's end-biomass
                    % save the opponents biomass too
                    save_biom(compets(battles),through_warriors) = X(2,end); % record b's end-biomass
                end
                FITS(through_warriors) = FITS(through_warriors) + wmut; % fitness = sum of fight biomasses
                % save the opponents biomass too
            end % through battles
        end % give each warrior its fitness
        WARR_save(:,:,round_c) = [WARR,TOXS,FITS]; % save [Top|Children|Mutations|Migrations] apart from first gen
        if round_c ~= numberofrounds % 3) !!!!NEW GENERATION!!!!
            % 3 A) keep the absolute top fraction top_f as they are
            fitscopy = FITS;
            for top_place = 1 : (top_f*n) % fill top_f*n top places
                [a, Index_max] = max(fitscopy);
                WARR_next(top_place,:) = WARR(Index_max,:);
                WARR_next(top_place,10) = WARR_next(top_place,10) + 1; % tagPool +1
                WARR_next(top_place,11) = WARR_next(top_place,11) + 1; % tagTop +1
                fitscopy(Index_max) = -inf;
            end %  fill top_f*n top places
            for mutant_place = (top_f*n + 1) : (n - n*frac_mig) % make n*(1 -top_f  - frac_mig) mutants
                % select the mother for mutatant children depending on fit
                mother_strat = WARR(randsample(n,1,true,FITS),:);
                which_sensor = (mother_strat([2,4,6,8]) ~= 0);
                which_gene_mut = randi(3); % which part of the sens stragegy mutates
                mut_strat_test = strat_mutate(which_sensor,which_gene_mut,mother_strat); % mutation happens
                light = test_strategy_range_qs(mut_strat_test,param_range);
                while light == 0 % redo mutation until it works
                    mut_strat_test = strat_mutate(which_sensor,which_gene_mut,mother_strat); % mutation happens
                    light = test_strategy_range_qs(mut_strat_test,param_range);
                end
                WARR_next(mutant_place,:) = mut_strat_test;
                WARR_next(mutant_place,10) = WARR_next(mutant_place,10) + 1; % tagPool +1
            end % make n*(1 -top_f - frac_mig) mutants
            
            % in first 20 generations only migration
            if round_c <= n_rounds_migration
                migrant_place_first = 1 ;
            else
                migrant_place_first = (n - n*frac_mig + 1) ;
            end
            % after that selection
            
            for migrant_place = migrant_place_first : n % make n*frac_mig migrants
                % first select sensing for this guy
                random_guy_sensors_i = [0 0 0 0];
                if sum(random_guy_sensors) >= 2
                    random_guy_sensors_i(randsample(find(random_guy_sensors ==1),1)) = 1 ;
                else
                    random_guy_sensors_i = random_guy_sensors;
                end
                % then create the random guy until it works
                mut_strat_test = create_random_strategy(param_range,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
                light = test_strategy_range_qs(mut_strat_test,param_range);
                while light == 0
                    mut_strat_test = create_random_strategy(param_range,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
                    light = test_strategy_range_qs(mut_strat_test,param_range);
                end
                WARR_next(migrant_place,1:9) = mut_strat_test;
                WARR_next(migrant_place,9+1) = 0; % tagPool reset
                WARR_next(migrant_place,9+2) = 0; % tagTop reset
            end % make n*frac_mig migrants
            WARR = WARR_next; % copy new generation and overwrite the old
        end % create new generation
        %% Time prediction unit
        if round_c == (startround + 1) && alg_runs == alg_num1
            looptimeXtimes = toc(loopStart);
            time_in_sec = looptimeXtimes/2 * ((numberofrounds - 2) + (numberofalgs-1)*numberofrounds); % expected time in seconds
            time_in_hms = secs2hms(time_in_sec); % expected time in hours,minutes,secs
            display1 = [datestr(now) '; Time left: ' time_in_hms]; % print expected time
            disp(display1);
        end
    end % through rounds of algorithm
    % save WARR_save forever
    savename = [rpath name name_sens num2str(alg_runs) 't' num2str(endround) '.mat']; % zero small much
    WARR_saveS = WARR_save;
    endroundS = endround;
    nS = n;
    frac_migS = frac_mig; n_fightsS = n_fights; top_fS = top_f;
    WARRS = WARR; 
    save(savename,'WARR_saveS','endroundS','nS','frac_migS','WARRS',...
        'n_fightsS','top_fS');
    % save again for R
    savename2 = [rpath name name_sens num2str(alg_runs) 't' num2str(endround) '_onlyWARR.mat']; % zero small much
    save(savename2,'WARR_saveS');
    
end
%%%%%%%%%%%%%% the end
