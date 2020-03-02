% Optimising an army of all 3 strategies strategies
% Invasion analysis of the ODE (Niehus) Toxin ODE
clc; clear all
% close all
tic;
%%%%%%%% get filename %%%%%%%%
p = mfilename('fullpath');
[pathstr, name, ext] = fileparts(p); % CHANGE IF YOU RUN IT IN ZOOLOGY
name = 'allevolvePNAS';
%pathstr = '/home/rene/Dropbox/ZoologyProject/ToxinRegulation2015';%carful, this is for SSH running
%name = 'Battle_bychone10_ODE_eitherSignal_A_02space';
plotting_0 = 1; % no plotting for Zoology server
rpath = [pathstr '/Results/'];
fpath = [pathstr '/graphs/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 5
    if jj==1
        random_guy_sensors = [1 0 0 0]; % [Nu TB TA QS]
    elseif jj==2
        random_guy_sensors = [0 1 0 0]; % [Nu TB TA QS]
    elseif jj==3
        random_guy_sensors = [0 0 0 1]; % [Nu TB TA QS]
    elseif jj==4
        random_guy_sensors = [0 0 0 0]; % [Nu TB TA QS]
    elseif jj==5
        random_guy_sensors = [1 1 0 1]; % [Nu TB TA QS]
    end
    name_sens_v = ["Nu","TB","TA","QS"] ;
    name_sens = char(name_sens_v(random_guy_sensors==1)) ; % give the name of the sensing
    if isempty(name_sens)
        name_sens = 'Z' ;
    end
    %
    name_sens = 'THREE' ;
    
    
    one_over_landing_rate = 0 ; % 0  1/10 1 ; % determines how long the first strains waits for second
    if one_over_landing_rate==0
        landing_tag='zero';
    elseif one_over_landing_rate==1/10
        landing_tag='short';
    elseif one_over_landing_rate==1
        landing_tag='long';
    end
    alg_num1 = 1;
    alg_num2 = 5;
    % evolution rounds
    startround = 1;
    endround = 50; % 100
    
    % Settings
    max_end_time = 1000;
    endtimewarning = 0;
    plotting = 0;
    % Run description
    %% set initial condition ODE and PDE
    Endtime = max_end_time; % given in hours
    TIME = [0,Endtime];
    dt = 0.1;       % time step
    param.Ca0 = 0.1; param.Cb0 = param.Ca0; Ta0 = 0; Tb0 = 0;  % define initial conditions
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
    param.kay = 30; %0.7; % how many cells are killed per unit toxin
    param.beta = 0; % food value of killing a cell
    param.l = 0; % loss of toxin
    param.D = 0.10; % 0.04; % diffusion loss of toxin in ODE
    param.Diff_T = 15*1000; % Toxin diffusion (time unit: hours, space unit: micrometers)
    param.Diff_N = 10000000;  % Nutrient diffusion (time unit: hours, space unit: micrometers)
    InitC = [param.Ca0;param.Cb0;Ta0;Tb0;param.N0];
    if 1 == 1 %% Optimise an army of strategies with only single sensors
        %%
        % ranges [fnull | fN  UN | fTB  UTB | fTA  UTA | fQS QS]
        space_on = 0; % space on or not
        if space_on == 1 % in space
            my_max_tox_guess = 0.2; % my guess for highets tox concentration
            my_max_nut_guess = 10/800;
        else % in well-mixed
            my_max_tox_guess = 0.01; % my guess for highets tox concentration
            my_max_nut_guess = 1;
        end
        give_tox_war = 0; % give toxin warning
        % present
        param_range = zeros(9,3); % give ranges and mutation steps
        % Baseline investment range
        % CAREFUL changed limits
        param_range(1,:) = [0   1/50 1]; % fNull_range
        % Sensing nutrients
        param_range(2,:) = [-1   1/50 1]; % fNut_range
        param_range(3,:) = [0.85   my_max_nut_guess/50    my_max_nut_guess]; % Nut threshold % keep the ranges now as they are.
        % Strategies should come in from a vast range of ranges, that is
        % fine. The exact range is up to me, and will of course make a
        % difference. But I decide to leave them like this.
        % Sensing opposing toxins
        param_range(4,:) = [-1   1/100 1]; % fTopp_range
        param_range(5,:) = [0.00   my_max_tox_guess/50    my_max_tox_guess]; % T range
        % Sensing own toxins
        param_range(6,:) = [-1   1/50 1]; % fown_range
        param_range(7,:) = [0.00   my_max_tox_guess/50    my_max_tox_guess]; % T range
        % Sensing QS
        param_range(8,:) = [-1   1/50 1]; % fOQ_range
        param_range(9,:) = [0.00   1/50 1.2]; % fOQ_range
        % useful things
        numberofalgs = alg_num2 - alg_num1 + 1;
        numberofrounds = endround - startround + 1; % number of rounds
        for alg_runs = alg_num1 : alg_num2
            % useful stuff for the the algorithm
            n = 60; %60 number of strategies in the pool
            n_fights = 60/n; %16/n; % fraction of population that cells have to meet
            top_f = 4/n; % 4/n fraction of best strategies that go straight into next gen
            child_f = 0/n; % HAS TO BE EVEN NUMBER fraction of cross-over children
            % mutation = 6/n
            frac_mig = 10/n; % 10/n fraction of entirly new guesses
            frac_mut = 46/n ; % 60-(4+10)
            % other tuning parameters
            f_all = 0; % fraction of new guesses that can sense all
            std_tiny = 0.005; % standard diviation for tiny steps
            std_small = 0.05; % standard diviation for small steps
            opp_tox_th_rang = 2;
            % load previous results
            loadname = [rpath name name_sens num2str(alg_runs) 't' num2str(startround - 1) landing_tag '.mat'];
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
                % completely random guess of toxin max concentration
                for seed1_WARR = 1 : n % seed the initial n warriors
                    % first select sensing for this guy
                    random_guy_sensors_i = [0 0 0 0];
                    if sum(random_guy_sensors) >= 2
                        random_guy_sensors_i(randsample(find(random_guy_sensors ==1),1)) = 1 ;
                    else
                        random_guy_sensors_i = random_guy_sensors;
                    end
                    % then create the random guy until it works
                    mut_strat_test = create_random_guy_range_qs(param_range,f_all,random_guy_sensors_i) ;
                    light = test_strategy_range_qs(mut_strat_test,param_range);
                    while light == 0 % while strategy gets red light try again
                        mut_strat_test = create_random_guy_range_qs(param_range,f_all,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
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
                save_Tav = NaN(n,n);
                save_Tmax = NaN(n,n);
                for through_warriors = 1 : n % give each warrior its fitness
                    FITS(through_warriors) = 0; % set warrior's fitness to zero
                    TOXS(through_warriors,:) = [0,0]; % set toxin measure to zero
                    % get the opponents
                    compets = randsample(n,(n*n_fights)); % vector of competitors, may include focal guy itself
                    for battles = 1 : numel(compets) % through battles
                        if ~isnan(save_biom(through_warriors,compets(battles))) % check if there is already a result saved
                            wmut =  save_biom(through_warriors,compets(battles));
                            avtox = save_Tav(through_warriors,compets(battles));
                            max_toxAB = save_Tmax(through_warriors,compets(battles));
                        else
                            if space_on == 1
                                [t,X] = my_own_pde40_sensALL(WARR(through_warriors,:),WARR(compets(battles),:),param,InitC,[0,Endtime],dt); % solve ODE b vs b
                            else
                                % decide which strain comes first and which
                                % 2nd
                                if rand < 0.5
                                    str_ord_delay = [1,2]; % strain arrival order
                                else
                                    str_ord_delay = [2,1]; % strain arrival order
                                end
                                InitC_delay = InitC; % update intiation conditions accordingly
                                InitC_delay(str_ord_delay(2)) = 0 ; % set initial biomass to zero
                                Endtime_delay = exprnd(one_over_landing_rate) ; % sample the delay distribution, mostly near zero
                                [t,X] = my_own_ode48QS_sensALL_noTD(WARR(through_warriors,:),WARR(compets(battles),:),param,InitC_delay,[0,Endtime_delay],dt,endtimewarning,plotting); % solve ODE b vs b
                                int_result_delay = X(:,end) ; % get the inbetween results
                                InitC_delay = int_result_delay(1:5) ; % transfer first species over, it's toxin and nutrients
                                InitC_delay(str_ord_delay(2)) = InitC(str_ord_delay(2)) ; % then put the second strain to it, of course it has no toxins
                                %
                                Endtime_delay = Endtime - Endtime_delay ; % do the rest of the competition
                                [t,X] = my_own_ode48QS_sensALL_noTD(WARR(through_warriors,:),WARR(compets(battles),:),param,InitC_delay,[0,Endtime_delay],dt,endtimewarning,plotting);
                            end
                            
                            wmut = X(1,end); % record a's end-biomass
                            avtox = sum(X(3,:),2); % sum/av toxin strain A
                            max_toxAB = max(X([3,4],:),[],2);
                            vmtoxconc = max(max_toxAB); % very maximum toxin concentration
                            if vmtoxconc > my_max_tox_guess && give_tox_war == 1
                                disp(['Increase max tox over ' num2str(vmtoxconc)])
                            end
                            save_biom(compets(battles),through_warriors) = X(2,end); % record b's end-biomass
                            save_Tav(compets(battles),through_warriors) = sum(X(4,:),2); % sum/av toxin strain B
                            save_Tmax(compets(battles),through_warriors) = max_toxAB(2); % tox max for strain B
                        end
                        FITS(through_warriors) = FITS(through_warriors) + wmut; % fitness = sum of fight biomasses
                        TOXS(through_warriors,1) = TOXS(through_warriors,1) + dt*avtox/(n*n_fights*Endtime); % add to av toxin
                        if max_toxAB(1) > TOXS(through_warriors,2)
                            TOXS(through_warriors,2) = max_toxAB(1);
                        end
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
                    for children_places = 1 : (child_f*n/2) % 3 c) make child_f*n cross-over children
                        the_parents_v = randsample(n,2,true,FITS); % select 2 parents based on fitness, double
                        % [fnull | fN UN | fTB UTB | fTA UTA]
                        %cross_point = randi(3); % select on of the 3 possible crossover-points
                        parent_1 = WARR(the_parents_v(1),:); % parent1
                        parent_2 = WARR(the_parents_v(2),:); % parent2
                        % get the parents sensors
                        %parent_1_sens_num_all = find(parent_1(2:end) ~= 0);
                        %parent_2_sens_num_all = find(parent_2(2:end) ~= 0);
                        child_1 = [parent_1(1) parent_2(2:end)];
                        child_2 = [parent_2(1) parent_1(2:end)];
                        child_tags = [child_1(9 + 1:9 + 2); child_2(9 + 1 : 9 + 2)]; % get both child tags
                        findrow = find(max(child_tags(:,1)) == child_tags(:,1));
                        if numel(findrow) == 1 % if one parent was older
                            child_1(9 + 1 : 9 + 2) = child_tags(findrow,:);
                            child_2(9 + 1 : 9 + 2) = child_tags(findrow,:);
                        end
                        light1 = test_strategy_range_qs(child_1,param_range);
                        light2 = test_strategy_range_qs(child_2,param_range);
                        while light1 == 0 && light2 == 0 % as long as both children die, try again
                            the_parents_v = randsample(n,2,true,FITS); % select 2 parents based on fitness, double
                            % [fnull | fN UN | fTB UTB | fTA UTA]
                            %cross_point = randi(3); % select on of the 3 possible crossover-points
                            parent_1 = WARR(the_parents_v(1),:); % parent1
                            parent_2 = WARR(the_parents_v(2),:); % parent2
                            % get the parents sensors
                            %parent_1_sens_num_all = find(parent_1(2:end) ~= 0);
                            %parent_2_sens_num_all = find(parent_2(2:end) ~= 0);
                            child_1 = [parent_1(1) parent_2(2:end)];
                            child_2 = [parent_2(1) parent_1(2:end)];
                            child_tags = [child_1(9 + 1:9 + 2); child_2(9 + 1 : 9 + 2)]; % get both child tags
                            findrow = find(max(child_tags(:,1)) == child_tags(:,1));
                            if numel(findrow) == 1 % if one parent was older
                                child_1(9 + 1 : 9 + 2) = child_tags(findrow,:);
                                child_2(9 + 1 : 9 + 2) = child_tags(findrow,:);
                            end
                            light1 = test_strategy_range_qs(child_1,param_range);
                            light2 = test_strategy_range_qs(child_2,param_range);
                        end
                        % if only one child if ok, make the other one the same
                        if light1 == 0
                            child_1 = child_2;
                        end
                        if light2 == 0
                            child_2 = child_1;
                        end
                        % plant children
                        child_place_1 = top_f*n+(children_places-1)*2 + 1;
                        child_place_2 = child_place_1 + 1;
                        WARR_next(child_place_1,:) = child_1;
                        WARR_next(child_place_1,10) = WARR_next(child_place_1,10) + 1; % tagPool +1
                        WARR_next(child_place_2,:) = child_2;
                        WARR_next(child_place_2,10) = WARR_next(child_place_2,10) + 1; % tagPool +1
                    end % make child_f*n cross-over children
                    for mutant_place = (top_f*n + child_f*n + 1) : (n - n*frac_mig) % make n*(1 -top_f - child_f - frac_mig) mutants
                        % select the mother for mutatant children depending on
                        % fitness
                        mother_strat = WARR(randsample(n,1,true,FITS),:);
                        which_genes_mutate = (mother_strat(1:9) ~= 0);
                        gene_mut = which_genes_mutate; % only mutate genes
                        %gene_mut = randi(7); % randi(7) SOME PARAMS FIXED?;
                        which_sens_func = 2*randi(4); % choose which sensor part mutates
                        rand_mut = rand;
                        % decide the type of mutation
                        if rand_mut < 0.8 % tiny tiny steps
                            mut_type = 1;
                        elseif rand_mut < 0.8 + 0.2  % small steps
                            mut_type = 2;
                        else % HGT dont let this happen
                            mut_type = 3;
                        end
                        mut_strat_test = create_mutationHGT_range_qs(gene_mut,mut_type,param_range,which_sens_func,mother_strat); % mutation happens
                        light = test_strategy_range_qs(mut_strat_test,param_range);
                        while light == 0 % redo mutation until it works
                            mut_strat_test = create_mutationHGT_range_qs(gene_mut,mut_type,param_range,which_sens_func,mother_strat); % mutation happens
                            light = test_strategy_range_qs(mut_strat_test,param_range);
                        end
                        WARR_next(mutant_place,:) = mut_strat_test;
                        WARR_next(mutant_place,10) = WARR_next(mutant_place,10) + 1; % tagPool +1
                    end % make n*(1 -top_f - child_f - frac_mig) mutants
                    for migrant_place = (n - n*frac_mig + 1) : n % make n*frac_mig migrants
                        % first select sensing for this guy
                        random_guy_sensors_i = [0 0 0 0];
                        if sum(random_guy_sensors) >= 2
                            random_guy_sensors_i(randsample(find(random_guy_sensors ==1),1)) = 1 ;
                        else
                            random_guy_sensors_i = random_guy_sensors;
                        end
                        % then create the random guy until it works
                        mut_strat_test = create_random_guy_range_qs(param_range,f_all,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
                        light = test_strategy_range_qs(mut_strat_test,param_range);
                        while light == 0
                            mut_strat_test = create_random_guy_range_qs(param_range,f_all,random_guy_sensors_i); % create_random_guy % CHANGE ALL PARAMS?
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
            savename = [rpath name name_sens num2str(alg_runs) 't' num2str(endround) landing_tag '.mat']; % zero small much
            WARR_saveS = WARR_save;
            endroundS = endround;
            nS = n;
            frac_migS = frac_mig; n_fightsS = n_fights; top_fS = top_f;
            WARRS = WARR; child_fS = child_f;
            save(savename,'WARR_saveS','endroundS','nS','frac_migS','WARRS',...
                'n_fightsS','top_fS','child_fS');
            % save again for R
            savename2 = [rpath name name_sens num2str(alg_runs) 't' num2str(endround) landing_tag '_onlyWARR.mat']; % zero small much
            save(savename2,'WARR_saveS');
            % Plot how the mean values evolve
            if plotting_0 == 1 % plotting
                % useful stuff
                gensize = 1:endround;
                save_average = mean(WARR_save(1:(n-n*frac_mig),1:9,:),1);
                save_average = reshape(save_average,[9,endround]);
                save_variance = var(WARR_save(1:(n-n*frac_mig),1:9,:),1);
                save_variance = reshape(save_variance,[9,endround]);
                % Defaults values
                width = 15;     % Width in cm
                height = 15;    % Height in inches
                alw = 0.75;    % AxesLineWidth
                fsz = 20;      % Fontsize
                lw = 2.5;      % LineWidth
                msz = 8;       % MarkerSize
                % create a white figure
                figure;
                pos = get(gcf, 'Position'); hold on
                set(gcf, 'Position', [pos(1) pos(2) width*39.37, height*39.37]); %<- Set size
                set(gca, 'FontSize', fsz, 'LineWidth', alw, 'FontName','Arial'); %<- Set properties
                %%%%%%%%%%% graph settings %%%%%%%%%%%%%%%%
                hold on % [fnull +-fN UN +-fTB UTB +-fTA UTA]
                fill([gensize fliplr(gensize)],  [(save_average(1,:) - save_variance(1,:)) fliplr(save_average(1,:)+ save_variance(1,:))],...
                    [0    0    0]/255,'EdgeColor','w','FaceAlpha', 0.3);
                plot1 = plot(gensize , save_average(1,:), 'LineWidth', lw,'Color',[0    0    0]/255); % Plots main line
                fill([gensize fliplr(gensize)],  [(save_average(2,:) - save_variance(2,:)) fliplr(save_average(2,:)+ save_variance(2,:))],...
                    [0    0    88]/255,'EdgeColor','w','FaceAlpha', 0.3);
                plot2 = plot(gensize , save_average(2,:), 'LineWidth', lw,'Color',[0    0    88]/255); % Plots main line
                fill( [gensize fliplr(gensize)],  [(save_average(4,:) - save_variance(4,:)) fliplr(save_average(4,:)+ save_variance(4,:))],...
                    [0    88    88]/255,'EdgeColor','w','FaceAlpha', 0.3);
                plot3 = plot(gensize , save_average(4,:), 'LineWidth', lw,'Color',[0    88    88]/255); % Plots main line
                fill( [gensize fliplr(gensize)],  [(save_average(6,:) - save_variance(6,:)) fliplr(save_average(6,:)+ save_variance(6,:))],...
                    [88    0    0]/255,'EdgeColor','w','FaceAlpha', 0.3);
                plot4 = plot(gensize , save_average(6,:), 'LineWidth', lw, 'Color',[88    0    0]/255); % Plots main line
                set(plot1,'DisplayName','f_{Null}');
                set(plot2,'DisplayName','f_{N}');
                set(plot3,'DisplayName','f_{Topp}');
                set(plot4,'DisplayName','f_{Tself}');
                % Limit axis
                %xlim([0 3]);
                %ylim([-0.4 0.4]);
                % Create xlabel and ylabel
                xlabel('Loops of algorithm');
                % Create legend
                if alg_runs == 1 % make legend only for 1st graph
                    legend1 = legend([plot1 plot2 plot3 plot4]);
                    set(legend1,'EdgeColor',[1 1 1],'Location','Best','box','off');
                end
                % Set Tick Marks
                %set(gca,'XTick',0:0.002:0.01);
                %set(gca,'YTick',0:2000:14000);
                % Here we preserve the size of the image when we save it.
                set(gcf,'InvertHardcopy','on');
                set(gcf,'PaperUnits','inches');
                papersize = get(gcf, 'PaperSize');
                left = (papersize(1)- width)/2;
                bottom = (papersize(2)- height)/2;
                myfiguresize = [left, bottom, width, height];
                set(gcf,'PaperPosition', myfiguresize);
                set(gcf, 'PaperPositionMode', 'auto');
                figname = [name num2str(alg_runs)];
                print(figname ,'-dpng')
                %
            end
        end
    end
end
%%%%%%%%%%%%%% the end

