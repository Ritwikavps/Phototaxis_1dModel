function [TrialCtr, dt_DivFactor, NonLinOnsetTimes_cell, Curr_NonLinOnsetTimeAct] = ...
                GetNonLinOnsetTimes_ForHtToWdthThresh(Sweep_Cmean, Sweep_Fc, Sweep_tau, Init_Tactual, Init_dt_DivFactor,Curr_HtToWdthThresh, TimeThreshToStop)

% Ritwika VPS, July 2025
% This function takes in an initialised vector of simulation time (Init_Tactual, in seconds) and corresponding scaling factors (Init_dt_DivFactor) to divide the default dimensionless 
% time increment (dt = 0.1) in the 2d simulations so that all simulations have roughly similar number of simulation steps, corresponding to a list of Cmean_ini, Fc, and 
% tau (Sweep_Cmean, Sweep_Fc, Sweep_tau; all non-dimensionalised); runs simulations iteratively to determine time to simulate the 2d sims for each parameter combo based on a specified
% ht-to-wdth non linearity threshold (Curr_HtToWdthThresh); and outputs the number of trials taken to iterate to get to a list of simulation times subject to the nonlinearity threshold 
% specified (TrialCtr), the cell array containing determined simulation times for each parameter combo for each iteration (NonLinOnsetTimes_cell), and the list of scale factors to
% divide the default dt by for the last set of actual simulation times for each parameter combo (dt_DivFactor) as well as the nonlinear onset times for the last set of simulations 
% carried out before the while loop is terminated (Curr_NonLinOnsetTimeAct). Makes use of parallelisation.
% 
% This function does the above recursively, and the idea is to 'converge' to simulation times for each parameter combo such that simulations are only carried out till a specified
% non-linearity finger ht-tp-wdth threshold *and* making sure that all parameter combos have roughly simlar number of simulation steps. This is because for sims that progress quickly, 
% each time step contains more evolution, such that sims that evolve more slowly have more fine grained front details (which becomes relevant for the FFT). The iterative procedure allows
% for these attempts to generate better estimates of the required simulation times based on a specified ht-to-wdh threshold because we can iterate through more precise estimates with 
% appropriately scaled division factors for the default dt. For instance, for a simulation that progresses really quickly, an initial estimate of when non linearity sets in will necessarily
% be coarse grained, because the default dt = 0.1 allows for each dt to contain a lot of time evolution. Once this initial estimate is done, we can scale dt to be smaller such that the
% amount of evolution is reduced per time ste. Repeating this subject to some specified tolerance (where the iteration is stopped for a given parameter combo if current and previous
% estimates of the nonlinearity onset time differs by less than a specified amount of time) allows to get better and better estimates of this non-linearity onset time.
% 
% Note that it is important to set the nonlinearoty onset threshold to higher than what you would use to cut off admitting fronts for FFT within the linear regime because we reliably want
% there to be some evolution beyond the linear regime in our final simulations so that we can get all fronts in the linear regime for the FFT vs the iterative determination of the 
% nonlinearity onset time resulting in some fronts that would have been deemed linear not being reached in the simulation (because th initialised Cmean has some randomness and no 2 sims
% are the same). So, if this iterative esimaton is based on the nonlinearity ht-to-wdth threshold that will be used for FFT cutoff, there is a chance that some trials may attain nonlinearity
% after the estimated cutoff. 

    tic 

    graphout = 0; %no graph output for 2d sims in parallel loop
    p = parpool(8); %initialise paralle ppool

    TimeAct_NonLinOnset = NaN*ones(numel(Sweep_Cmean),1); %initalise arrays/vectors to store non-linearity onset time for each iteration of the sweep
    
    %Do the first sweep to get initial numbers
    parfor i_sweep = 1:numel(Sweep_Cmean) %go through parameter combos for the sweep
        [ZeroContouryCell,TimeActVecCell,XvecCell,~, ~, ~, ~, ~] = ...
                                    GetPhotoFronts_w_SlimeDecay(Sweep_Cmean(i_sweep),Sweep_Fc(i_sweep),Sweep_tau(i_sweep), ...
                                    graphout, Init_Tactual(i_sweep), Init_dt_DivFactor(i_sweep));
        
        %get time of nonlinearity onset for the specified finger ht-to-wdth thresolhd.
        [TimeAct_NonLinOnset(i_sweep), ~, ~] = Get_HtToWdthThreshCutOffTime(ZeroContouryCell, XvecCell, TimeActVecCell, Curr_HtToWdthThresh);
    end
    
    
    % %% %Plot to visualise the times decreasing: 2 sets of colours for the two cut off thresholds used.
    % red = [1, 0.1, 0.1]; %red
    % pink = [0.1, 1, 0.1]; %green
    % colors_p_1 = [linspace(red(1),pink(1),20)', linspace(red(2),pink(2),20)', linspace(red(3),pink(3),20)'];
    % colors_p_1 = [colors_p_1; colors_p_1(end:-1:1,:); colors_p_1; colors_p_1(end:-1:1,:)]; %stack this back to front a few times in case there are many iteration
    % 
    % figure('Color',[1 1 1]); hold on %to plot the consecutive simulation times
    % %plot initial siumulation time (600000 s)
    % scatter(1:numel(TimeAct_NonLinOnset),Init_Tactual*ones(size(TimeAct_NonLinOnset)), 40, 'filled','MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6);
    % %first set of non linear onset times
    % scatter(1:numel(TimeAct_NonLinOnset),TimeAct_NonLinOnset, 40, 'filled','MarkerFaceColor', colors_p_1(1,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6); 
    % 
    % 
    % %% %Now resume compuatation after plotting
    RefMaxTime = max(TimeAct_NonLinOnset); %The max value of the list of nonlinearity onset times. This serves as a refernce to scale the simulation time increments (dt) for
    %all parameter combos such that there are roughly similar number of simulation steps for all parameter combos.
    
    Curr_NonLinOnsetTimeAct = TimeAct_NonLinOnset; %set the outputed nonlinearlity onset time as the previous time vector
    Prev_NonLinOnsetTimeAct = Init_Tactual; %Set the initial T_actual value as the previous nonlinear onset time vector (whiile this is not fully accurate, 
    % we refer the Curr_NonLinOnsetTimeAct with this previous list of times to determine if there has been a huge update in the simulation time (at each update, we use the most 
    % recent nonlinear onset time as the full simulation time for the next trial for that parameter combo).
    clear TimeAct_NonLinOnset %and clear the variable
    
    TrialCtr = 1; %count the number of trials/iterations 
    NonLinOnsetTimes_cell{TrialCtr} = Curr_NonLinOnsetTimeAct; %store the vectors of non linearity onset
    
    while true %repeat till break condition
    
        TrialCtr = TrialCtr + 1 %keep track of number of trials
        Curr_Prev_Tact_diff = abs(Curr_NonLinOnsetTimeAct - Prev_NonLinOnsetTimeAct); %get the difference between the cureent non linear onset time and the previous iteration
        MaxTime_ToCurrTimeRatio = round(RefMaxTime./Curr_NonLinOnsetTimeAct,1); %divide the refernce max time by the list of non linearity onset times (and round to oen decimal place)
        dt_DivFactor = round(MaxTime_ToCurrTimeRatio);%get the factor to divide dt by so that we have roughly the same number of timesteps for all sims (with
        % RefMaxTime as the refernce)
    
        TimeAct_NonLinOnset = NaN*ones(numel(Sweep_Cmean),1); %initialise
        Prev_NonLinOnsetTimeAct = Curr_NonLinOnsetTimeAct; %recaste the 'current' vector (from the last iteration) as previous
    
        parfor i_sweep = 1:numel(Sweep_Cmean) %go through indiced for fc_nondim values for the sweep
            %Each structure for a given Cmean-tau combo is indexed by fc values (SimRunStruct(i_sweep), where i_sweep indexes fc). For each fc value, there are results from 8 trials.
            % That is, SimRunStruct(i_sweep).YContours{i_trial}, etc. Note that each Ycontours{i_trial} is itself a cell array storing the y values of the fronts at each simulation time
            % point, and thgerefore, each Ycontours{i_trial} is indexed by the simulation time point: Ycontours{i_trial}{i_simtimepoint}. All other fields are either vectors or scalars
            % (SimRunStruct(i_sweep).XvedGridVals{i_trial}, etc are vectors, while SimRunStruct(i_sweep).Fc_nondim are scalars and are not indexed further with i_trial).
        
            %Get the difference between the current T_actual and the previous T_actual for this parameter combo. 
            
            if Curr_Prev_Tact_diff(i_sweep) > TimeThreshToStop %if the difference is more than a set threhold, do it again; note that this threhsold may have to be adapted for other
                %parameter combos currently not being invetsigated. The threhsold number.
                [ZeroContouryCell,TimeActVecCell,XvecCell,~, ~, ~, ~, ~] = ...
                                            GetPhotoFronts_w_SlimeDecay(Sweep_Cmean(i_sweep),Sweep_Fc(i_sweep),Sweep_tau(i_sweep), ...
                                            graphout, Curr_NonLinOnsetTimeAct(i_sweep),dt_DivFactor(i_sweep));%TimeAct_NonLinOnset_05_2(i_sweep)); %0 for graphout
            
            
                [TimeAct_NonLinOnset(i_sweep), ~, ~] =... %NonLinOnsetIndex(i_sweep)
                            Get_HtToWdthThreshCutOffTime(ZeroContouryCell, XvecCell, TimeActVecCell,Curr_HtToWdthThresh);
    
            elseif Curr_Prev_Tact_diff(i_sweep) <= TimeThreshToStop %if the difference is less
                TimeAct_NonLinOnset(i_sweep) = Curr_NonLinOnsetTimeAct(i_sweep); %copy the same value in; this is not getting modified in this iteration
            end
        end
    
        % %plot most recent non linear onset times.
        % scatter(1:numel(TimeAct_NonLinOnset),TimeAct_NonLinOnset, 40, 'filled','MarkerFaceColor', colors_p_1(TrialCtr,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6); 
    
        RefMaxTime = round(max(TimeAct_NonLinOnset(:)),-2); %round to nearest 100
        Curr_NonLinOnsetTimeAct = TimeAct_NonLinOnset; %set the outputed nonlinearlity onset time as the current time vector
    
        %store vector of non linear onset times as well as corresponding ht-to-wdth thresholds used in computtaion.
        NonLinOnsetTimes_cell{TrialCtr} = Curr_NonLinOnsetTimeAct; %store the vectors of non linearity onset
    
        %if the difference between the current and previous nonlinearity onset vectors are less than a set threshold for every parameter combo, we break.
        if (all(abs(Curr_NonLinOnsetTimeAct-Prev_NonLinOnsetTimeAct) < TimeThreshToStop*ones(size(TimeAct_NonLinOnset))))
        disp(['This is just a check that we have reached the break condition for ht-to-wdth threshold' num2str(Curr_HtToWdthThresh)])
            break
        end 
    end
    
    delete(p);
    
    t_elap = toc;
    disp(['Computation took: ' num2str(t_elap/60/60) ' hrs for ht-to-wdth threshold' num2str(Curr_HtToWdthThresh)])
    input(['Please make sure that all text printed to console makes physical sense. If yes, press any key to proceed. The console (but not workspace memort) ' ...
        'will be cleared once you provide an input'],'s')
    clc