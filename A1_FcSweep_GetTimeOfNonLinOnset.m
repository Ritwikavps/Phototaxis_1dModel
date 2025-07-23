clearvars; clc

%Ritwika VPS; Jule 2025
% This script recursively analyses front profiles for various parameter combos to determine time corresponding to non-linearity onset to determine appropriare simulation times so that
% roughly equal number of fronts are available for Fourier analysis and growth rate extraction for all parameter combos. A challenge here is that even as actual simulated time decreases,
% the time increment (dt, dimensionless) does not change, which means that we need to scale dt such that all simulations are run for roughly the same number of time steps. This is 
% necessary because time evolution happens quicker for higher fc values (for instance) which means that 'more' evolution happens within a time step for higher fc values vs lower fc values.
% Therefore, to get roughly the same amount of information for different simulations, we first determine (recursivly) the non-linearity onset with a generous peak height-to-width 
% threshold (no lower than 0.3), scale the time increment such that simulations for the corresponding non-linearity onset times for all parameter combos have approximately the same
% number of time steps.
% 
% The proposed algorithm is as follows:
%   - initialise the set of simulations (for parameter value combos in Sweep_Cmean, Sweep_tau, and Sweep_Fc, which corresponds to regimes verified for finger formation) with 
%       the default total time simulated (600000 s) and default dimensionless time increment (dt = 0.1)
%   - determine non-linearity onset time with a finger height-to-width threhsold of 0.5 as the cutoff, for each parameter combo
%   - set the new total time simulated (say, T_actual_iter) as the set of non-linearity onset times estimated
%   - set the max value of the above set of new total time simulated (T_actual_iter) as the 'standard' to estimate number of time steps to scale time increments for all other
%       parameter combos
%   - repeat till the difference in estimated non-linearity onset times for two consecutive sweps is less than 20 minutes (since ramp times tend to be in hours; see Fig, 7 in the
%       Motility paper, Ursell etal, 2013)
%   - Optionally, once this goal reached for Ht-to-Width cutoff = 0.5, repeat with lower thresholds to get better estimates.
%   - get associated time step increment and actual time (in seconds) to simulate
% 
% Note that the threhsolds set in this script may have to be adapted for other parameter combos currently not being invetsigated. These numbers have been set based on looking at
% relevant time scales genreated, etc.
% 
% This script takes roughly 1.5 hours or so.

%Here, we use a sweep over fc values (non-dimensionalised) and 3 sets of Cmean and tau values (also non-dimensionalised) that were tested and confirmed as corresponding to finger
% formation.
paramsweep_Fc_nondim = [1 2 3 4 5 6]; %values of fc (non-dimensionalised) to sweep through
Cmean_nondim_Vals = [2^-6 2^-6 2^-7]; tau_nondim_Vals = [2^7 2^6 2^5]; %different pairs of cmean_nondim and tau_nondim to have distinct sweeps of fc

% Put these values into vectors so that parfor can more efficiently parse. In this version, a parameter value combo is indexed as a set, i.e., Cmean_nondim_Vals(i) and 
% tau_nondim_Vals(i) go together (see for loop below).
Sweep_Fc = [paramsweep_Fc_nondim paramsweep_Fc_nondim paramsweep_Fc_nondim];
Sweep_Cmean = []; Sweep_tau = [];
for i = 1:numel(Cmean_nondim_Vals)
    Sweep_Cmean = [Sweep_Cmean Cmean_nondim_Vals(i)*ones(size(paramsweep_Fc_nondim))];
    Sweep_tau = [Sweep_tau tau_nondim_Vals(i)*ones(size(paramsweep_Fc_nondim))];
end

%set values for the initial parfor to generate the first set of non linearity onset times (with default dt and T_actual).
Init_Tactual = 600000*ones(numel(Sweep_Cmean),1); %Initial actual simulated time
%dt = 0.1; %dimensionless time increment in sims (default)
Init_dt_DivFactor = 1*ones(numel(Sweep_Cmean),1); %Default factor by which the time increment is divided (this will later be amended to get appopriate 
% simulation step numbers for various parameter combos)
graphout = 0;
Curr_HtToWdthThresh = 0.4; %Height-to-width threhsold for non-linearity onset initialised
TimeThreshToStop = 200; %seconds; if non linear onset times from two consecutive iterative trials differ by less than this amount, then that parameter combo is not tested again

%Do for ht-to-wdth thresh = 0.4
[TrialCtr_0_4, dt_DivFactor_0_4, NonLinOnsetTimes_cell_0_4, Curr_NonLinOnsetTimeAct_0_4] = ...
                GetNonLinOnsetTimes_ForHtToWdthThresh(Sweep_Cmean, Sweep_Fc, Sweep_tau, Init_Tactual, Init_dt_DivFactor,Curr_HtToWdthThresh, TimeThreshToStop);


%repeat for threshold 0.3 (with a more generous time difference tolerance, cuz 0.3 is a lower threhold for ht-to-wdth).
Curr_HtToWdthThresh = 0.3; %Height-to-width threhsold for non-linearity onset initialised
TimeThreshToStop = 600; %seconds; if non linear onset times from two consecutive iterative trials differ by less than this amount, then that parameter combo is not tested again

[TrialCtr_0_3, dt_DivFactor_0_3, NonLinOnsetTimes_cell_0_3, Curr_NonLinOnsetTimeAct_0_3] = ...
                GetNonLinOnsetTimes_ForHtToWdthThresh(Sweep_Cmean, Sweep_Fc, Sweep_tau, Curr_NonLinOnsetTimeAct_0_4, dt_DivFactor_0_4, Curr_HtToWdthThresh, TimeThreshToStop);

%write results to file
Destination_Path = '~/Desktop/GoogleDriveFiles/research/phototaxis/FreshAttempt2025/'; cd(Destination_Path)
Fname = 'NonLinOnsetTimes_Iterated_NonLinHtToWdthThresh_0_4__0_3.mat';

%Put parameter values in the structure.
NonLinOnsetTimeStruct(1).Fc_nondim = Sweep_Fc;  NonLinOnsetTimeStruct(1).Tau_nondim = Sweep_tau; NonLinOnsetTimeStruct(1).Cmean_nondim = Sweep_Cmean; 
NonLinOnsetTimeStruct(2).Fc_nondim = Sweep_Fc;  NonLinOnsetTimeStruct(2).Tau_nondim = Sweep_tau; NonLinOnsetTimeStruct(2).Cmean_nondim = Sweep_Cmean; 

%results for HtToWdthThresh = 0.4
NonLinOnsetTimeStruct(1).HtToWdthThresh = 0.4; 
NonLinOnsetTimeStruct(1).NonLinOnsetTimesCell_secs = NonLinOnsetTimes_cell_0_4; 
NonLinOnsetTimeStruct(1).EstimatedFinalSimTime = Curr_NonLinOnsetTimeAct_0_4; 

NonLinOnsetTimeStruct(2).HtToWdthThresh = 0.3;
NonLinOnsetTimeStruct(2).NonLinOnsetTimesCell_secs = NonLinOnsetTimes_cell_0_3;
NonLinOnsetTimeStruct(2).EstimatedFinalSimTime = Curr_NonLinOnsetTimeAct_0_3; 

save(Fname,'NonLinOnsetTimeStruct')

%% Below, code from before functionalising the parallel iterative simulation time estimation based on ht-to-width non linearity threshold.
% p = parpool(8); %initialise paralle ppool
% TimeAct_NonLinOnset = NaN*ones(numel(Sweep_Cmean),1); %initalise arrays/vectors to store non-linearity onset time for each iteration of the sweep
% 
% %Do the first sweep to get initial numbers
% parfor i_sweep = 1:numel(Sweep_Cmean) %go through parameter combos for the sweep
%     [ZeroContouryCell,TimeActVecCell,XvecCell,~, ~, ~, ~, ~] = ...
%                                 GetPhotoFronts_w_SlimeDecay(Sweep_Cmean(i_sweep),Sweep_Fc(i_sweep),Sweep_tau(i_sweep), ...
%                                 graphout, Init_Tactual, Init_dt_DivFactor);%TimeAct_NonLinOnset_05_2(i_sweep)); %0 for graphout
% 
% 
%     [TimeAct_NonLinOnset(i_sweep), ~] =... %NonLinOnsetIndex(i_sweep)
%                 Get_HtToWdthThreshCutOffTime(ZeroContouryCell, XvecCell, TimeActVecCell,Curr_HtToWdthThresh);
% end
% 
% 
% % %% %Plot to visualise the times decreasing: 2 sets of colours for the two cut off thresholds used.
% % red = [1, 0.1, 0.1]; %red
% % pink = [0.1, 1, 0.1]; %green
% % colors_p_1 = [linspace(red(1),pink(1),20)', linspace(red(2),pink(2),20)', linspace(red(3),pink(3),20)'];
% % colors_p_1 = [colors_p_1; colors_p_1(end:-1:1,:); colors_p_1; colors_p_1(end:-1:1,:)]; %stack this back to front a few times in case there are many iteration
% % 
% % figure('Color',[1 1 1]); hold on %to plot the consecutive simulation times
% % %plot initial siumulation time (600000 s)
% % scatter(1:numel(TimeAct_NonLinOnset),Init_Tactual*ones(size(TimeAct_NonLinOnset)), 40, 'filled','MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6);
% % %first set of non linear onset times
% % scatter(1:numel(TimeAct_NonLinOnset),TimeAct_NonLinOnset, 40, 'filled','MarkerFaceColor', colors_p_1(1,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6); 
% % 
% % 
% % %% %Now resume compuatation after plotting
% RefMaxTime = max(TimeAct_NonLinOnset); %The max value of the list of nonlinearity onset times. This serves as a refernce to scale the simulation time increments (dt) for
% %all parameter combos such that there are roughly similar number of simulation steps for all parameter combos.
% 
% Curr_NonLinOnsetTimeAct = TimeAct_NonLinOnset; %set the outputed nonlinearlity onset time as the previous time vector
% Prev_NonLinOnsetTimeAct = Init_Tactual*ones(numel(Sweep_Cmean),1); %Set the initial T_actual value as the previous nonlinear onset time vector (whiile this is not fully accurate, 
% % we refer the Curr_NonLinOnsetTimeAct with this previous list of times to determine if there has been a huge update in the simulation time (at each update, we use the most 
% % recent nonlinear onset time as the full simulation time for the next trial for that parameter combo).
% clear TimeAct_NonLinOnset %and clear the variable
% 
% TrialCtr = 1; %count the number of trials/iterations of the sweep
% NonLinOnsetTimes_cell{TrialCtr} = Curr_NonLinOnsetTimeAct; %store the vectors of non linearity onset
% HtToWdthThresh_Vec(TrialCtr) = Curr_HtToWdthThresh;
% 
% while true %repeat till break condition
% 
%     TrialCtr = TrialCtr + 1 %keep track of number of trials
%     Curr_Prev_Tact_diff = abs(Curr_NonLinOnsetTimeAct - Prev_NonLinOnsetTimeAct); %get the difference between the cureent non linear onset time and the previous iteration
%     MaxTime_ToCurrTimeRatio = round(RefMaxTime./Curr_NonLinOnsetTimeAct,1); %divide the refernce max time by the list of non linearity onset times (and round to oen decimal place)
%     dt_DivFactor = round(MaxTime_ToCurrTimeRatio);%get the factor to divide dt by so that we have roughly the same number of timesteps for all sims (with
%     % RefMaxTime as the refernce)
% 
%     TimeAct_NonLinOnset = NaN*ones(numel(Sweep_Cmean),1); %initialise
%     Prev_NonLinOnsetTimeAct = Curr_NonLinOnsetTimeAct; %recaste the 'current' vector (from the last iteration) as previous
% 
%     parfor i_sweep = 1:numel(Sweep_Cmean) %go through indiced for fc_nondim values for the sweep
%         %Each structure for a given Cmean-tau combo is indexed by fc values (SimRunStruct(i_sweep), where i_sweep indexes fc). For each fc value, there are results from 8 trials.
%         % That is, SimRunStruct(i_sweep).YContours{i_trial}, etc. Note that each Ycontours{i_trial} is itself a cell array storing the y values of the fronts at each simulation time
%         % point, and thgerefore, each Ycontours{i_trial} is indexed by the simulation time point: Ycontours{i_trial}{i_simtimepoint}. All other fields are either vectors or scalars
%         % (SimRunStruct(i_sweep).XvedGridVals{i_trial}, etc are vectors, while SimRunStruct(i_sweep).Fc_nondim are scalars and are not indexed further with i_trial).
% 
%         %Get the difference between the current T_actual and the previous T_actual for this parameter combo. 
% 
%         if Curr_Prev_Tact_diff(i_sweep) > 200 %if the difference is more than 600 s (10 min), do it again; note that this threhsold may have to be adapted for other
%             %parameter combos currently not being invetsigated. The threhsold number.
%             [ZeroContouryCell,TimeActVecCell,XvecCell,~, ~, ~, ~, ~] = ...
%                                         GetPhotoFronts_w_SlimeDecay(Sweep_Cmean(i_sweep),Sweep_Fc(i_sweep),Sweep_tau(i_sweep), ...
%                                         graphout, Curr_NonLinOnsetTimeAct(i_sweep),dt_DivFactor(i_sweep));%TimeAct_NonLinOnset_05_2(i_sweep)); %0 for graphout
% 
% 
%             [TimeAct_NonLinOnset(i_sweep), ~] =... %NonLinOnsetIndex(i_sweep)
%                         Get_HtToWdthThreshCutOffTime(ZeroContouryCell, XvecCell, TimeActVecCell,Curr_HtToWdthThresh);
% 
%         elseif Curr_Prev_Tact_diff(i_sweep) <= 200 %if the difference is less
%             TimeAct_NonLinOnset(i_sweep) = Curr_NonLinOnsetTimeAct(i_sweep); %copy the same value in; this is not getting modified in this iteration
%         end
%     end
% 
%     % %plot most recent non linear onset times.
%     % scatter(1:numel(TimeAct_NonLinOnset),TimeAct_NonLinOnset, 40, 'filled','MarkerFaceColor', colors_p_1(TrialCtr,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.6); 
% 
%     RefMaxTime = round(max(TimeAct_NonLinOnset(:)),-2); %round to nearest 100
%     Curr_NonLinOnsetTimeAct = TimeAct_NonLinOnset; %set the outputed nonlinearlity onset time as the current time vector
% 
%     %store vector of non linear onset times as well as corresponding ht-to-wdth thresholds used in computtaion.
%     NonLinOnsetTimes_cell{TrialCtr} = Curr_NonLinOnsetTimeAct; %store the vectors of non linearity onset
%     HtToWdthThresh_Vec(TrialCtr) = Curr_HtToWdthThresh;
% 
%     %if the difference between the current and previous nonlinearity onset vectors are less than 600 s (10 min) for every parameter combo, we break.
%     if (all(abs(Curr_NonLinOnsetTimeAct-Prev_NonLinOnsetTimeAct) < 600*ones(size(TimeAct_NonLinOnset))))
%         disp('hello')
%         break
%     end 
% end
% 
% hold off %remove hold on plot
% delete(p);
% 
% t_elap = toc;
% disp(['Computation took: ' num2str(t_elap/60/60) ' hrs'])

