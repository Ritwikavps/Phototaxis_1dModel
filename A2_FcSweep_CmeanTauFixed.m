clear; clc

%Ritwika VPS, Jul 2025
%Parameter sweep for 2d simulations to test predictions about the bias force. Here, we do 3 sets of combos for the values of the slime decay constant and the initial total bacterial
% concentration (tau and C_mean_ini, both non-dimensionalised). For each such combo, we have 6 values of the bias force, fc (also non-dimensionalised), which is the parameter we
% are doing the sweep over. For each combo of Cmean, fc, and tau, we perform 8 trials so that noise in the extracted growth curve (growth rate, lambda, vs wave numer k) can be averaged
% out. These parameter value combos were chosen by visually inspecting simulation outputs for entire simulation length for a larger set of parameters. Selection criteria involved finger 
% formation for the all fc values in the set of fc values below (which themselves are centered around the fc value used in Fig 9 in the Motility enhancement paper, Ursell et al, 2013).
% 
% Because we are interested in the linear regime where fingers haven't grown very large, we utilise non-linearity onset times estimated for a finger ht-to-width threshold of 0.3 to set 
% both the amount of time each simulation (for each parameter combo) runs (in s). We also scale the simulation time step (default = 0.1 dimensionless time units) such that all sims
% for all parameter combos investigated here have roughly similar number of time steps so that the growth rate fitting for each fourier mode is done using roughly the same number of data
% points.

paramsweep_Fc_nondim = [1 2 3 4 5 6]; %values of fc (non-dimensionalised) to sweep through
Cmean_nondim_Vals = [2^-6 2^-6 2^-7]; tau_nondim_Vals = [2^7 2^6 2^5]; %different pairs of cmean_nondim and tau_nondim to have distinct sweeps of fc
%These are indexed as a set, i.e., Cmean_nondim_Vals(i) and tau_nondim_Vals(i) go together (see for loop below_.

SavePath = '~/Desktop/GoogleDriveFiles/research/phototaxis/FreshAttempt2025/'; %path to save files; also has non linearity onset time data (see FcSweep_GetTimeOfNonLinOnset.m)
cd(SavePath) 

TempStruct = load('NonLinOnsetTimes_Iterated_NonLinHtToWdthThresh_0_4__0_3.mat'); %load structure with the non linear onset times (index 1 corresponds to a ht-to-width 
% threshold of 0.4, and index 2 corresponds to a threshold of 0.3)
NLt_struct_0_3 = TempStruct.NonLinOnsetTimeStruct([TempStruct.NonLinOnsetTimeStruct.HtToWdthThresh] == 0.3); %pick out the substruct corresponding to 0.3 threshold

%We choose to use the 0.3 threshold times by adding an hour of extra simulation time to these times (as a buffer).
SimTimes = NLt_struct_0_3.NonLinOnsetTimesCell_secs{end}; %get times corresponding to the last iteration (See FcSweep_GetTimeOfNonLinOnset.m for details) 
dt_DivVec = ceil(max(SimTimes)./(SimTimes)); %get the most conservative factor to divide dt by (we add the hour to SimTimes later, because getting the dt division factor 
% gives a slightly higher division factor)
SimTimes = SimTimes + 3600; %add the extra hour buffer
Inds = 1:numel(SimTimes); %get indices for teh SimTimes vector (to match for the correct fc, cmean, and tau values)

Num_ParForLoops = 8; %number of loops for parfor (different trials for each fc_nondim + fixed Cmean_nondim/tau_nondim pairs in a sweep))
p = parpool(Num_ParForLoops); %initialise parallel pool

for i_nonsweep = 1:numel(Cmean_nondim_Vals) %go through the indices for the fixed pairs of cmean_nondim and tau_nondim for each fc_nondim sweep
    %Results from each i_nonsweep iteration are saved as one structure file. That is, each distinct structure file (saved as .mat files) has a single value of Cmean and a single value
    % of tau.

    i_nonsweep

    for i_sweep = 1:numel(paramsweep_Fc_nondim) %go through indiced for fc_nondim values for the sweep
        %Each structure for a given Cmean-tau combo is indexed by fc values (SimRunStruct(i_sweep), where i_sweep indexes fc). For each fc value, there are results from 8 trials.
        % That is, SimRunStruct(i_sweep).YContours{i_trial}, etc. Note that each Ycontours{i_trial} is itself a cell array storing the y values of the fronts at each simulation time
        % point, and thgerefore, each Ycontours{i_trial} is indexed by the simulation time point: Ycontours{i_trial}{i_simtimepoint}. All other fields are either vectors or scalars
        % (SimRunStruct(i_sweep).XvedGridVals{i_trial}, etc are vectors, while SimRunStruct(i_sweep).Fc_nondim are scalars and are not indexed further with i_trial).

        ZeroContouryCell = cell(Num_ParForLoops,1); %initialise cell arrays to store results in for each trial (to pass to parfor)
        TimeVecCell = cell(Num_ParForLoops,1);
        XvecCell = cell(Num_ParForLoops,1);
        YvecCell = cell(Num_ParForLoops,1);
        Final_C_Cell = cell(Num_ParForLoops,1);

        %find index for the SimTimes vector corresponding to the fc, cmean, and tau values being investigated.
        SimTime_index = Inds(NLt_struct_0_3.Fc_nondim == paramsweep_Fc_nondim(i_sweep) & ...
                    NLt_struct_0_3.Cmean_nondim == Cmean_nondim_Vals(i_nonsweep) & NLt_struct_0_3.Tau_nondim == tau_nondim_Vals(i_nonsweep));

        parfor i_par = 1:Num_ParForLoops %do sims for each trial using parfor
            [ZeroContouryCell{i_par},TimeVecCell{i_par},XvecCell{i_par},YvecCell{i_par},Final_C_Cell{i_par}, ~, ~, ~] = ...
                GetPhotoFronts_w_SlimeDecay(Cmean_nondim_Vals(i_nonsweep),paramsweep_Fc_nondim(i_sweep),tau_nondim_Vals(i_nonsweep), ...
                                                0, SimTimes(SimTime_index), dt_DivVec(SimTime_index)); %0 for graphout
    
        end

        %store results for the sweep just performed in a structure
        SimRunStruct(i_sweep).YContours = ZeroContouryCell;
        SimRunStruct(i_sweep).TimeVals = TimeVecCell;
        SimRunStruct(i_sweep).XvecGridVals = XvecCell;
        SimRunStruct(i_sweep).YvecGridVals = YvecCell;
        SimRunStruct(i_sweep).FinalConcProfile = Final_C_Cell;
    
        %store relevant parameters
        SimRunStruct(i_sweep).Fc_nondim = paramsweep_Fc_nondim(i_sweep);
        SimRunStruct(i_sweep).Cmean_nondim = Cmean_nondim_Vals(i_nonsweep);
        SimRunStruct(i_sweep).tau_nondim = tau_nondim_Vals(i_nonsweep);
    end

    tdata = datetime; %get date and time info to save files if needed
    fdate = [date '_' num2str(tdata.Hour) '-' num2str(tdata.Minute)];
    FileNameToSave = ['Photo2dSimsSweepFcNondim_CmeanNondim_' num2str(Cmean_nondim_Vals(i_nonsweep)) ...
        '_TauNondim_' num2str(tau_nondim_Vals(i_nonsweep)) '_8Trials_dtAndTactScaled_' fdate '.mat'];
    save(FileNameToSave,'SimRunStruct','-v7.3') %large files

    clear SimRunStruct
end

delete(p);





