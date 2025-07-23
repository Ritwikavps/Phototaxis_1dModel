clear; clc

%Ritwika VPS; 2025
%This script loads .mat files (containing structs) from specified simulation runs, estimates non-linearity onset time(s) based on set thresholds for feature height-to-width
% ratio thresholds in the fronts, does FFT and computes one-sided amplitude spectra for all fronts (and removes the DC component), and saves all this info. 

%In thsi iteration, where fc is teh parameter swept over, For each Struct_w_Data, we have the same values of tau and Cmean_ini (both nondim) constant. Struct_w_Data is indexed
% by fc (nondim). Within each substruct (corresponding to a single combo of fc, tau, and Cmean_ini), we have 8 trials. 

tic

DataPath = '~/Desktop/GoogleDriveFiles/research/phototaxis/FreshAttempt2025';
cd(DataPath)

Dir_FcSweep = dir('Photo2dSimsSweepFcNondim_*.mat'); %read in all relevant .mat files (which are results of sweeping over fc)

p = parpool(8); %set up parallel pool: same number of workers as there are trials for each combo of fc x tau x Cmean_ini (all non dim)

for i1_struct = 1:numel(Dir_FcSweep) %go through list of mat files

    CurrStruct = load(Dir_FcSweep(i1_struct).name); %load current mat file
    Struct_w_Data = CurrStruct.SimRunStruct; %get the structure
    % Struct_w_data is indexed by fc, and each element of Struct_w_data is a structure containing details of 8 trials for the same combo of fc, Cmean_ini, and tau. 
    % So, within each Struct_w_data(i), we have cell arrays indexed by the trial number with details about each trial sim run for that parameter combo.

    %CHECKS!!
    u_tauVals = unique([Struct_w_Data.tau_nondim]);
    u_CmeanVals = unique([Struct_w_Data.Cmean_nondim]);
    u_fcVals = unique([Struct_w_Data.Fc_nondim]);
    
    %Because each .mat file corresponds to a single sweep across fc values with tau and Cmean_ini values the same across the entire sweep, the number of unique tau vals and the
    % number of unique Cmean_ini vals should be 1, while the number of unique fc values should be the same as the number of fc values. Note, the struct is indexed by the fc value.
    if (numel(u_tauVals) ~= 1) || (numel(u_CmeanVals) ~= 1) || (numel(u_fcVals) ~= numel(Struct_w_Data))
        error('Number of unique parameters is not as expected for the sim rin. JAIL!')
    end

    for i_fc = 1:numel(Struct_w_Data)  %go through each fc value
        
        Substruct_i = Struct_w_Data(i_fc); %get the substrcutre corresponding to the i-th fc value
        NumTrials = numel(Substruct_i.TimeVals); %get the number of trials; TimeVals is a cell array with each element being a vector with all time values in 
        %that trial sim run

        if NumTrials ~= 8 %check for number of trials
            warning(['Number of trials for fc = ' num2str(Substruct_i.Fc_nondim) ', Cmean_ini = '  num2str(Substruct_i.Cmean_nondim) ...
                ', tau = ' num2str(Substruct_i.tau_nondim) ' (all nondim) is not equal to 8'])
        end
    
        %Inigialise cell arrays to store outputs for FFT (upto the last time point in simulation, for all fronts that have non-empty FFT amp peaks).
        All_kVals_wPks_temp = cell(NumTrials,1);
        AmpPks_cell_Processed_temp = cell(NumTrials,1);
        Time_wPks_s_temp = cell(NumTrials,1);
    
        %Inigialise vectors to store time and index corresponding to non linearity onset for various finger height-to-width thresholds.
        TimeAct_NonLinOnset_0_3_temp = zeros(NumTrials,1); 
        NonLinOnsetIndex_0_3_temp = zeros(NumTrials,1);
        TimeAct_NonLinOnset_0_2_temp = zeros(NumTrials,1); 
        NonLinOnsetIndex_0_2_temp = zeros(NumTrials,1);
        TimeAct_NonLinOnset_0_25_temp = zeros(NumTrials,1); 
        NonLinOnsetIndex_0_25_temp = zeros(NumTrials,1);
    
        parfor j_trial  = 1:NumTrials %parallelise across trials
    
            Xvec = Substruct_i.XvecGridVals{j_trial}; %get the x grid vec
            Fronts_Y = Substruct_i.YContours{j_trial}; %get the cell array with fronts for all times in sim
            TimeVec = Substruct_i.TimeVals{j_trial}; %get vector with all times in sim (in s)
    
            %get FFT outpits
            [All_kVals_wPks_temp{j_trial},AmpPks_cell_Processed_temp{j_trial},Time_wPks_s_temp{j_trial}] = GetFourierSpectra(Xvec,Fronts_Y,TimeVec);
    
            %get non linearity onset details for various thresholds
            [TimeAct_NonLinOnset_0_3_temp(j_trial), ~, ~] = Get_HtToWdthThreshCutOffTime(Fronts_Y,Xvec,TimeVec,0.3);
            [TimeAct_NonLinOnset_0_25_temp(j_trial), ~, ~] = Get_HtToWdthThreshCutOffTime(Fronts_Y,Xvec,TimeVec,0.25);
            [TimeAct_NonLinOnset_0_2_temp(j_trial), ~, ~] = Get_HtToWdthThreshCutOffTime(Fronts_Y,Xvec,TimeVec,0.2);
            [TimeAct_NonLinOnset_0_15_temp(j_trial), ~, ~] = Get_HtToWdthThreshCutOffTime(Fronts_Y,Xvec,TimeVec,0.15);
            [TimeAct_NonLinOnset_0_1_temp(j_trial), ~, ~] = Get_HtToWdthThreshCutOffTime(Fronts_Y,Xvec,TimeVec,0.1);
            
        end
    
        %Store in output structure indexed by fc: FFT outputs
        FFTStruct(i_fc).All_kVals_wPks = All_kVals_wPks_temp; 
        FFTStruct(i_fc).AmpPks_cell_Processed = AmpPks_cell_Processed_temp;
        FFTStruct(i_fc).Time_wPks_s = Time_wPks_s_temp;
        FFTStruct(i_fc).AllTimeVals = Substruct_i.TimeVals; %copy all time values (including ones for which fronts do not have non-empty FFT amp peaks) from Substrcut_i
    
        %Store in output structure indexed by fc: non linearity onset outputs
        FFTStruct(i_fc).TimeAct_NonLinOnset_0_3 = TimeAct_NonLinOnset_0_3_temp;
        FFTStruct(i_fc).TimeAct_NonLinOnset_0_25 = TimeAct_NonLinOnset_0_25_temp;
        FFTStruct(i_fc).TimeAct_NonLinOnset_0_2 = TimeAct_NonLinOnset_0_2_temp;
        FFTStruct(i_fc).TimeAct_NonLinOnset_0_15 = TimeAct_NonLinOnset_0_15_temp;
        FFTStruct(i_fc).TimeAct_NonLinOnset_0_1 = TimeAct_NonLinOnset_0_1_temp;

        %Store in output structure indexed by fc: parameters
        FFTStruct(i_fc).Fc_nondim = Substruct_i.Fc_nondim;
        FFTStruct(i_fc).Cmean_nondim = Substruct_i.Cmean_nondim;
        FFTStruct(i_fc).tau_nondim = Substruct_i.tau_nondim;
    end

    tdata = datetime; %get date and time info to save files if needed
    fdate = [date '_' num2str(tdata.Hour) '-' num2str(tdata.Minute)];
    FileNameToSave = ['FftOneSidedAmpPksDcRem_Photo2dSimsSweepFcNondim_CmeanNondim_' num2str(Substruct_i.Cmean_nondim) ...
        '_TauNondim_' num2str(Substruct_i.tau_nondim) '_8Trials_dtAndTactScaled_' fdate '.mat'];
    save(FileNameToSave,'FFTStruct','-v7.3') %large files
    clear FFTStruct
end

delete(p)
t_elap = toc;
disp(['Time elapsed = ' num2str(t_elap/60) ' mins'])


