clear; clc

%Ritwika VPS, July 2025
%Processing code for results from sweep over fc values for 3 sets of tau and Cmean_ini values (from FcSweepSims.m).
%Code is broke into sections so that the data structure need only be loaded once. Recommended to run things in each required section so that the data structure is preserved after initial
%loading into the workspace.

DataPath = '~/Desktop/GoogleDriveFiles/research/phototaxis/FreshAttempt2025/';
cd(DataPath)

DataDir = dir('Photo2dSimsSweepFcNondim_*.mat'); %dir fc sweep results (might have to change the string in function call later)
for i = 1:numel(DataDir)
    TempDat = load(DataDir(i).name);
    DataStruct(i).SweepResults = TempDat.SimRunStruct; %put all data into a structure

    %Now, DataStruct has as many elements as the number of .mat files that were loaded in, and has one field name SweepResults (which is itself a structure).
    % So, DataStruct(1).SweepResults will access data from the first .mat file loaded, etc.
    % SweepResults is a 1x6 structure (as of now, for the 6 fc values investigated), with each element corresponding to one combo of fc_nondim, with values of tau_nondim and 
    % Cmean_nondim remaining the same for each SweepResults. So:
    % DataStruct (indexed by the number of pairs of tau_nondim and Cmean_nondim investigated; DataStruct(i).SweepResults):
    %       SweepResults (indexed by number of fc_nondim values investigated within the current tau_nondim-Cmean_nondim pair; DataStruct(i).SweepResults(i_fc).YContours):
    %             YContours, etc (indexed by the trial number for the current fc_nondim-Cmean_nondim-tau_nondim combo; DataStruct(i).SweepResults(i_fc).YContours{i_trial}):
    %                   Each YContours{i_trial} is itself a cell arry, so needs to be indexed as YContours{i_trial}{i_time}, since each i_time-th element is the front profile a time 
    %                   point in the simulation. All other non-parameter fields of sweep results are cell arrays (XvecGridVals, TimeVals, etc) and each of them contains vectors. 
    %                       So, for example, XvecGridVals does not need to be indexed further
    %                   The files of sweep results that are parameter values are scalars and do not need to be indexed by the i_time: DataStruct(i).SweepResults(i_fc).Fc_nondim
    %                       (and similarly for Cmean_nondim and tau_nondim).                    
end

%% %Exploratory plots: Testing how finger size varies with fc.

clearvars -except DataStruct %Clearing DataStruct will reuiqre re-loading it.

FrontsTimeEvoPlot = true; %to determine if we want to plot time evolution of fronts for each trial for each fc-tau-Cmean combo
inds = 1; %do one index (indexing the combo of tau and Cmean) at a time
TempDat = DataStruct(inds).SweepResults; %get the structure corresponding to the relevant combo

% %Uncomment below loop as necessaryFor loop just to confirm the indexing scheme.
% for i = 1:numel(TempDat)
%     [TempDat(i).Cmean_nondim TempDat(i).Fc_nondim TempDat(i).tau_nondim]
% end

%plot contour time evolution (and plot all trials for a given fc-tau-Cmean combo to ensure that the trials all show similar patterns for the same parameter combo).
if FrontsTimeEvoPlot
    for i = 1%:numel(TempDat) %TempDat is the structure with results from fc sweep (8 trials) for one combo of tau and Cmean.
        %i indexes the fc value.

        figure('Color',[1 1 1]); %Each figure has a single fc value
        for j = 1:numel(TempDat(i).YContours) %j indexes the trial for a given fc-tau-Cmean combo
            subplot(3,3,j); hold all %plot each trial as its own subplot
            for k = 1:300:numel(TempDat(i).YContours{j}) %only plot every 300th front
                plot(TempDat(i).XvecGridVals{1},TempDat(i).YContours{j}{k},'k')
            end

            if j == 2
                title(['\tau (non-dim) = ' num2str(TempDat(1).tau_nondim) ', C_{mean, ini} (non-dim) = ' num2str(TempDat(1).Cmean_nondim) ', fc = ' num2str(TempDat(i).Fc_nondim)])
            end
            axis tight
        end
    end
end

%plot final fronts only: how finger size varies with fc.
figure('Color',[1 1 1]); 

%colour gradient
red = [1, 0, 1]; %magenta
pink = [0, 1, 1]; %cyan
colors_p = [linspace(red(1),pink(1),numel(DataStruct(1).SweepResults))',...
    linspace(red(2),pink(2),numel(DataStruct(1).SweepResults))',... 
    linspace(red(3),pink(3),numel(DataStruct(1).SweepResults))'];

for i = 1:numel(DataStruct) %loop through each combo of Cmean and tau (one subplot for each)
    subplot(1,3,i); hold all
    for j_fc = 1:numel(DataStruct(i).SweepResults) %indexes fc value %plotting all final fronts for a given Cmean-tau combo in one subplot)
        plot(DataStruct(i).SweepResults(j_fc).XvecGridVals{1},DataStruct(i).SweepResults(j_fc).YContours{1}{end},'Color',colors_p(j_fc,:),'LineWidth',2)
        %This plots the last conc profile for the first trial (YContours{1}) for each combo of Fc_nondim (j) and tau_nondim-Cmean_nondim (i).
        %XvecGridVals is the same for all indices (cuz it is the just the grid values for x grid) so just index with 1. Ycountour is indexed with 1 because within each 
        %DataStruct(i).SweepResults(j_fc), Ycontour{ind} is such that ind indexes the trial number. Then pick the last contour in the YContours cell (YSontours{1}{end}).
        axis tight
        legendCell{j_fc} = ['fc (non-dim) = ' num2str(DataStruct(i).SweepResults(j_fc).Fc_nondim)];
    end
    %tau and Cmean are the same within each DataStruct(i).
    title(['\tau (non-dim) = ' num2str(DataStruct(i).SweepResults(1).tau_nondim) ', C_{mean, ini} (non-dim) = ' num2str(DataStruct(i).SweepResults(1).Cmean_nondim)])
    legend(legendCell); clear legendCell
    grid on; grid minor
end

%% %Exploratory plots: Effects of different non-linearity cut-offs

clearvars -except DataStruct; %Clearing DataStruct will reuiqre re-loading it.

HtToWdthThresh = [1 0.5 0.4 0.3 0.2 0.1]; %set various thresholds for non-linearity cut off based on peak height to width ratio

figure('Color',[1 1 1]); hold all
%colour gradient
red = [1, 0, 1]; %magenta
pink = [0, 0, 0]; %cyan
colors_p = [linspace(red(1),pink(1),numel(HtToWdthThresh))',...
    linspace(red(2),pink(2),numel(HtToWdthThresh))',... 
    linspace(red(3),pink(3),numel(HtToWdthThresh))'];

FcInds = [1 3 6]; %for lowest, middle and highest value of fc_nondim: to show the effect of different ht-to-width thresholds

%(Testing only for the first tau-Cmean combo in DataStruct)
%There is a separate sub-plot for each fc value in the first tau-Cmean combo DataStruct, and within each subplot, different colours for each nonlinearity threshold.
for i_fc = 1:numel(FcInds) %go through vairous fc values we want to explore (indexed by FcInds)
    subplot(1,3,i_fc); hold all %one subplot for each fc index

    fc_ind = FcInds(i_fc); %get the index corresponding to the intended fc value

    for ii = 1:numel(HtToWdthThresh) %go through non-linearity threshold values
        %get nonlinearity onset flag vector for the i_fc-th simulation and ii-th non-linearity cut off (we use the first of the 8 trials (YContours{1}) for the current fc value; 
        % the tau-Cmean combo is set as that for the first data structure in the series (DataStruct(1))).
        clear NonLinOrNo
        for i = 1:numel(DataStruct(1).SweepResults(fc_ind).YContours{1}) %for each front (YContours{1}{i}) in the current simulation
            NonLinOrNo(i) = FrontNonLinTest_PkWdthsAndHts(DataStruct(1).SweepResults(fc_ind).XvecGridVals{1}, ...
                                                DataStruct(1).SweepResults(fc_ind).YContours{1}{i}, HtToWdthThresh(ii));
           
        end

        NonLinOnsetIndex = NaN; %Default

        %going backwards from final front, stop at first linear front.
        for i = numel(NonLinOrNo):-1:1
            if NonLinOrNo(i) == 0
                NonLinOnsetIndex = i;
                break
            end
        end

        %Store the NonLinOrNo vector in a cell array: row indices are for the fc value, column indices are for the height-to-width non-linear threshold cutoff value.
        % Because we are only looking at a subset of the fc values, for future retrieval, I am indexing it corresponding to each fc value in the original fc vec.
        % So, while FcInds = 1 corresponds to fc = 1, and will be indexed by a row index of 1. But FcInds = 3 corresponds to the third fc value (fc = 3), and therefore, will be
        % indexed by a row index of 3, and FcInds = 6 will be indexed a row index of  (corresponding to fc = 6).
        % The column indices are similarly indexed but since we are exploring all values of HtToWdthThresh, this indexing is more straightforward.
        % But, essentially, to get the NonLinOrNo vector for a given fc and HtToWdthThresh, the NonLinOrNo_Cell needs to be indexed by the indices for fc and HtToWdthThresh
        % as they correspond to the entire fc and HtToWdthThresh vectors.
        NonLinOrNo_Cell{fc_ind,ii} = NonLinOrNo; 
        NonLinOnsetIndex_Array(fc_ind,ii) = NonLinOnsetIndex; %similarly, store the NonlinOnsetIndex in an array, same indexing as NonLinOrNo_Cell

        plot(DataStruct(1).SweepResults(fc_ind).XvecGridVals{1}, DataStruct(1).SweepResults(fc_ind).YContours{1}{NonLinOnsetIndex},'Color',colors_p(ii,:),'LineWidth',2) %plot
        if i_fc == 1 %only do legend if this is the first subplot
            LegendCell{ii} = ['Ht-to-Wdth < ' num2str(HtToWdthThresh(ii))]; %get legend entries
        end
    end
    
    if i_fc == 1 %only do legend if this is the first subplot
        legend(LegendCell)
    end
    title(['\tau (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).tau_nondim) ', C_{mean, ini} (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).Cmean_nondim) ...
            ', f_c (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).Fc_nondim)])
    axis tight; grid on; grid minor  
end


%% Exploratory plots: demo of initial fronts getting tagged as non linear for fc = 1 vs fc = 6, when HtToWdthThresh = 0.1, 0.5.

clearvars -except DataStruct NonLinOrNo_Cell NonLinOnsetIndex_Array; %Clearing DataStruct will reuiqre re-loading it (Same for the others).

%Still only looking at the first tau-Cmean combo (DataStruct(1))

%NonLinOrNo_Cell and NonLinOnsetIndex_Vec are indexed such that the row index indexes the fc value, and the column index indexes the HtToWdthThresh value.
HtToWdthThresh = [1 0.5 0.4 0.2 0.1]; %full vector
HtToWdthThresh_candidates = [0.5 0.1]; %required values
HtToWdthThresh_reqinds = [2 5]; %indices corresponding to required values

FcInds = [1 6]; %for lowest and highest value of fc_nondim: to show the effect of different ht-to-width thresholds
%For HtToWdthThresh, we want values 0.5, 0.1. So that is indices 2, 5. For fc, we want 1 and 6, which are indices 1 and 6.

figure('Color',[1 1 1]); hold all

%plot to show that a lot of initial fronts also get tagged as non-linear for low height-to-width nonlinearity thresholds.
%Working with fronts tagged as non-linear ONLY for fc = 1 vs fc = 6, HtToWdthThresh = 0.1, 0.5 
Ctr = 0; %counter variable for subplots
for i_fc = 1:numel(FcInds) %go through vairous fc values
    
    fc_ind = FcInds(i_fc); %get the index corresponding to the intended fc value

    for i_nonlin = 1:numel(HtToWdthThresh_candidates)

        Ctr = Ctr + 1;
        subplot(2,2,Ctr); hold all
        nonlin_ind = HtToWdthThresh_reqinds(i_nonlin); %get the index corresponding to the intended threshold value
        
        %get required vector of non linear flags and index forresponding to nonlinearity onset (see note in previous section indexing NonLinOrNo_Cell and NonLinOnsetIndex_Vec).
        CurrNonLinVec = NonLinOrNo_Cell{fc_ind,nonlin_ind};
        CurrNonLinOnsetInd = NonLinOnsetIndex_Array(fc_ind,nonlin_ind);

        %set increment to plot fronts depending on how mnay fronts to plot. Set increment such that only 100 fronts are plotted (the hope is that some fronts tagged as nonlinear
        % and some tagged as linear, before the final determined non-linearity onset will be represented.
        PlotInc = round(CurrNonLinOnsetInd/100); 

        for i = 1:PlotInc:CurrNonLinOnsetInd %go through the last vector of nonlinearity of no till the nonlinearity onset index
            if CurrNonLinVec(i) == 0
                %plot red if front not tagged as non-linear.
                plot(DataStruct(1).SweepResults(fc_ind).XvecGridVals{1},DataStruct(1).SweepResults(fc_ind).YContours{1}{i},'r','LineWidth',2)
            elseif CurrNonLinVec(i) == 1
                %plot grey if tagged as non linear
                plot(DataStruct(1).SweepResults(fc_ind).XvecGridVals{1},DataStruct(1).SweepResults(fc_ind).YContours{1}{i},'Color',[0.5 0.5 0.5])
            end
        end

        %plot last front before non-linearity onset per the cutoff threshold.
        plot(DataStruct(1).SweepResults(fc_ind).XvecGridVals{1},DataStruct(1).SweepResults(fc_ind).YContours{1}{CurrNonLinOnsetInd},'r','LineWidth',2)

        title(['\tau (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).tau_nondim) ', C_{mean, ini} (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).Cmean_nondim) ...
            ', f_c (non-dim) = ' num2str(DataStruct(1).SweepResults(fc_ind).Fc_nondim) ', height to width ratio cut-off = ' num2str(HtToWdthThresh(nonlin_ind))])
        axis tight
    end
end

clear all
        



