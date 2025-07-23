function [GrowthRate, GrowthRatePreFac, FitRsq, NumPtsFit] = GetGrowthRates(All_kVals_wPks,AmpPks_cell_Processed,Time_wPks_s,NonLinOnsetTime)

%This function takes in the vector with all k values that have peaks in the given dataset (All_kVals_wPks), the cell array containing vectors with FFT amp peaks for each
% simulation time point (AmpPks_cell_Processed; corresponds to All_kVals_wPks, with NaN for k values that do not have a peak in the current front), and the time values in 
% seconds (Time_wPks_s) corresponding to the fronts represented by the row vectors in AmpPks_cell_Processed, and computes the growth rate for each k value for the entire time
% represented in Time_wPks_s as well as the time of nonlinearity onset (NonLinOnsetTime).
% 
% Outputs: - GrowthRate: Vector of growth rates correspinding to k values in All_kVals_wPks (provided there are at least 10 time points to fit. If not, the corresponding
                       % growth rate is NaN)
         % - GrowthRatePreFac: the prefactor for the fitted growth rate ('a' in a exp(bx)); same as above re: NaN and non-NaN values
         % - FitRsq: r squared value of fit; same as above re: NaN vs non-NaN values
         % - NumPtsFit: number of non-NaN amplitude values that in the amp vector for each k value
         % (All outputs are row vectors).

Amps_mat = cell2mat(AmpPks_cell_Processed); %plop all amplitudes for all times into a mat; this would be a matric of dimensions (numel(Time_wPks_s)) x (numel(All_kVals_wPks)). 
% The rows are indexed by the simulation time point, while the columsns are indexed by the k value. That is, each row i_t of Amps_all carries amplitude info for time i_t and 
% each column i_k carries the amplitude info for wave vector value k. 

%Checks!
if height(Amps_mat) ~=  numel(Time_wPks_s)
    error('Num. rows in the cell2mat-d FFT peak amplitudes is not the same as the number of time points')
end

if width(Amps_mat) ~= numel(All_kVals_wPks)
    error('Num. cols in the cell2mat-d FFT peak amplitudes is not the same as the number of k values that have FFT amp peaks in the data.')
end

if NonLinOnsetTime > Time_wPks_s(end)
    error('Provided onset of nonlinearity time is greater than the total simulation time')
end

TimeInds = 1:numel(Time_wPks_s); %get indices for the time vector

%subset everything only up to nonlinearity onset
if NonLinOnsetTime ~= Time_wPks_s(end)
    TimeInds = TimeInds(Time_wPks_s <= NonLinOnsetTime);
    Time_wPks_s = Time_wPks_s(TimeInds);
    Amps_mat = Amps_mat(TimeInds,:); 
end

%Ibitialise vectors to store fitting outputs
GrowthRate = NaN*ones(size(All_kVals_wPks)); 
GrowthRatePreFac = NaN*ones(size(All_kVals_wPks));
FitRsq = NaN*ones(size(All_kVals_wPks));
NumPtsFit = NaN*ones(size(All_kVals_wPks));

for i_k = 1:numel(All_kVals_wPks) %plot growth profiles for each k
    Amp_ForCurrk = Amps_mat(:,i_k); %get amoplitude profiles for the current k value

    TimesForFit = Time_wPks_s(~isnan(Amp_ForCurrk)); %get times for non-nan amps
    AmpsForFit = Amp_ForCurrk(~isnan(Amp_ForCurrk)); %get non-nan amps

    %check that inputs to fit func are column vectors and convert if not
    if (height(AmpsForFit) == 1) && numel(AmpsForFit > 1)
        AmpsForFit = AmpsForFit';
    end

    if (height(TimesForFit) == 1) && numel(TimesForFit > 1)
        TimesForFit = TimesForFit';
    end

    NumPtsFit(i_k) = numel(AmpsForFit);

    if numel(AmpsForFit) > 10 %if there are at least 20 non-zero points to fit
        [f1,gof] = fit(TimesForFit,AmpsForFit,'exp1'); %do an exponential fit (a little more computationally intensive but more accurate, esp if there are varying 
        %error profiles across the amplitude as a fucntion of time

        GrowthRate(i_k) = f1.b; %get growth rate for the current k value
        GrowthRatePreFac(i_k) = f1.a;
        FitRsq(i_k) = gof.rsquare; %r sq of fit
    end   
end