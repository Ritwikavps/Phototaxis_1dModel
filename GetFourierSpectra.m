function [All_kVals_wPks,AmpPks_cell_Processed,Time_wPks_s] = GetFourierSpectra(Xvec,Fronts_Y,TimeVec)

%This function takes in y coords of fronts (Fronts_Y) detected from a given simulation run along with a time vector (Timevec; in seconds) corresponding to the fronts and the X grid 
% vector (Xvec) correspoding to the fronts and does FFT on each front at each time point to get the following:
    % - All_kVals_wPks: Vector of all k values that have FFT amp peaks in the (one sided FFT amp after removing DC amp) data. That is, if a k value has an amp peak in at least 
                      % one front, that k value will be in this vector, in sorted order. DC k value is removed.
    % - AmpPks_cell_Processed: the cell array with one sided FFT amplitude spectra peaks (after DC amplitude removed) for each front, corresponding to the k values in All_kVals_wPks.
                             % If a front does not have a peak corresponding to a k value in All_kVals_wPks, then that position has an NaN value. All elements in AmpPks_cell_Processed
                             % are row vectors of the same size, and AmpPks_cell_Processed is a column cell array. Note that if a front has no FFT amp peaks, it is not included in
                             % AmpPks_cell_Processed.
    % - Time_wPks_s: time vector (in s) corresponding to the fronts that are represented in AmpPks_cell_Processed. 


%% Get FFT for all fronts

delta = min(diff(Xvec)); %sampling bin width
fkc = 1/delta; %sampling frequency; also sets the Nyquist frequency as fkc/2
% For the fourier analysis, the sampling in the data is at roughly 5 um (at least for simulation with decay when the bias force is the parameter being varied; this changes
% if m0 or k_s are varied), while the features we are interested in are around ~100 to 200 um, so my hunch is that we don't need to interpolate for a higher sampling frequency 
% and we are safe from worrying about alisaing (ie., we are sampling at above the Nyquist frequency). We can revisit this if necessary.

clear kvals_forPks_cell Amps_forPks_cell OneSidedAmp_DcRem_cell k_DcRem_cell GrowthRate

NonEmptyPks_Ctr = 0;

for i = 1:numel(Fronts_Y)

    % We are going to concatenate the i/p signal a few times to extend the sampling window so that frequency bin resolution is increased and to decrease spectral leakage. To do this,
    % I am assuming the fronts are provided as row vectors. Check for that and get the transpose if necessary.
    if size(Fronts_Y{i},1) ~= 1
        Ydata = (Fronts_Y{i})';
    else
        Ydata = Fronts_Y{i};
    end

    %Ydata = smooth(Ydata);
    DataForFFt = [Ydata Ydata Ydata Ydata Ydata Ydata]; %concatenate 6x; this means that the number of samples is always even
    fft_y = fft(DataForFFt); %get FFT; This is a row vector

    N = numel(DataForFFt); %number of sample points

    %Get one sided spectra
    Inds = 1:(N/2 + 1); %1st value is the DC component (0 freq) and the last one is the Nyquist
    %For odd N: Inds = 1:((N-1)/2 + 1); %Nyquist freq is not included, so the doubling for amplitude scaling should only exclude the first freq
    %(but we remove that anyway because we do not want to get the DC signal)
    OneSidedSpec = fft_y(Inds);
    OneSidedAmp = abs(OneSidedSpec)/N; %scale
    OneSidedAmp(2:end-1) = 2*OneSidedAmp(2:end-1); %double everything excepyt the Nyquist and DC; This is a row vector

    %Get freqs
    k = (0:N-1)*fkc/N; %This is for the whole 2-sided spectrum; note that this is not quite accurate because the second half would be negative frequencies, but we are ignoring
    %that for now, because we will only be taking the positive side anyway
    k = k(Inds); %Get only the one-sided (+ve freqs)

    OneSidedAmp_DcRem = OneSidedAmp(2:end); %remove DC comp; This is a row vector
    k_DcRem = 2*pi*k(2:end);%remove DC comp and multiply by 2pi for wave vector value; This is a row vector

    % Both OneSidedAmp_DcRem_Cell and k_DcRem_Cell are column cell arrays and each element is a row vector (Nx1 cell array, where N = number of time points for which
    % we have fronts): for OneSidedAmp_DcRem_Cell, the i-th vector in the cell array has the one sided FFT amplitude after removing the DC component element for the i-th
    % simulation time point. Similarly, for k_DcRem_Cell, the i-th vector has all k values corresponding to the one sided FFT amplitude after removing the DC k value. Note that
    % because we intentionally have an even numbered input sample, this k vector includes k = half of fkc (the sampling frequency). This is the Nyquist frequency for the
    % sampling rate used.
    OneSidedAmp_DcRem_cell{i,1} = OneSidedAmp_DcRem; %This is a column cell array and each element is a row vector
    k_DcRem_cell{i,1} = k_DcRem; %remember taht the first k after removing DC corresponds to the size of the system, but here the system has been 
    %concatenated 6 times (in DataForFFt)

    [pks,locs] = findpeaks(OneSidedAmp_DcRem); %get FFT amplitude peaks; pks and locs are both row vectors
    if (~isempty(pks)) && (~isempty(locs)) %if pks and locs are non-empty, get corresponding k values and store both in cell arryas; 
        NonEmptyPks_Ctr = NonEmptyPks_Ctr + 1; %increment counter

        % Both kvals_forPks_cell and Amps_forPks_cell are column cell arrays and each element is a row vector (Nx1 cell array, where N = number of time points for which
        % we have non empty vectors for FFT peaks): For both these cell arrays, each element is a row vector. However, these row vectors may not have the same size because
        % not all time points may have all frequencies present as peaks.
        kvals_forPks_cell{NonEmptyPks_Ctr,1} = k_DcRem(locs);  %k values corresponding to peaks in FFT amp spectra (DC k removed)
        Amps_forPks_cell{NonEmptyPks_Ctr,1} = pks; %amplitudes for peaks in FFT amp spectra
        Time_wPks_s(NonEmptyPks_Ctr,1) = TimeVec(i); %time corresponding to the non empty FFT peak vector; in seconds
    end   
end


%% Checks!!

if unique(OneSidedAmp_DcRem_cell{1}) ~= 0 %The front for the first time point should be a stright line and therefore, only have a DC component (which we removed), 
    % so all other FFT amplitudes should be 0
    error('The FFT amplitude spectra (after removing DC comp) for the t = 0 front is not empty.')
end

if unique(diff(Fronts_Y{1})) ~= 0 %The front for the first time point should be a stright line 
    error('The front at t = 0 is not a flat, straight line.')
end


%% Process peaks and corresponding k values: 

% the goal is to get the set of unique k values that have peaks in the FFT amplitude for at least one front, and then for each front, assign corresponding amplitude values 
% to k values that have peaks for that front, and NaN values otherwise. So, say at time t = 1, k = [1, 2, 3, 4] has peaks, and for t = 3, k = [1, 3, 10] has peaks. Then, 
% kvals_forPks_cell{1} = [1 2 3 4] and kvals_forPks_cell{3} = [1 3 10] (note that these are not necessarily actual k values in the data, this is simply meant as an example).
% Then, in this processing, we get the total set of unique k values = [1 2 3 4 10] (All_kVals_wPks), and for each time point, we have corresponding amplitude vectors.
% So, AmpPks_cell_Processed{1} = [amp_1 amp_2 amp_3 amp_4 amp_10=NaN] and we have AmpPks_cell_Processed{3} = [amp_1 amp_2=NaN amp_3 amp_4=NaN amp_10].
% This way, we get a cell array for the amplitude where each element of the cell array is a row vector of the same size.

All_kVals_wPks = unique(cell2mat(kvals_forPks_cell')); %get all unique k values in the entire data corresponding to peaks; transpose because kvals_forPks_cell is a column cell array
%with each element being a row vector, but all these row vectors are not necessarily the same size
Num_ukVals_wPks = numel(All_kVals_wPks); %get number of unique k values

for i = 1:numel(Amps_forPks_cell) %go through cell array storing FFT amp peaks for each time point that has a non-empty vector for these peaks
    
    %get ith peaks and kvals vectors
    Pks_i = Amps_forPks_cell{i}; kvals_i = kvals_forPks_cell{i};
    
    %Below, in doing intersect this way (and this order in the intersect function call is important), we get the indices of the k values in kvals_i as they appear in
    % All_kVals_wPks (inds_All_k).
    [~,~,inds_All_k] = intersect(kvals_i,All_kVals_wPks,'stable');

    NewPks_temp = NaN*ones(1,Num_ukVals_wPks); %NaN row vector with the same number of elements as there are unique k values corresponding to all peaks in the current sim run.
    NewPks_temp(inds_All_k) = Pks_i; %Put the peaks for the current peak amplitude vector into the vector corresponding to all unique k values in the data. This way,
    %if the current front does not have a peak for a specific k, that will remain NaN.
    AmpPks_cell_Processed{i,1} = NewPks_temp;

    %get test vectors to make sure that the indexing through intersect is correct. If done correctly, these test vectors after removing the NaN values should be equal to 
    % the ith (unprocessed) amp and kvals vectors.
    kvals_Test = All_kVals_wPks(~isnan(NewPks_temp)); Pks_Test = NewPks_temp(~isnan(NewPks_temp));

    if ~isequal(kvals_Test,kvals_forPks_cell{i}) || ~isequal(Pks_Test,Amps_forPks_cell{i})
        disp('Indexing with intersect incorrect. JAIL!')
    end
end

%More checks: Number of unique k values across the entire dataset (All_kVals_wPks) should be equal to the length of every row vector in AmpPks_cell_Processed. AND the number
% of elements in AmpPks_cell_Processed should be the same as the number of time points in Time_wPks_s.
if (numel(All_kVals_wPks) ~= unique(cellfun(@numel,AmpPks_cell_Processed))) || (numel(AmpPks_cell_Processed) ~= numel(Time_wPks_s))
    error('Inconsistent numel() results. JAIL!')
end

% %%
% %Include a check for empty peak vectors etc
% for i = 3:numel(Pks_k) %check if all the peaks are at the same location 
%     %exclude the first front cuz that is a flat front (at least when based
%     %on Tristan's saving code cuz the first time step is saved
%     if ~isequal(Pks_k{2},Pks_k{i})
%         i
%     end
% end
% 
% figure;
% hold all
% for i = 2:20:numel(TimeVec)
%     stem(k_DcRem_Cell{i},OneSidedAmp_DcRem_Cell{i})
% end
% 
% Amps = cell2mat(Pks(2:end));
% kvals = Pks_k{2};
% RelTimes = TimeVec(2:end);
% 
% red = [1, 0, 1]; %magenta
% pink = [0, 1, 1]; %cyan
% colors_p = [linspace(red(1),pink(1),numel(kvals))', linspace(red(2),pink(2),numel(kvals))', linspace(red(3),pink(3),numel(kvals))'];
% 
% figure; hold all
% for i = 1:numel(kvals)
%     Amp_wTime_ForCurrk = Amps(:,i);
% 
%     if ~any(Amp_wTime_ForCurrk) %if there are no non-zero values
% 
%         Amp_wTime_ForCurrk(Amp_wTime_ForCurrk == 0) = NaN;
% 
%         [f1,gof] = fit(RelTimes',log(Amp_wTime_ForCurrk),'m*x + c');
%         GrowthRate(i) = f1.m;
% 
%         plot(RelTimes,Amp_wTime_ForCurrk,'Color',colors_p(i,:))
%     end
% end
% 
% figure;
% hold all
% plot(kvals,GrowthRate)



