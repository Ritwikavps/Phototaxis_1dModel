function [TimeAct_NonLinOnset, NonLinOnsetIndex, NonLinOrNo] = Get_HtToWdthThreshCutOffTime(FrontsCell,Xvec,TimeActVec,HtToWdthThresh)

%Ritwika VPS, July 2025
%This function evaluates front profiles for a given simulation to estimate the actual time (in s) of non-linearity onset.
% Inputs: - FrontsCell:  cell array containing the fronts as a function of actual time elapsed in the simulation (TimActVec)
%         - Xvec: the vector of X grid points
%         - TimeActVec: the vector of actual time points corresponding to the fronts
%         - HtToWdithThresh: the height to width threshold for peaks in the front profiles to estimate non-linearity
% 
% Outputr: - TimeAct_NonLinOnset: time (in s) of nonlinearity onset
%          - NonLinOnsetIndex: Index of nonlinearity onset 
%          - NonLinOrNo: vector mapping each front as non linear or not.
    
    %Determine if each front profile is deemed nonlinear or not based on a provided nonlinearity height-to-width threshold
    % (see relevant user-defined function).
    for i = 1:numel(FrontsCell)
        NonLinOrNo(i) = FrontNonLinTest_PkWdthsAndHts(Xvec, FrontsCell{i}, HtToWdthThresh);
    end

    %Determine the onset of non-linearity as the last index (counted from the end of the nonlinearity flag vector to the beginning) that is tagged as 'linear' under
    % the set height-to-width threshold value.

    NonLinOnsetIndex = NaN; %if every front is non linear, this is set to NaN by default
    
    for i = numel(NonLinOrNo):-1:1
        if NonLinOrNo(i) == 0
            NonLinOnsetIndex = i;
            break
        end
    end

    %get the actual time (in seconds) corresponding to the non-linearity onset
    TimeAct_NonLinOnset = TimeActVec(NonLinOnsetIndex);
end



