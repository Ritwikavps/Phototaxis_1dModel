clear; clc

%Ritwika VPS, June 2025
% Script to generate figures for various combinations of scaled bias force (fc_nondim), scaled initial mean conc (Cmean_nondim), and slime decay constant (tau_nondim) to demonstrate
% that they match with phase diagrms in Figs 6 and 9 in the Motility enhancement paper.

graphout = 0; %whether simulations generate figures showing simulation progress or not
Tactual = 600000; %in seconds; total time to simulate; this is the baseline time set in the Motility enhancement paper, Ursell et al, 2013).
%-------------------------------------------------------------------------------------------------------------------------------------------------------------
%Sims without slime decay: in order of points 1, 2, 3, 4, 5, 6 in the phase space diagram in Fig 6 (numbers sourced from TU's simulation code).
%-------------------------------------------------------------------------------------------------------------------------------------------------------------
Cmean_nondim_F6_2exp = [-6, -2, -3, -10, 0, -12]; %exponents for base 2 for Cmean (nondimensionalised)
Cmean_nondim_F6 = 2.^Cmean_nondim_F6_2exp; %get the values by exponentiating

fc_nondim_F6 = [5.12, 5.12, 1.28, 1.28, 0.64, 0.16]; %values
fc_nondim_F6_2exp = log10(fc_nondim_F6)/log10(2);  %Get the exponents by taking log-to-the-base 2, to compare against Fig 6

tau_nondim_F6 = Inf*ones(size(fc_nondim_F6)); %Fig 6 is simulations without slime decay

%Sims for Fig. 6
[Final_C_6A,Final_S_6A] = GetProofOfConceptProfiles_Par(Cmean_nondim_F6,fc_nondim_F6,tau_nondim_F6,graphout, Tactual);


%Plot Concentration profiles
fig6C = figure('Color',[1 1 1]); %Initialise figures to plot slime and bacterial conc

for i = 1:numel(Final_C_6A)
    axes_i_6A_C = subplot(2,3,i,'Parent',fig6C); hold(axes_i_6A_C,'on')
    h = surf(Final_C_6A{i});
    h.EdgeColor = 'none'; axis tight; view(2);
    title(['Fig. 6A' num2str(i) ', (C): Scaled C_{mean} = 2^{' num2str(Cmean_nondim_F6_2exp(i)) ...
                '}, scaled f_c = 2^{' num2str(fc_nondim_F6_2exp(i)) '}, scaled \tau = ' num2str(tau_nondim_F6(i))])
    colorbar
    %I tried to scale the colourbar so that there was a common colourbar for all figures, but that does not represent the various profiles very well.
    hold(axes_i_6A_C,'off'); set(axes_i_6A_C,'YTickLabel',[],'XTickLabel',[],'FontSize',16)
end


%Plot Slime profiles
fig6S = figure('Color',[1 1 1]);  %Initialise figures to plot slime and bacterial conc

for i = 1:numel(Final_S_6A)
    axes_i_6A_S = subplot(2,3,i,'Parent',fig6S); hold(axes_i_6A_S,'on')
    h = surf(Final_S_6A{i});
    h.EdgeColor = 'none'; axis tight; view(2);
    title(['Fig. 6A' num2str(i) ', (S): Scaled C_{mean} = 2^{' num2str(Cmean_nondim_F6_2exp(i)) ...
                '}, scaled f_c = 2^{' num2str(fc_nondim_F6_2exp(i)) '}, scaled \tau = ' num2str(tau_nondim_F6(i))])
    colorbar
    %I tried to scale the colourbar so that there was a common colourbar for all figures, but that does not represent the various profiles very well.
    hold(axes_i_6A_S,'off'); set(axes_i_6A_S,'YTickLabel',[],'XTickLabel',[],'FontSize',16)
end


%-------------------------------------------------------------------------------------------------------------------------------------------------------------
%Sims with slime decay: in order of points 1-6 in the phase space diagram in Fig 9 (numbers sourced from TU's simulation code).
%-------------------------------------------------------------------------------------------------------------------------------------------------------------
Cmean_nondim_F9_2exp = [0, 0, -1, -5, 3, -10]; 
Cmean_nondim_F9 = 2.^Cmean_nondim_F9_2exp;

fc_nondim_F9 = 4*ones(size(Cmean_nondim_F9)); %source: Motility enhancement paper, pg 10

tau_nondim_F9_2exp = [1, 0, 0, 2, -1, -2];
tau_nondim_F9 = 2.^tau_nondim_F9_2exp;

%Sims for Fig. 6
[Final_C_9B,Final_S_9B] = GetProofOfConceptProfiles_Par(Cmean_nondim_F9,fc_nondim_F9,tau_nondim_F9,graphout, Tactual);


%Plot Concentration profiles
fig9C = figure('Color',[1 1 1]); 

for i = 1:numel(Final_C_9B)
    axes_i_9B_C = subplot(2,3,i,'Parent',fig9C); hold(axes_i_9B_C,'on')
    h = surf(Final_C_9B{i});
    h.EdgeColor = 'none'; axis tight; view(2);
    title(['Fig. 9B' num2str(i) ', (C): Scaled C_{mean} = 2^{' num2str(Cmean_nondim_F9_2exp(i)) ...
                '}, scaled \tau = 2^{' num2str(tau_nondim_F9_2exp(i)) '}, scaled f_c = ' num2str(fc_nondim_F9(i))])
    colorbar
    hold(axes_i_9B_C,'off'); set(axes_i_9B_C,'YTickLabel',[],'XTickLabel',[],'FontSize',16)
end


%Plot Slime profiles
fig9S = figure('Color',[1 1 1]);  %Initialise figures to plot slime and bacterial conc

for i = 1:numel(Final_S_9B)
    axes_i_9B_S = subplot(2,3,i,'Parent',fig9S); hold(axes_i_9B_S,'on')
    h = surf(Final_S_9B{i});
    h.EdgeColor = 'none'; axis tight; view(2);
    title(['Fig. 9B' num2str(i) ', (S): Scaled C_{mean} = 2^{' num2str(Cmean_nondim_F9_2exp(i)) ...
                '}, scaled \tau = 2^{' num2str(tau_nondim_F9_2exp(i)) '}, scaled f_c = ' num2str(fc_nondim_F9(i))])
    colorbar
    hold(axes_i_9B_S,'off'); set(axes_i_9B_S,'YTickLabel',[],'XTickLabel',[],'FontSize',16)
end


%-------------------------------------------------------------------------------------------------------------------------------------------------------------
%Parallel for loop fn
%-------------------------------------------------------------------------------------------------------------------------------------------------------------
% This function parallelises computing the final bacterial and slime concentration profiles for the different cases demonstrated in the phase diagrams in the Motility
% enhancement paper.
% 
% Inputs: CmeanVec, FcVec, TauVec: Vectors with values for scaled initial mean concetration, bias force, and slime decay time constant, in the order of 1-6 in the
%                                  phase diagram for both figs (all vectors the same size)
%         graphoput: graph out logical to pass to the simulation function
%         Tactual: total actual time to simulate, set as an input because non-linearity onset times differ for different parameter combos.
% 
% Outputs: Final_C, Final_S: cell arrays with final bacterical (C) and slime (S) concentration profiles for each combo of parameters for the different cases demo-d in the phase
%                            space diagrams.
function [Final_C,Final_S] = GetProofOfConceptProfiles_Par(CmeanVec,FcVec,TauVec,graphout, Tactual)

    N = numel(CmeanVec); %get number of elements

    %Error check
    if (numel(CmeanVec) ~= numel(FcVec)) || (numel(CmeanVec) ~= numel(TauVec))
        error('Input vectors should have the same length')
    end

    Final_C = cell(1,N); Final_S = cell(1,N); %Initialise outputs

    p = parpool(N); %open parallel pool

    parfor i = 1:numel(CmeanVec) %do the simulations, parallelised
        [~,~,~,~,Final_C{i}, Final_S{i}, ~, ~] = GetPhotoFronts_w_SlimeDecay(CmeanVec(i),FcVec(i),TauVec(i),graphout, Tactual, 1); %1 for the factor to divide 
        % dt (dimensionless time increment) so all sims carried out have roughly similar number of simulation time steps; see GetPhotoFronts_w_SlimeDecay.m.
    end

    delete(p) %delete parallel pool
end

