function NonLinOrNo = FrontNonLinTest_PkWdthsAndHts(Xvec, FrontsY, WdthToHtRatio_Lt)

%Ritwika VPS, JUne 2025
%This function gets peak widths and heights for each 'peak' (feature) in a given bacterial front to test for non-linearity onset. The peak widths are computed by estimating distance
% between the troughs bracketing each peak (see below for details). The peak heights are computed by taking the max height of the height of the peak from the troughs on either 
% side (see below for details). 
% 
% Then, non-linearity is assessed by checkig if there are features that have height exceeding some factor of the peak width
%
% Inputs: - Xvec, FrontsY: front x and y values
%           - WdthToHtRatioLt: ratio of peak width to height to determine non-linearity onset. A ratio less than this set ratio being present in the front will trigger NonLinOrNo to 1
%                              indicating non-linearity onset
% 
% Outputs: NonLinOrNo: logical determining non linearity onset. 1 for yes, 0 for no.

NonLinOrNo = 0; %set to zero by default

%get peak values and corresponding indices for peaks and troughs (the latter by getting the negative of the curve and then getting peaks).
[pks, locs] = findpeaks(FrontsY);
[pks_neg, locs_neg] = findpeaks(-FrontsY);

PkWdth = NaN*ones(numel(pks),1); %Initialise the peak width vector

% There are going to be, broadly, two possibilities each for the first and last 'feature' in the fronts: a peak or a trough. Therefore, there are, broadly, two cases for these curve
% profiles: i) a peak as the first feature and a trough as the last feature
%           ii) a peak as the first feature AND a peak as the last feature
%           iii) a trough as the first feature and a peak as the last feature
%           iv) a trough as the first feature AND a trough as the last feature
% (Note that findpeaks does not identify a max value at the boundaries of the front (ie, first and last values of Xvec) as a peak, AND that there can't be two consecutive peaks
% without a trough in between, so the analysis above works. In cases (i) and (iii) above, there will be an equal number of peaks and troughs. For case (ii), there will be one 
% more peak than the number of troughs while for case (iv), there will be one less peal than the number of troughs. Utilising these details, we can get the height (from the closest 
% trough) as well as the width of peaks.

%Error checks
if (isempty(pks)) && (numel(pks_neg) > 1)
    error('Peaks cannot be empty when there is more than one trough') 
    %This is because there can be a possibility where there the entire front is a single trough.
end

if (numel(pks) > 1) && (isempty(pks_neg))
    error('Troughs cannot be empyt when thereis more than one peak')
    %This is because there can be a possibility where there the entire front is a single peak.
end

if ~isempty(pks) && ~isempty(pks_neg) %only proceed if there are peaks and troughs in the data
    if numel(pks) == numel(pks_neg) %cases (i) (peak as first and trough as last feature) or (iii) (trough as first and peak as last feature):
    
        if Xvec(locs(1)) < Xvec(locs_neg(1))  %case (i): first peak is before first trough
    
            %   P1                      P3
            %   /\         P2           /\          
            %  /  \        /\          /  \     
            % /    \      /  \        /    \    /  
            %/      \    /    \      /      \  /     
            %        \  /      \    /        \/      
            %         \/        \  /         T3
            %         T1         \/
            %                    T2
    
            PkWdth(2:end) = Xvec(locs_neg(2:end))-Xvec(locs_neg(1:end-1)); %all peak widths except the first one is simple the difference between the x values of 
            %the i-th and (i-1)-th troughs. So, becuase the first trough is *after* the first peak, we can index as above.
            PkWdth(1) = Xvec(locs_neg(1)) - Xvec(1); %the width of the first peak is the differnece between the x value of the first trough and the first x value 
    
            % To get peak heights, we take the height on the left and the right side, and then get the max value of each pair of left and right heights. The left height is simply
            % teh y value of the peak minus the y value of the trough to the left (and similarly for the left height). Note that since the first peak occurs before the first 
            % trough, the first left height is simply the difference between the y value of the first peak and the first y value in the front. Here, we get the y values to be subtracted
            % from the y value of the peak, for both left and right heights.
            PkHt_LeftSubtract = [FrontsY(1) FrontsY(locs_neg(1:end-1))];
            PkHt_RightSubtract = [FrontsY(locs_neg)];
        
            %Get difference from the peak y value and pick the max of left and right heights, for each peak.
            for i = 1:numel(pks)
                PkHt(i) = max(abs(pks(i) - PkHt_LeftSubtract(i)),abs(pks(i) - PkHt_RightSubtract(i)));
            end
    
        elseif Xvec(locs(1)) > Xvec(locs_neg(1)) %case (iii): first trough is before first peak
            
            %                    P2     
            %       P1           /\       P3   
            %       /\          /  \      /\
            %      /  \        /    \    /  \      
            %\    /    \      /      \  /    \
            % \  /      \    /        \/     
            %  \/        \  /         T3       
            %  T1         \/                   
            %             T2 
               
            PkWdth(1:end-1) = Xvec(locs_neg(2:end))-Xvec(locs_neg(1:end-1)); %all peak widths except the last one is simple the difference between the x values of 
            %the (i+1)th and i-th troughs. So, becuase the last peak is *after* the last trough, we can index as above.
            PkWdth(end) = Xvec(end) - Xvec(locs_neg(end)); %the width of the last peak is the differnece between last the first x value and the x value of the last trough 
    
            % Get left and right y values to subtract from the peak y values, to get left and right heights for each peak. Note that since the last peak occurs after the last 
            % trough, the last right height is simply the difference between the y value of the last peak and the last y value in the front. 
            PkHt_LeftSubtract = [FrontsY(locs_neg)];
            PkHt_RightSubtract = [FrontsY(locs_neg(2:end)) FrontsY(end)];
        
            %Get difference from the peak y value and pick the max of left and right heights, for each peak.
            for i = 1:numel(pks)
                PkHt(i) = max(abs(pks(i) - PkHt_LeftSubtract(i)),abs(pks(i) - PkHt_RightSubtract(i)));
            end
        end
        
    elseif numel(pks) == numel(pks_neg) + 1 %case (ii): peak as first AND last feature
    
            %   P1                      P3
            %   /\         P2           /\       P4   
            %  /  \        /\          /  \      /\
            % /    \      /  \        /    \    /  \
            %/      \    /    \      /      \  /    \ 
            %        \  /      \    /        \/      \
            %         \/        \  /         T3
            %         T1         \/
            %                    T2
        
        PkWdth(2:end-1) = Xvec(locs_neg(2:end))-Xvec(locs_neg(1:end-1)); %all peak widths except the first and last ones are simply the difference between the x values of 
        %the i-th and (i-1)-th troughs. 
        PkWdth(1) = Xvec(locs_neg(1)) - Xvec(1); %first peak width is the x value of first trough minues first x value of the front
        PkWdth(end) = Xvec(end) - Xvec(locs_neg(end)); %last peak width is the last x value of the front minus the x value of last trough 
    
        % Get left and right y values to subtract from the peak y values, to get left and right heights for each peak. Note that since there is one peak before and one peak after the 
        % first and last trough respectively, the first left height is the difference between the y value of the first peak and the first y value in the front. The last right height 
        % is the difference between the y value of the last peak and the last y value in the front. 
        PkHt_LeftSubtract = [FrontsY(1) FrontsY(locs_neg)];
        PkHt_RightSubtract = [FrontsY(locs_neg) FrontsY(end)];
    
        %Get difference from the peak y value and pick the max of left and right heights, for each peak.
        for i = 1:numel(pks)
            PkHt(i) = max(abs(pks(i) - PkHt_LeftSubtract(i)),abs(pks(i) - PkHt_RightSubtract(i)));
        end
           
    elseif numel(pks) == numel(pks_neg) - 1 %case (iv): trough as first AND last feature
    
        %                    P2     
        %       P1           /\       P3   
        %       /\          /  \      /\
        %      /  \        /    \    /  \      
        %\    /    \      /      \  /    \    /
        % \  /      \    /        \/      \  /
        %  \/        \  /         T3       \/
        %  T1         \/                   T4
        %             T2
    
        PkWdth = Xvec(locs_neg(2:end))-Xvec(locs_neg(1:end-1));%all peak widths are simply the difference between the x values of the (i+1)-th and i-th troughs. 
        
        % Get left and right y values to subtract from the peak y values, to get left and right heights for each peak. Note that since there is one trough before and one trough after 
        % the first and last peak respectively, the left and right y values to get left and right heights are indexed as below. The left y values is the set of y values for all but 
        % the last trough, and the right y values is the set of y values of all but teh first trough.
        PkHt_LeftSubtract = [FrontsY(locs_neg(1:end-1))];
        PkHt_RightSubtract = [FrontsY(locs_neg(2:end))];
    
        %Get difference from the peak y value and pick the max of left and right heights, for each peak.
        for i = 1:numel(pks)
            PkHt(i) = max(abs(pks(i) - PkHt_LeftSubtract(i)),abs(pks(i) - PkHt_RightSubtract(i)));
        end
    
    else
        error('Front profile with peak-and-trough pattern that has not been accounted for in the code!')
    end
    
    if numel(PkHt) ~= numel(PkWdth)
        error('Peak height and width vectors have different number of elements')
    end
    
    %check to make sure that these are the same size (either all row or column vectors!).
    if ~isequal(size(PkHt), size(PkWdth))
        PkHt = PkHt';
    end

    Inds = 1:numel(PkHt);
    NonLinInds = Inds(PkWdth./PkHt < WdthToHtRatio_Lt);
    if ~isempty(NonLinInds)
        NonLinOrNo = 1;
    end
end

%% Below, some commented out code to plot peak widths, and left and right peak heights for a case(iii) front (this will have to be tweaked accordingly for the other cases):
% 
% %                    P2     
% %       P1           /\       P3   
% %       /\          /  \      /\
% %      /  \        /    \    /  \      
% %\    /    \      /      \  /    \
% % \  /      \    /        \/     
% %  \/        \  /         T3       
% %  T1         \/                   
% %             T2
% %
% figure; hold all
% plot(Xvec,FrontsY); axis tight; grid on %plot front
% plot(Xvec(locs_neg),FrontsY(locs_neg),'Marker','v','MarkerFaceColor','b','MarkerSize',8,'LineStyle','none') %plot peak locations.
% plot(Xvec(locs),FrontsY(locs),'Marker','v','MarkerFaceColor','k','MarkerSize',8,'LineStyle','none') %plot trough locations
% 
% % Plot peak widths (because all peaks are preceded by a trough, we can start each peak width using the trough X index).
% for i = 1:numel(PkWdth)
%     % X value of trough on the left and add peakwidth for the second x value
%     % For y values, we just need the y value of the current trough (cuz y value does not change).
%     plot([Xvec(locs_neg(i)) Xvec(locs_neg(i))+PkWdth(i)], [FrontsY(locs_neg(i)) FrontsY(locs_neg(i))],'k') 
% end
% %
% % Plot left peak heights (because all peaks are preceded by a trough, we can start each peak width using the trough X index).
% for i = 1:numel(PkHt)
%     % Both X values are the same: x value of the trough on the left
%     % Y values start with the y value of the left trough and for the second y value for the height, add the left peak height.
%     plot([Xvec(locs_neg(i)) Xvec(locs_neg(i))], [FrontsY(locs_neg(i)) FrontsY(locs_neg(i))+ abs(pks(i)-PkHt_LeftSubtract(i))],'g','LineWidth',2,'LineStyle','--')
% end
% %
% % Plot left peak heights (because one peak does not have a trough to the right, we have to index accordinlgly). 
% for i = 1:numel(PkHt)-1
%     % For all but the last peak, the x value for the right height is the x value of the (i+1)th trough
%     % The first y value is the y value of the (i+1)th trough, and for the second y value, add the height (for which, we subtract the i-th right subtract value from the y value
%     % of the ith peak, because the (i+1) trough corresponds to the ith peak).
%     plot([Xvec(locs_neg(i+1)) Xvec(locs_neg(i+1))], [FrontsY(locs_neg(i+1)) FrontsY(locs_neg(i+1))+ abs(pks(i)-PkHt_RightSubtract(i))],'b','LineWidth',2,'LineStyle','--')
% end
% % Get the right height for the last peak (without a trough after): x values are the last x values in the front.
% plot([Xvec(end) Xvec(end)], [FrontsY(end) FrontsY(end)+ abs(pks(end)-PkHt_RightSubtract(end))],'b','LineWidth',2,'LineStyle','.-')

end