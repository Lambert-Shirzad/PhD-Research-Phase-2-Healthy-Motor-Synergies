% Phase 1-B: Analysis step 6

% In this step, I will match the synergies identified in step 5 within
% subjects using the similarity value (0.89) found in step 2. 

% 20160704 Written by Navid Lambert-Shirzad

function Match_Synergies_Within_Subjects
    
    DOF = 10;
    ndim_Global = 3;
    similarity_range_DP = 0.91;
    similarity_range_VAF = 89.7;
    similarity_range_DOF_VAF = 85.4;
    IDs = [1:15];%[16:18]
    Hand_Dominance = [1,1,1,0,1,1,1,1,1,0,1,1,1,1,1];%[0,1,1]%1;%
    RangeofMotion = [180, 180, 180, 220, 180, 110, 130, 180, 35, 140]; 

 
    CurrentDirectory = cd;
    CurrentDirectoryUp = strrep(CurrentDirectory,'6 match synergies within subjects',''); 
    DataFolder1 = strcat(CurrentDirectoryUp, '\5 extract ndim_global vectors per subject\');
    DataFolder2 = strcat(CurrentDirectoryUp, '\1 identifying synergies\');

    for SubjCount=1:size(IDs,2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %load the original synergies
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if IDs(SubjCount) < 10
            SubjID = strcat('0', num2str(IDs(SubjCount)));
        else
            SubjID = num2str(IDs(SubjCount));
        end
        load(strcat(DataFolder1,'Y', SubjID, '_Synergies_ndimGlobal' )); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate dot products between all synergy vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DP_All = zeros(ndim_Global, ndim_Global); %DP for dot product
        for i = 1:ndim_Global 
            for j = 1:ndim_Global
                DP_All(i,j) = Synergy_Dom(i,:)*Synergy_NonDom(j,:)'; %Synergies comes from the loaded data
            end
        end      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %match synergy vectors, i.e. which ones have the highest DP 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %matching on the rows
        NumMatched = 0;
        NumDPSimilar(SubjCount, 1) = 0;
        while NumMatched < ndim_Global
            for RowNum = 1:ndim_Global
                %find the max DP in the row
                [MaxDPRow, IndRow] = max(DP_All(RowNum,:));
                %check to see if the max in row is also max in column
                %if so it is a match
                [MaxDPCol, IndCol] = max(DP_All(:,IndRow));
                if MaxDPRow ~= 0 %this row is not previously matched
                    if IndCol == RowNum
                        %this is a match
                        NumMatched = NumMatched+1;
                        DPvalue_Matched_Syn(SubjCount, NumMatched) = MaxDPRow;
                        if DPvalue_Matched_Syn(SubjCount, NumMatched) >= similarity_range_DP
                            NumDPSimilar(SubjCount, 1) = NumDPSimilar(SubjCount, 1) + 1;
                        end
                        matched_ind(SubjCount, NumMatched,:) = [RowNum, IndRow, MaxDPRow];
                        DP_All(RowNum,:) = zeros(size(DP_All(RowNum,:)));
                        DP_All(:,IndRow) = zeros(size(DP_All(:,IndRow)));
                        %RowNum
                        %IndRow                                               
                    end
                    if NumMatched == ndim_Global
                        break
                    end
                end
                
            end    
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %reconstruct motion data of a limb based on synergies of the other
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %load the data (joint motion data)
        load(strcat(DataFolder2,'FullSet_', 'Y', SubjID, '.mat'))
        FullSet(:,2) = FullSet(:,2) + sign(mean(FullSet(:,2)))*mean(FullSet(:,2));
        FullSet(:,3) = FullSet(:,3) + sign(mean(FullSet(:,3)))*mean(FullSet(:,3));
        FullSet(:,4) = FullSet(:,4) + sign(max(FullSet(:,4)))*max(FullSet(:,4));
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12) ]; %RightFullSet = [FullSet(:,2:6) FullSet(:,8:10)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)];%LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:18)]; %t=0s is not included in FullSet (first row is t=1/30s)

        %make sure data is non-negative
        ProcessedRightSide = repmat([90 90 90 90 0 90 0 90 10 70],size(RightFullSet,1),1) + RightFullSet; %[90 90 90 90 0 90 0 90 ] from OpenSim model, abs(lower bound) of each DOF
        ProcessedLeftSide = repmat([90 90 90 90 0 90 0 90 10 70],size(LeftFullSet,1),1) + LeftFullSet;
        
        %normalize the data
% % % % % %         maxRight = max(RightFullSet);
% % % % % %         maxLeft = max(LeftFullSet);       
% % % % % %         for i=1:DOF %assuming right and left side have the same size
% % % % % %             RightFullSet(:,i) = 100*RightFullSet(:,i)/maxRight(1,i);
% % % % % %             LeftFullSet(:,i) = 100*LeftFullSet(:,i)/maxLeft(1,i);
% % % % % %             
% % % % % %         end        
        for i=1:DOF
            ProcessedRightSide(:,i) = 100*ProcessedRightSide(:,i) / RangeofMotion (i);
            ProcessedLeftSide(:,i) = 100*ProcessedLeftSide(:,i) / RangeofMotion (i);
        end
        
        %reconstruct and calculate VAF & DOF_VAF
        if Hand_Dominance(SubjCount) == 1
            Recons_Dom = (ProcessedRightSide / Synergy_NonDom) * Synergy_NonDom;
            Recons_NonDom = (ProcessedLeftSide / Synergy_Dom) * Synergy_Dom;
            VAF_Dom_Recons(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_Dom).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            VAF_NonDom_Recons(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_NonDom).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            DOF_VAF_Dom_Recons(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_Dom).^2,1) ./ sum((ProcessedRightSide).^2,1));
            DOF_VAF_NonDom_Recons(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_NonDom).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            Avg_DOF_VAF_Dom_Recons(SubjCount, 1) = mean(DOF_VAF_Dom_Recons(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons(SubjCount, :));
        else
            Recons_Dom = (ProcessedLeftSide / Synergy_NonDom) * Synergy_NonDom;
            Recons_NonDom = (ProcessedRightSide / Synergy_Dom) * Synergy_Dom;
            VAF_Dom_Recons(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_Dom).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            VAF_NonDom_Recons(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_NonDom).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            DOF_VAF_Dom_Recons(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_Dom).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            DOF_VAF_NonDom_Recons(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_NonDom).^2,1) ./ sum((ProcessedRightSide).^2,1));
            Avg_DOF_VAF_Dom_Recons(SubjCount, 1) = mean(DOF_VAF_Dom_Recons(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons(SubjCount, :));
        end
        VAF_Dom_Original(SubjCount, 1) = VAF_Dom;
        VAF_NonDom_Original(SubjCount, 1) = VAF_NonDom;
        DOF_VAF_Dom_Original(SubjCount, :) = DOF_VAF_Dom;
        DOF_VAF_NonDom_Original(SubjCount, :) = DOF_VAF_NonDom;
        Avg_DOF_VAF_Dom_Original(SubjCount, 1) = mean(DOF_VAF_Dom);
        Avg_DOF_VAF_NonDom_Original(SubjCount, 1) = mean(DOF_VAF_NonDom);    
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Graph the results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure()
    subplot(1,2,1)
    boxplot([VAF_Dom_Original VAF_Dom_Recons VAF_NonDom_Original VAF_NonDom_Recons],'Colors','k')
    ylabel('VAF (%)', 'FontSize',11)
    hold on
    p1 = plot([0.5 4.5],[similarity_range_VAF similarity_range_VAF],'k-');
    Leg1 = 'VAF Similarity Limit = 89.7%';
    legend(p1, Leg1)
    axis([0.5 4.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed' 'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    %title('Dominant Limb')
    
    subplot(1,2,2)
    boxplot([Avg_DOF_VAF_Dom_Original Avg_DOF_VAF_Dom_Recons Avg_DOF_VAF_NonDom_Original Avg_DOF_VAF_NonDom_Recons],'Colors','k')
    ylabel('Average DOF VAF (%)', 'FontSize',11)
    hold on
    p3 = plot([0.5 4.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
    Leg3 = 'DOF VAF Similarity Limit = 85.4%';
    legend(p3, Leg3)
    axis([0.5 4.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed' 'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    title('Within Participant Reconstruction of Joint Motions')
    
    figure()
    
    subplot(3,2,1)
    boxplot([VAF_Dom_Original VAF_Dom_Recons],'Colors','k')
    ylabel('VAF (%)', 'FontSize',11)
    hold on
    p1 = plot([0.5 2.5],[similarity_range_VAF similarity_range_VAF],'k-');
    Leg1 = 'VAF Similarity Limit = 89.7%';
    legend(p1, Leg1)
    axis([0.5 2.5 88.5 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    title('Dominant Limb')
    
    subplot(3,2,2)
    boxplot([VAF_NonDom_Original VAF_NonDom_Recons],'Colors','k')
    hold on
    p2 = plot([0.5 2.5],[similarity_range_VAF similarity_range_VAF],'k-');
    Leg2 = 'VAF Similarity Limit = 89.7%';
    legend(p2, Leg2)
    axis([0.5 2.5 88.5 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    title('Non-Dominant Limb')
    
    subplot(3,2,3)
    boxplot([Avg_DOF_VAF_Dom_Original Avg_DOF_VAF_Dom_Recons],'Colors','k')
    ylabel('Average DOF VAF (%)', 'FontSize',11)
    hold on
    p3 = plot([0.5 2.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
    Leg3 = 'DOF VAF Similarity Limit = 85.4%';
    legend(p3, Leg3)
    axis([0.5 2.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    
    subplot(3,2,4)
    boxplot([Avg_DOF_VAF_NonDom_Original Avg_DOF_VAF_NonDom_Recons],'Colors','k')
    hold on
    p4 = plot([0.5 2.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
    Leg4 = 'DOF VAF Similarity Limit = 85.4%';
    legend(p4, Leg4)
    axis([0.5 2.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Original' 'Reconstructed'} ...
        ,'fontsize',11)
    
    subplot(3,2,5)
    bar([mean(DOF_VAF_Dom_Original); mean(DOF_VAF_Dom_Recons)]', 'grouped')
    axis([0 DOF+1 0 135])
    ylabel('DOF VAF (%) by Muscles', 'FontSize',11)
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',11)
    legend('Original','Reconstructed', 'fontsize',11)
    colormap(gray)
    
    subplot(3,2,6)
    bar([mean(DOF_VAF_NonDom_Original); mean(DOF_VAF_NonDom_Recons)]', 'grouped')
    axis([0 DOF+1 0 135])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',11)
    legend('Original','Reconstructed', 'fontsize',11)
    
    figure()
    boxplot(sort(DPvalue_Matched_Syn, 2, 'descend'),'Colors','k')
    ylabel('Dot Product of Matched Synergy Vectors', 'FontSize',11)
    hold on
    p = plot([0.5 ndim_Global+0.5],[similarity_range_DP similarity_range_DP],'k-');
    Leg = 'Dot Product Similarity Limit = 0.91';
    legend(p, Leg,'fontsize',11)
    axis([0.5 ndim_Global+0.5 0.89 1])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'1st pair' '2nd pair' '3rd pair'} ,'fontsize',11)
    
    for SubjCount = 1:size(IDs,2)
        if IDs(SubjCount) < 10
            SubjID = strcat('0', num2str(IDs(SubjCount)));
        else
            SubjID = num2str(IDs(SubjCount));
        end
        load(strcat(DataFolder1,'Y', SubjID, '_Synergies_ndimGlobal' ));    
        
        temp(:,:) = matched_ind(SubjCount,:,:);
        sorted_matched_ind = sortrows(temp, 3); %sort based on DP values
        figure()
        colormap(gray)
        for i = 1:ndim_Global
            subplot(ndim_Global,1,ndim_Global-i+1)
            bar([Synergy_Dom(sorted_matched_ind(i,1),:); Synergy_NonDom(sorted_matched_ind(i,2),:)]','grouped')
            legend('Dominant Limb', 'Non-dominant Limb')
            text(0.05,0.81,strcat('Dot Product Value = ', num2str(sorted_matched_ind(i,3))), 'fontsize',11)
            ylabel(strcat('Synergy Vector #', num2str(ndim_Global-i+1)))
            axis([0 DOF+1 0 1.05])
            xt = get(gca, 'XTick');
            set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
                    'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',11) 
            if i == ndim_Global
                title(strcat('Participant ID# ', num2str(IDs(SubjCount))))
                savefig(strcat('matched_synergy_vectors_Y', SubjID,'.fig'))
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Save the Symilarity Results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sorted_DPvalue_Matched_Syn = sort(DPvalue_Matched_Syn, 2, 'descend');
    save(strcat('Within_Subjects_Similarity.mat'), ...
            'VAF_Dom_Original', 'VAF_NonDom_Original', 'DOF_VAF_Dom_Original', ...     
            'DOF_VAF_NonDom_Original', 'Avg_DOF_VAF_Dom_Original', 'Avg_DOF_VAF_NonDom_Original', ...
            'VAF_Dom_Recons', 'VAF_NonDom_Recons', 'DOF_VAF_Dom_Recons', 'DOF_VAF_NonDom_Recons', ...
            'Avg_DOF_VAF_Dom_Recons', 'Avg_DOF_VAF_NonDom_Recons', ...
            'sorted_DPvalue_Matched_Syn', 'NumDPSimilar', 'matched_ind');        
        
end %function
    
    
    
  