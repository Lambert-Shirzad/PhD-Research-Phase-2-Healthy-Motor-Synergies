% Phase 1-B: Analysis step 4
% ndim_Global = 4 from Step 3 of the analysis 
% Factorize motion data of healthy participants into ndim_Global =1:DOF
% synergies for each of the dominant and non-dominant sides of the body 
% using NNMF and save the goodness of fit stats to finalize that 
% ndim_Global = 4 is valid. 

% The function batch-processes the data of participants.
% SubjectIDs: ID of participants that this step of analysis will be run for
% try [1:12] to run participants 1-12
% HandDominance: Which hand is dominant for each subject in SubjectIDs, 1
% for right hand and 0 for left hand. Example: [1,1,0]

%  Establish_ndim_Globally_2([1,2,5,6,9,10,11,12,14,15,17,18,19,20,21],ones(1,15))

% All_DOF_Synergies_Stats_Dom.mat: for dominant hand of all subjects, has
%                                 VAF, DOF_VAF, Delta_VAF for ndim=1:DOF
% All_DOF_Synergies_Stats_Dom.mat: for non-dominant hand of all subjects, 
%                                  has VAF, DOF_VAF, Delta_VAF for ndim=1:DOF


% 20160628 Written by Navid Lambert-Shirzad

function Establish_ndim_Globally_2(SubjectIDs, HandDominance)
    
    ndim_Global = 4;
    DOF = 8;
    %% Train the 1:DOF Synergies for Each Subject    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end
      
        RightVAF = zeros(DOF,1); LeftVAF = zeros(DOF,1); %VAF for each number of synergies
        DeltaRightVAF = zeros(DOF,1); DeltaLeftVAF = zeros(DOF,1); %how much VAF changes when a synergy vector is added
        RightDOF_VAF = zeros(DOF,DOF); LeftDOF_VAF = zeros(DOF,DOF); %how VAF of each DOF changes as more synergies are added for each fold. first DOF is the synergies included, the second is the DOF being observed, third is the fold number
        
        tic
        %load the data (EMG data)
        CurrentDirectory = cd;
        CurrentDirectoryUp = strrep(CurrentDirectory,'4 establish ndim globally 2',''); 
        OrigSynergiesFolder = strcat(CurrentDirectoryUp, '\1 identifying synergies\');
        load(strcat(OrigSynergiesFolder,'Processed_Subj_', SubjID, '_Left.mat'))
        load(strcat(OrigSynergiesFolder,'Processed_Subj_', SubjID, '_Right.mat'))
        RightFullSet = ProcessedRightSide(:,2:9); %1st column is time
        LeftFullSet = ProcessedLeftSide(:,2:9); 

        %make sure data is non-negative: EMG data has been preprocessed
        %to make sure it is non-negative (in fact normalized to each
        %channel's max and saved as a percentage). So, yes, check!

        %train ndim=1:DOF synergies and record all the data so you can generate plots
        GoodTrainR = 0; GoodTrainL=0; 
        numSynergy = 1;
        subjectcounter
        while numSynergy < DOF
            numSynergy
            %perform NNMF on data
            true = 0; count = 1; Right = 0
            while true == 0 && count < 10
                [ScoresRightTrTemp, SynergiesRightTrAll] = nnmf(RightFullSet, numSynergy); 
                if rank(SynergiesRightTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
                end
                count = count + 1;
                if count == 10
                    %[ScoresRightTrTemp, SynergiesRightTrAll] = nnmf(RightFullSet, numSynergy - 1);
                    rand;
                    count = 1;
                    numSynergy
                    Right = 1
                end
            end
            true = 0; count = 1; Left = 0
            while true == 0 && count < 10
                [ScoresLeftTrTemp, SynergiesLeftTrAll] = nnmf(LeftFullSet, numSynergy); 
                if rank(SynergiesLeftTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
                end
                count = count + 1;
                if count == 10
                    %[ScoresLeftTrTemp, SynergiesLeftTrAll] = nnmf(LeftFullSet, numSynergy - 1);
                    rand;
                    count = 1;
                    numSynergy
                    Left = 1
                end
            end

            RightApprox = ScoresRightTrTemp * SynergiesRightTrAll;
            RightVAF(numSynergy,1) = 100*(1 - (sum(sum((RightFullSet - RightApprox).^2,2),1)) / (sum(sum((RightFullSet).^2,2),1))); %1-SSE/SST
            LeftApprox = ScoresLeftTrTemp * SynergiesLeftTrAll;
            LeftVAF(numSynergy,1) = 100*(1 - (sum(sum((LeftFullSet - LeftApprox).^2,2),1)) / (sum(sum((LeftFullSet).^2,2),1))); %1-SSE/SST
            if numSynergy ~= 1
                DeltaRightVAF(numSynergy,1)=RightVAF(numSynergy,1)-RightVAF(numSynergy-1,1);
                DeltaLeftVAF(numSynergy,1)=LeftVAF(numSynergy,1)-LeftVAF(numSynergy-1,1);
            else
                DeltaRightVAF(numSynergy,1)=RightVAF(numSynergy,1);
                DeltaLeftVAF(numSynergy,1)=LeftVAF(numSynergy,1);
            end
            RightDOF_VAF(numSynergy,:) = 100*(1 - sum((RightFullSet - RightApprox).^2,1) ./ sum((RightFullSet).^2,1));
            LeftDOF_VAF(numSynergy,:) = 100*(1 - sum((LeftFullSet - LeftApprox).^2,1) ./ sum((LeftFullSet).^2,1));

            numSynergy = numSynergy+1;
        end          

        timeElapsed = toc;
        
        if HandDominance(subjectcounter) == 1
            SubjDomHand(subjectcounter,1) = 'R';
            SubjNonDomHand(subjectcounter,1) = 'L';
            VAF_Dom_All(1:DOF,subjectcounter) = RightVAF(1:DOF,1);
            DeltaVAF_Dom_All(1:DOF,subjectcounter) = DeltaRightVAF(1:DOF,1);
            DOF_VAF_Dom_All(1:DOF,1:DOF,subjectcounter) = RightDOF_VAF(1:DOF,1:DOF);
            VAF_NonDom_All(1:DOF,subjectcounter) = LeftVAF(1:DOF,1);
            DeltaVAF_NonDom_All(1:DOF,subjectcounter) = DeltaLeftVAF(1:DOF,1);
            DOF_VAF_NonDom_All(1:DOF,1:DOF,subjectcounter) = LeftDOF_VAF(1:DOF,1:DOF);
        else
            SubjDomHand(subjectcounter,1) = 'L';
            SubjNonDomHand(subjectcounter,1) = 'R';
            VAF_NonDom_All(1:DOF,subjectcounter) = RightVAF(1:DOF,1);
            DeltaVAF_NonDom_All(1:DOF,subjectcounter) = DeltaRightVAF(1:DOF,1);
            DOF_VAF_NonDom_All(1:DOF,1:DOF,subjectcounter) = RightDOF_VAF(1:DOF,1:DOF);
            VAF_Dom_All(1:DOF,subjectcounter) = LeftVAF(1:DOF,1);
            DeltaVAF_Dom_All(1:DOF,subjectcounter) = DeltaLeftVAF(1:DOF,1);
            DOF_VAF_Dom_All(1:DOF,1:DOF,subjectcounter) = LeftDOF_VAF(1:DOF,1:DOF);
        end 
    end
    VAF_NonDom_All(DOF,:) = 100;  
    DOF_VAF_NonDom_All(DOF,:,:) = 100;
    VAF_Dom_All(DOF,:) = 100;
    DOF_VAF_Dom_All(DOF,:,:) = 100;
    
    %% Save the Goodness of Fit Stats and Plot the Results
    save('All_DOF_Synergies_Stats_Dom.mat', ...
                'VAF_Dom_All', 'DeltaVAF_Dom_All', 'DOF_VAF_Dom_All', 'SubjDomHand');
    save('All_DOF_Synergies_Stats_NonDom.mat', ...
                'VAF_NonDom_All', 'DeltaVAF_NonDom_All', 'DOF_VAF_NonDom_All', 'SubjNonDomHand'); 
    
    figure()
    subplot(3,2,1)
    boxplot(VAF_Dom_All','Colors','k')
    hold on
    plot(median(VAF_Dom_All'),'k-')
    ylabel('VAF (%)', 'FontSize',12)
    axis([0.5 DOF+0.5 70 105])
    title('Dominant Limb','FontSize',12)
    
    subplot(3,2,2)
    boxplot(VAF_NonDom_All','Colors','k')
    hold on
    plot(median(VAF_NonDom_All'),'k-')
    ylabel('VAF (%)', 'FontSize',12)
    axis([0.5 DOF+0.5 70 105])
    title('Non-Dominant Limb','FontSize',12)
    
    subplot(3,2,3)
    boxplot(DeltaVAF_Dom_All','Colors','k')
    hold on
    plot([2:DOF], median(DeltaVAF_Dom_All(2:DOF,:)'),'k-')
    ylabel('Delta VAF (%)', 'FontSize',12)
    axis([0.5 DOF+0.5 0 20])
    
    subplot(3,2,4)
    boxplot(DeltaVAF_NonDom_All','Colors','k')
    hold on
    plot([2:DOF], median(DeltaVAF_NonDom_All(2:DOF,:)'),'k-')
    ylabel('Delta VAF (%)', 'FontSize',12)
    axis([0.5 DOF+0.5 0 20])
    
    subplot(3,2,5)
    hold on
    for i = 1:DOF
        plot(mean([DOF_VAF_Dom_All(:,i,1) DOF_VAF_Dom_All(:,i,2) ... %DOF i
            DOF_VAF_Dom_All(:,i,3) DOF_VAF_Dom_All(:,i,4) ...
            DOF_VAF_Dom_All(:,i,5) DOF_VAF_Dom_All(:,i,6) ...
            DOF_VAF_Dom_All(:,i,7) DOF_VAF_Dom_All(:,i,8) ...
            DOF_VAF_Dom_All(:,i,9) DOF_VAF_Dom_All(:,i,10) ...
            DOF_VAF_Dom_All(:,i,11) DOF_VAF_Dom_All(:,i,12) ...
            DOF_VAF_Dom_All(:,i,13) DOF_VAF_Dom_All(:,i,14) ...
            DOF_VAF_Dom_All(:,i,15)], 2))
    end
    xlabel('Number of Synergies', 'FontSize',12)
    ylabel('DOF VAF (%)', 'FontSize',12)
    %axis([0.5 8.5 0 105])
    axis([0.5 DOF+0.5 60 105])
    MuscleVector={'DelAnt' 'DelMed' 'DelPos' ...
        'Biceps' 'TriLong' 'TriLat' 'Brachi' 'PectMaj'};
	legend(MuscleVector,'location','southeast')
    
    subplot(3,2,6)
    hold on
    for i = 1:DOF
        plot(mean([DOF_VAF_NonDom_All(:,i,1) DOF_VAF_NonDom_All(:,i,2) ... %DOF i
            DOF_VAF_NonDom_All(:,i,3) DOF_VAF_NonDom_All(:,i,4) ...
            DOF_VAF_NonDom_All(:,i,5) DOF_VAF_NonDom_All(:,i,6) ...
            DOF_VAF_NonDom_All(:,i,7) DOF_VAF_NonDom_All(:,i,8) ...
            DOF_VAF_NonDom_All(:,i,9) DOF_VAF_NonDom_All(:,i,10) ...
            DOF_VAF_NonDom_All(:,i,11) DOF_VAF_NonDom_All(:,i,12) ...
            DOF_VAF_NonDom_All(:,i,13) DOF_VAF_NonDom_All(:,i,14) ...
            DOF_VAF_NonDom_All(:,i,15)], 2))
    end
    xlabel('Number of Synergies', 'FontSize',12)
    ylabel('DOF VAF (%)', 'FontSize',12)
    axis([0.5 DOF+0.5 60 105])
    legend(MuscleVector,'location','southeast')
    
    %% Train ndim_Global = 4 synergy vectors and save them for comparison
    
    
%         if HandDominance(subjectcounter) == 1
%             SubjDomHand = 'R';
%             save(strcat('Y', SubjID, '_Dom_Synergies_Globalndim.mat'), ...
%                 'VAF_R', 'DeltaVAF_R', 'DOF_VAF_R', 'Synergies_R', ...
%                 'Scores_R', 'ndim_R', 'AvgReconsErr_R', 'RegCoeff_R', ...
%                 'nCommonAll', 'SubjDomHand', 'timeElapsed');
% 
%             SubjNonDomHand = 'L';
%             save(strcat('Y', SubjID, '_NonDom_Synergies_Globalndim.mat'), ...
%                 'VAF_L', 'DeltaVAF_L', 'DOF_VAF_L', 'Synergies_L', ...
%                 'Scores_L', 'ndim_L', 'AvgReconsErr_L', 'RegCoeff_L', ...
%                 'nCommonAll', 'SubjNonDomHand', 'timeElapsed');
%         else
%             SubjDomHand = 'L';
%             save(strcat('Y', SubjID, '_Dom_Synergies_Globalndim.mat'), ...
%                 'VAF_L', 'DeltaVAF_L', 'DOF_VAF_L', 'Synergies_L', ...
%                 'Scores_L', 'ndim_L', 'AvgReconsErr_L', 'RegCoeff_L', ...
%                 'nCommonAll', 'SubjDomHand', 'timeElapsed');
% 
%             SubjNonDomHand = 'R';
%             save(strcat('Y', SubjID, '_NonDom_Synergies_Globalndim.mat'), ...
%                 'VAF_R', 'DeltaVAF_R', 'DOF_VAF_R', 'Synergies_R', ...
%                 'Scores_R', 'ndim_R', 'AvgReconsErr_R', 'RegCoeff_R', ...
%                 'nCommonAll', 'SubjNonDomHand', 'timeElapsed');
% 
%         end 
    
    
end