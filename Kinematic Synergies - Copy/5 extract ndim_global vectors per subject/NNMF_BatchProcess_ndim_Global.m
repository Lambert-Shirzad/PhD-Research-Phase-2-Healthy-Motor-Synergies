% Phase 1-B: Analysis step 5
% ndim_Global = 4 from the previous two steps of the analysis 
% Factorize motion data of healthy participants into ndim_Global =4 
% synergies for each of the dominant and non-dominant sides of the body 
% using NNMF and save the synergies for the next steps of analysis.

% The function batch-processes the data of participants.
% SubjectIDs: ID of participants that this step of analysis will be run for
% try [1:12] to run participants 1-12
% HandDominance: Which hand is dominant for each subject in SubjectIDs, 1
% for right hand and 0 for left hand. Example: [1,1,0]

% NNMF_BatchProcess_ndim_Global([1,2,5,6,9,10,11,12,14,15,17,18,19,20,21],ones(1,15))

% OUTPUTS:

% Y##_Synergies_ndimGlobal.mat: for both hands of subject##, has
%            VAF, DOF_VAF, Delta_VAF, Synergies, Synergy_Scores
%            for ndim = ndim_Global


% 20160704 Written by Navid Lambert-Shirzad

function NNMF_BatchProcess_ndim_Global(SubjectIDs, HandDominance)
    
    ndim_Global = 3;
    DOF = 10;
     RangeofMotion = [180, 180, 180, 220, 180, 110, 130, 180, 35, 140]; 
    %% Train the 1:DOF Synergies for Each Subject    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end
      
        RightVAF = zeros(DOF,1); LeftVAF = zeros(DOF,1); %VAF for each number of synergies
        DeltaRightVAF = zeros(DOF,1); DeltaLeftVAF = zeros(DOF,1); %how much VAF changes when a synergy vector is added
        RightDOF_VAF = zeros(ndim_Global,DOF); LeftDOF_VAF = zeros(ndim_Global,DOF); %how VAF of each DOF changes as more synergies are added for each fold. first DOF is the synergies included, the second is the DOF being observed, third is the fold number
        
        %load the data 
        CurrentDirectory = cd;
        CurrentDirectoryUp = strrep(CurrentDirectory,'5 extract ndim_global vectors per subject',''); 
        OrigSynergiesFolder = strcat(CurrentDirectoryUp, '\1 identifying synergies\');
        %load the data (joint motion data)
        load(strcat(OrigSynergiesFolder,'FullSet_', 'Y', SubjID, '.mat'))
        FullSet(:,2) = FullSet(:,2) + sign(mean(FullSet(:,2)))*mean(FullSet(:,2));
        FullSet(:,3) = FullSet(:,3) + sign(mean(FullSet(:,3)))*mean(FullSet(:,3));
        FullSet(:,4) = FullSet(:,4) + sign(max(FullSet(:,4)))*max(FullSet(:,4));
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12) ]; %RightFullSet = [FullSet(:,2:6) FullSet(:,8:10)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)];%LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:18)]; %t=0s is not included in FullSet (first row is t=1/30s)

        %make sure data is non-negative
        RightFullSet = repmat([90 90 90 90 0 90 0 90 10 70],size(RightFullSet,1),1) + RightFullSet; %[90 90 90 90 0 90 0 90 ] from OpenSim model, abs(lower bound) of each DOF
        LeftFullSet = repmat([90 90 90 90 0 90 0 90 10 70],size(LeftFullSet,1),1) + LeftFullSet;
        
        %normalize the data
% % % % % %         maxRight = max(RightFullSet);
% % % % % %         maxLeft = max(LeftFullSet);       
% % % % % %         for i=1:DOF %assuming right and left side have the same size
% % % % % %             RightFullSet(:,i) = 100*RightFullSet(:,i)/maxRight(1,i);
% % % % % %             LeftFullSet(:,i) = 100*LeftFullSet(:,i)/maxLeft(1,i);
% % % % % %             
% % % % % %         end        
        for i=1:DOF
            RightFullSet(:,i) = 100*RightFullSet(:,i) / RangeofMotion (i);
            LeftFullSet(:,i) = 100*LeftFullSet(:,i) / RangeofMotion (i);
        end


        %make sure data is non-negative: EMG data has been preprocessed
        %to make sure it is non-negative (in fact normalized to each
        %channel's max and saved as a percentage). So, yes, check!

        %train ndim=ndim_Global synergies and record all the data 
        GoodTrainR = 0; GoodTrainL=0; 
        numSynergy = 1; %ndim_Global;
        while numSynergy < ndim_Global+1
            %numSynergy
            %perform NNMF on data
            true = 0; 
            while true == 0
                [ScoresLeftTrTemp, SynergiesLeftTrAll] = nnmf(LeftFullSet, numSynergy); 
                if rank(SynergiesLeftTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
                end
            end
            true = 0; 
            while true == 0
                [ScoresRightTrTemp, SynergiesRightTrAll] = nnmf(RightFullSet, numSynergy, 'h0', SynergiesLeftTrAll); %+ 0.1*rand(size(SynergiesLeftTrAll))
                if rank(SynergiesRightTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
                end
            end
             true = 0; 
            while true == 0
                [ScoresLeftTrTemp, SynergiesLeftTrAll] = nnmf(LeftFullSet, numSynergy, 'h0', SynergiesRightTrAll); %+ 0.1*rand(size(SynergiesRightTrAll))
                if rank(SynergiesLeftTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
                end
            end
            true = 0; 
            while true == 0
                [ScoresRightTrTemp, SynergiesRightTrAll] = nnmf(RightFullSet, numSynergy, 'h0', SynergiesLeftTrAll); %+ 0.1*rand(size(SynergiesLeftTrAll))
                if rank(SynergiesRightTrAll) == numSynergy
                    true = 1; %not underfitting or stuck in local minima
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
            
            %check if the goodness of fit criteria is upheld
            if numSynergy == ndim_Global
                if RightVAF(numSynergy,1)>90 & DeltaRightVAF(numSynergy,1)<5 & RightDOF_VAF(numSynergy,:)>65
                    ndim_R = ndim_Global;
                    VAF_R = RightVAF(numSynergy,1);
                    DeltaVAF_R = DeltaRightVAF(numSynergy,1);
                    DOF_VAF_R = RightDOF_VAF(numSynergy,:);
                    GoodTrainR = 1;
                    Synergies_R = SynergiesRightTrAll;
                    Scores_R = ScoresRightTrTemp;
                end
                if LeftVAF(numSynergy,1)>90 & DeltaLeftVAF(numSynergy,1)<5 & LeftDOF_VAF(numSynergy,:)>65
                    ndim_L = ndim_Global;
                    VAF_L = LeftVAF(numSynergy,1);
                    DeltaVAF_L = DeltaLeftVAF(numSynergy,1);
                    DOF_VAF_L = LeftDOF_VAF(numSynergy,:);
                    GoodTrainL = 1;
                    Synergies_L = SynergiesLeftTrAll;
                    Scores_L = ScoresLeftTrTemp;
                end
                if GoodTrainR == 0 || GoodTrainL == 0
                    numSynergy = numSynergy - 1; 
                else
                    numSynergy = ndim_Global + 1; %to terminate while loop
                end
            else
                numSynergy = numSynergy + 1;
                
            end

        end   
        
   
        if HandDominance(subjectcounter) == 1
            Hand_Dom = 'R';
            Hand_NonDom = 'L';
            
            VAF_Dom = VAF_R;
            DeltaVAF_Dom = DeltaVAF_R;
            DOF_VAF_Dom = DOF_VAF_R;
            Synergy_Dom = Synergies_R;
            Scores_Dom = Scores_R;
            
            VAF_NonDom = VAF_L;
            DeltaVAF_NonDom = DeltaVAF_L;
            DOF_VAF_NonDom = DOF_VAF_L;
            Synergy_NonDom = Synergies_L;
            Scores_NonDom = Scores_L;
        else
            Hand_Dom = 'L';
            Hand_NonDom = 'R';
            
            VAF_NonDom = VAF_R;
            DeltaVAF_NonDom = DeltaVAF_R;
            DOF_VAF_NonDom = DOF_VAF_R;
            Synergy_NonDom = Synergies_R;
            Scores_NonDom = Scores_R;
            
            VAF_Dom = VAF_L;
            DeltaVAF_Dom = DeltaVAF_L;
            DOF_VAF_Dom = DOF_VAF_L;
            Synergy_Dom = Synergies_L;
            Scores_Dom = Scores_L;
        end 
        save(strcat('Y', SubjID, '_Synergies_ndimGlobal.mat'), ...
            'Hand_Dom', 'Hand_NonDom', ...     
            'VAF_NonDom', 'DeltaVAF_NonDom', 'DOF_VAF_NonDom', 'Synergy_NonDom', ...
            'Scores_NonDom', 'VAF_Dom', 'DeltaVAF_Dom', 'DOF_VAF_Dom', ...
            'Synergy_Dom', 'Scores_Dom')
    end 
    
end