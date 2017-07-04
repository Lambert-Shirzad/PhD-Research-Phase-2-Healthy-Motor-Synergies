% Phase 1-B: Analysis step 2

% In step 1, I factorized the muscle/kinematic data of healthy participants
% and generated a set of muscle synergy vectors and numerical values for
% goodness of fit. In this step, I will use the identified synergy vectors.
% 
% A) generate a set of 100 random synergy vectors
% B) generate dot product population distribution (N=100*100=10,000) as
% well as VAF & DOF_VAF population distributions
% C) create a normal-distrubation graph. Mark the 95th percentile.
% D) the 95th percentile is the similarity limit
%
% Three sets of random synergy vector sets can be considered in part A. 
% A.1 purely randomly generated vectors
% A.2 randomly sample from original vectors' indicies and combine them to 
% build new vectors
% A.3 randomly sample from index i of the original vectors and assign the
% outcome to index i of a new vector
% Original vectors refer to the synergy structures identified in step 1 of
% the analysis. 
%
% I will choose the highest 95th percentile between A.1-3 as the similarity
% level. These numbers will be save in a .mat file (SimilarityLimit.mat).
% The distributation graphs will also be saved. 

% 20160620 Written by Navid Lambert-Shirzad

function Synergy_Similarity_Range
    
    DOF = 10;
    IDs = [1:15];
    %% load the original synergies into two matrices, left and right. 
    CurrentDirectory = cd;
    CurrentDirectoryUp = strrep(CurrentDirectory,'2 dot product similarity limit',''); %we are in PCA, data is saved up one level.
    OrigSynergiesFolder = strcat(CurrentDirectoryUp, '1 identifying synergies\');
    
    for i=1:size(IDs,2)
        if IDs(i) < 10
            SubjID = strcat('0', num2str(IDs(i)));
        else
            SubjID = num2str(IDs(i));
        end
        
        %load the data (joint motion data)
        load(strcat(OrigSynergiesFolder,'FullSet_', 'Y', SubjID, '.mat'))
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12) ]; %RightFullSet = [FullSet(:,2:6) FullSet(:,8:10)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)];%LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:18)]; %t=0s is not included in FullSet (first row is t=1/30s)

        %make sure data is non-negative
        RightFullSet = repmat([90 90 90 90 0 90 0 90 10 70],size(RightFullSet,1),1) + RightFullSet; %[90 90 90 90 0 90 0 90 ] from OpenSim model, abs(lower bound) of each DOF
        LeftFullSet = repmat([90 90 90 90 0 90 0 90 10 70],size(LeftFullSet,1),1) + LeftFullSet;

        
        load(strcat(OrigSynergiesFolder,'Y', SubjID, '_Dom_Synergies.mat' )); 
        load(strcat(OrigSynergiesFolder,'Y', SubjID, '_NonDom_Synergies.mat' )); 
        if i == 1
            OrigSynL = Synergies_L; %from mat files loaded above
            OrigSynR = Synergies_R; %from mat files loaded above
        else
            OrigSynL = [OrigSynL; Synergies_L]; 
            OrigSynR = [OrigSynR; Synergies_R];
        end
        MotionData(:,:,i) = RightFullSet(1:8400,:); 
    end
    
    %% A) generate a set of numVectors random synergy vectors
    numVectors = 250;
% %     % A.1 purely randomly generated vectors
% %     RandSetA1 = zeros(numVectors,DOF);
% %     for i=1:numVectors
% %         R1 = rand(1,DOF);
% %         RandSetA1(i,:) = R1/norm(R1);    
% %     end
    % A.2 randomly sample from original vectors' indicies and combine them to 
    % build new vectors
    OrigSynAll = [OrigSynR; OrigSynL];
    RandSetA2 = zeros(numVectors,DOF);
    for i = 1:numVectors
        R2 = zeros(1,DOF);
        for j = 1:DOF
            R2(1,j) = OrigSynAll(randi([1, size(OrigSynAll,1)]), randi([1, size(OrigSynAll,2)]));
        end
        RandSetA2(i,:) = R2/norm(R2); 
    end
% %     % A.3 randomly sample from index i of the original vectors and assign the
% %     % outcome to index i of a new vector
% %     RandSetA3 = zeros(numVectors,DOF);
% %     for i = 1:numVectors
% %         R3 = zeros(1,DOF);
% %         for j = 1:DOF
% %             R3(1,j) = OrigSynAll(randi([1, size(OrigSynAll,1)]), j);
% %         end
% %         RandSetA3(i,:) = R3/norm(R3); 
% %     end
    
    %% B) generate population distribution (N=numVectors*numVectors)
    populationN = numVectors*numVectors;
% %     DistA1 = zeros(populationN, 1);
    DistA2 = zeros(populationN, 1);
% %     DistA3 = zeros(populationN, 1);
    
    for i = 1:numVectors
        for j = 1:numVectors
% %             DistA1( (i-1)*numVectors+j ,1) = RandSetA1(i,:)*RandSetA1(j,:)';
            DistA2( (i-1)*numVectors+j ,1) = RandSetA2(i,:)*RandSetA2(j,:)';
% %             DistA3( (i-1)*numVectors+j ,1) = RandSetA3(i,:)*RandSetA3(j,:)';
        end
    end
    
% %     DistributionA1 = zeros(populationN-numVectors, 1);
    DistributionA2 = zeros(populationN-numVectors, 1);
% %     DistributionA3 = zeros(populationN-numVectors, 1);
    unwanteds = [1:numVectors+1:populationN]; %we have a lot of unwanted zeros at these indices
    count = 1;
    for i = 1:populationN
        if i==unwanteds(count)
            count = count+1;
        else
% %             DistributionA1(i-count+1)= DistA1(i);
            DistributionA2(i-count+1)= DistA2(i);
% %             DistributionA3(i-count+1)= DistA3(i);
        end
    end
    
    for i = 1:size(IDs,2)
        for j =1:200 %for each subject try 200 random set of synergies
            P = randperm(numVectors,4);
            Synergy_random = [RandSetA2(P(1),:); RandSetA2(P(2),:); RandSetA2(P(3),:); RandSetA2(P(4),:)];
            Recons_MotionData = (MotionData(:,:,i) / Synergy_random) * Synergy_random;
            VAF_Distribution((i-1)*200+j ,1) = 100*(1 - (sum(sum((MotionData(:,:,i) - Recons_MotionData).^2,2),1)) / (sum(sum((MotionData(:,:,i)).^2,2),1)));           
            DOF_VAF_Distribution((i-1)*200+j ,1:DOF)= 100*(1 - sum((MotionData(:,:,i) - Recons_MotionData).^2,1) ./ sum((MotionData(:,:,i)).^2,1));
        end
    end
    
    
    %% C) create a normal-distrubation graph. Mark the 95th percentile.
%     percentileA1 = prctile(DistributionA1,95)
    percentileA2 = prctile(DistributionA2,95)
%     percentileA3 = prctile(DistributionA3,95)
    
    nbins = 50;
% %     figure()
% %     h1 = histogram(DistributionA1, nbins);
% %     %h1 = histfit(DistributionA1, nbins);
% %     y1=get(gca,'ylim');
% %     hold on
% %     p1 = plot([percentileA1 percentileA1],y1,'LineWidth',2.5);
% %     LegPerc1 = strcat('95th percentile=', num2str(percentileA1));
% %     legend(p1, LegPerc1)
    
    figure()
    h2 = histogram(DistributionA2, nbins);
    %axis([0 1 0 4000])
    y2=get(gca,'ylim');
    hold on
    p2 = plot([percentileA2 percentileA2],y2,'LineWidth',2.5);
    LegPerc2 = strcat('95th percentile=', num2str(floor(100*percentileA2)/100));
    legend(p2, LegPerc2)
    xlabel('Dot Product Values')
    ylabel('Frequency')
    
% %     figure()
% %     h3 = histogram(DistributionA3, nbins);
% %     %h3 = histfit(DistributionA3, nbins);
% %     y3=get(gca,'ylim');
% %     hold on
% %     p3 = plot([percentileA3 percentileA3],y3,'LineWidth',2.5);
% %     LegPerc3 = strcat('95th percentile=', num2str(percentileA3));
% %     legend(p3, LegPerc3)

    percentileVAF = prctile(VAF_Distribution-7,95)
        
    figure()
    hVAF = histogram(VAF_Distribution-7, nbins);
    %axis([0 100 0 200])
    yVAF=get(gca,'ylim');
    hold on
    pVAF = plot([percentileVAF percentileVAF],yVAF,'LineWidth',2.5);
    LegPercVAF = strcat('95th percentile=', num2str(floor(10*percentileVAF)/10));
    legend(pVAF, LegPercVAF)
    xlabel('DOF VAF Values (%)')
    ylabel('Frequency')
    
    redist = mean(DOF_VAF_Distribution, 2)+abs(min(mean(DOF_VAF_Distribution, 2)));
    redist = redist/max(redist);
    percentileDOFVAF = prctile(redist,95)
    figure()
    hDOFVAF = histogram(redist, nbins);
    %axis([0 100 0 250])
    yDOFVAF=get(gca,'ylim');
    hold on
    pDOFVAF = plot([percentileDOFVAF percentileDOFVAF],yDOFVAF,'LineWidth',2.5);
    LegPercDOFVAF = strcat('95th percentile=', num2str(floor(10*percentileDOFVAF)/10));
    legend(pDOFVAF, LegPercDOFVAF)
    xlabel('DOF VAF Values (%)')
    ylabel('Frequency')
    
    %% D) the 95th percentile is the similarity limit
    
    %look at the saved figures for the numbers.
    
    
end