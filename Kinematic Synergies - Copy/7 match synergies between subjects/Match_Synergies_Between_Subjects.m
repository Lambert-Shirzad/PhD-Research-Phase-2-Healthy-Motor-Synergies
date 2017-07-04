% Phase 1-B: Analysis step 7

% In this step, I will match the synergies identified in step 5 between
% subjects using the similarity value found in step 2. 

% 20160710 Written by Navid Lambert-Shirzad

function Match_Synergies_Between_Subjects
    
    DOF = 10;
    ndim_Global = 3;
    similarity_range_DP = 0.91;
    similarity_range_VAF = 89.7;
    similarity_range_DOF_VAF = 85.4;
    IDs = [1,3,4,6,7,8,9,10,12,13,14,15]; %  [1:18]; %
    Hand_Dominance = [1,1,0,1,1,1,1,0,1,1,1,1];%[1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1];
    RangeofMotion = [180, 180, 180, 220, 180, 110, 130, 180, 35, 140]; 
    
    CurrentDirectory = cd;
    CurrentDirectoryUp = strrep(CurrentDirectory,'7 match synergies between subjects',''); 
    OrigSynergiesFolder = strcat(CurrentDirectoryUp, '\5 extract ndim_global vectors per subject\');
    
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find the average synergy vectors for each subjects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for SubjCount=1:size(IDs,2) 
        if IDs(SubjCount) < 10
            SubjID = strcat('0', num2str(IDs(SubjCount)));
        else
            SubjID = num2str(IDs(SubjCount));
        end
        load(strcat(OrigSynergiesFolder,'Y', SubjID, '_Synergies_ndimGlobal' )); 
        %match vectors and calculate the average
        DP_All = zeros(ndim_Global, ndim_Global); %DP for dot product
        for i = 1:ndim_Global 
            for j = 1:ndim_Global
                DP_All(i,j) = Synergy_Dom(i,:)*Synergy_NonDom(j,:)'; %Synergies comes from the loaded data
            end
        end   
        NumMatched = 0;
        while NumMatched < ndim_Global
            for RowNum = 1:ndim_Global
                %find the max DP in the row
                [MaxDPRow, IndRow] = max(DP_All(RowNum,:));
                %check to see if the max in row is also max in column
                %if so it is a match
                [~, IndCol] = max(DP_All(:,IndRow));
                if MaxDPRow ~= 0 %this row is not previously matched
                    if IndCol == RowNum
                        %this is a match
                        NumMatched = NumMatched+1;
                        tempAvg = (Synergy_Dom(RowNum,:)+Synergy_NonDom(IndRow,:))/2;
                        %normalize it
                        Subj_Avg_Synergy(SubjCount,NumMatched,1:DOF) = tempAvg/norm(tempAvg);
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
    end %calculating the average synergy vector of a subject
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %is Healthy synergy template dependent on the choice of first subject to
    %whom the other subjects will be matched and averaged
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Healthy_Avg_Synergy = Subj_Avg_Synergy;
    for SubjCount1=1:size(IDs,2)
        subjOne(:,:) = Subj_Avg_Synergy(SubjCount1,1:ndim_Global,1:DOF);  %match with the original, what if we match with the average?          
        temp_Avg(:,:) = Healthy_Avg_Synergy(SubjCount1,1:ndim_Global,1:DOF);
        temp_DP = zeros(ndim_Global,1);
        for SubjCount2=1:size(IDs,2)
            subjTwo(:,:) = Subj_Avg_Synergy(SubjCount2,1:ndim_Global,1:DOF);

            DP_All = zeros(ndim_Global, ndim_Global); %DP for dot product
            for i = 1:ndim_Global 
                for j = 1:ndim_Global
                    DP_All(i,j) = subjOne(i,:)*subjTwo(j,:)'; 
                end
            end
            %matching on the rows
            NumMatched = 0;
            NumSimilar = 0;
            while NumMatched < ndim_Global
                for RowNum = 1:ndim_Global
                    %find the max DP in the row
                    [MaxDPRow, IndRow] = max(DP_All(RowNum,:));
                    %check to see if the max in row is also max in column
                    %if so it is a match
                    [~, IndCol] = max(DP_All(:,IndRow));
                    if MaxDPRow ~= 0 %this row is not previously matched
                        if IndCol == RowNum
                            %this is a match
                            if MaxDPRow >= similarity_range_DP
                                NumSimilar = NumSimilar + 1;
                            end
                            NumMatched = NumMatched+1;                       
                            temp_Avg(RowNum,:) = temp_Avg(RowNum,:)+subjTwo(IndRow,:);
                            temp_DP(RowNum) = temp_DP(RowNum)+MaxDPRow;
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
            
        end
        for i=1:ndim_Global
            Healthy_Avg_Synergy(SubjCount1,i,1:DOF) = temp_Avg(i,:)/norm(temp_Avg(i,:));
        end
        Avg_DP_Subj(:,SubjCount1) = temp_DP/size(IDs,2);
        %plot templates coming from each of the subjects as the first
        %subject to start the matching procedure
%         figure()
%         for i=1:ndim_Global
%             subplot(ndim_Global, 1, i)
%             bar(temp_Avg(i,1:DOF)/norm(temp_Avg(i,1:DOF)))
%             axis([0 DOF+1 0 1])
%         end

    end
    %Avg_DP_Subj
    % if you calculate the DP of each of these prelim templates you will
    % see ALL of the DPs are above the similarity range. So, the healthy
    % template does not depend on the choice of subject to start the
    % averaging process. Great!
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1. build the healthy population synergy template 
    %2. re-arrange the subj average synergies to match the order of the vectors
    %in the healthy synergy template (makes comparison a bit more clear as
    %they will not be shuffled anymore)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Healthy_Synergy_Template(1:ndim_Global,1:DOF) = Healthy_Avg_Synergy(1,1:ndim_Global,1:DOF);
    for SubjCount2= 1:size(IDs,2)   
        subjTwo(:,:) = Subj_Avg_Synergy(SubjCount2,1:ndim_Global,1:DOF); 
        DP_All = zeros(ndim_Global, ndim_Global); %DP for dot product
        for i = 1:ndim_Global 
            for j = 1:ndim_Global
                DP_All(i,j) = Healthy_Synergy_Template(i,:)*subjTwo(j,:)'; 
            end
        end      

        %matching on the rows
        NumMatched = 0;
        NumSimilar = 0;
        while NumMatched < ndim_Global
            for RowNum = 1:ndim_Global
                %find the max DP in the row
                [MaxDPRow, IndRow] = max(DP_All(RowNum,:));
                %check to see if the max in row is also max in column
                %if so it is a match
                [~, IndCol] = max(DP_All(:,IndRow));
                if MaxDPRow ~= 0 %this row is not previously matched
                    if IndCol == RowNum
                        %this is a match
                        if MaxDPRow >= similarity_range_DP
                            NumSimilar = NumSimilar + 1;
                        end
                        Healthy_Synergy_Template(RowNum,:) = Healthy_Synergy_Template(RowNum,:)+subjTwo(IndRow,:);
                        Subj_Avg_Synergy(SubjCount2,RowNum,1:DOF) = subjTwo(IndRow,:);
                        DP_All(RowNum,:) = zeros(size(DP_All(RowNum,:)));
                        DP_All(:,IndRow) = zeros(size(DP_All(:,IndRow)));
                        %RowNum
                        %IndRow
                        NumMatched = NumMatched+1;                       
                    end
                    if NumMatched == ndim_Global
                        break
                    end
                end

            end    
        end
        %BetweenSubjSimilarity(SubjCount2) = NumSimilar; %they are all similar
%         figure()
        for i=1:ndim_Global
            temp(1,:) = Subj_Avg_Synergy(SubjCount2,i,1:DOF);
            Subj_Avg_Synergy(SubjCount2,i,1:DOF) = temp(1,:)/norm(temp(1,:));
%             subplot(ndim_Global, 1, i)
%             bar(temp(1,:)/norm(temp(1,:)))
%             axis([0 DOF+1 0 1])
        end
    end    
    for i=1:ndim_Global
        Healthy_Synergy_Template(i,1:DOF) = Healthy_Synergy_Template(i,:)/norm(Healthy_Synergy_Template(i,:));
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compare subjects on DP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for SubjCount1=1:size(IDs,2)
        subjOne(:,:) = Healthy_Synergy_Template;%Subj_Avg_Synergy(SubjCount1,1:ndim_Global,1:DOF);  %match with the original, what if we match with the average?          
        temp_DP = zeros(ndim_Global,size(IDs,2));
        for SubjCount2=1:size(IDs,2)
            subjTwo(:,:) = Subj_Avg_Synergy(SubjCount2,1:ndim_Global,1:DOF);
            DP_All = zeros(ndim_Global, ndim_Global); %DP for dot product
            for i = 1:ndim_Global 
                for j = 1:ndim_Global
                    DP_All(i,j) = subjOne(i,:)*subjTwo(j,:)'; 
                end
            end
            %matching on the rows
            NumMatched = 0;
            NumSimilar = 0;
            while NumMatched < ndim_Global
                for RowNum = 1:ndim_Global
                    %find the max DP in the row
                    [MaxDPRow, IndRow] = max(DP_All(RowNum,:));
                    %check to see if the max in row is also max in column
                    %if so it is a match
                    [~, IndCol] = max(DP_All(:,IndRow));
                    if MaxDPRow ~= 0 %this row is not previously matched
                        if IndCol == RowNum
                            %this is a match
                            if MaxDPRow >= similarity_range_DP
                                NumSimilar = NumSimilar + 1;
                            end
                            NumMatched = NumMatched+1;                       
                            temp_DP(RowNum, SubjCount2) = MaxDPRow;%temp_DP(RowNum)+MaxDPRow;
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
            
        end
        Avg_DP_Subj(:,SubjCount1) = median(temp_DP,2);%temp_DP/size(IDs,2);
        
        %plot each subject's average synergy and median of DP with respect 
        %to the matched values
%         if SubjCount1 == 1
%             figure()
%         end
%         for i=1:ndim_Global
%             subplot(size(IDs,2),ndim_Global, ndim_Global*(SubjCount1-1)+i)
%             bar(subjOne(i,:))
%             axis([0 DOF+1 0 1])
%             title(strcat('DP = ', num2str(Avg_DP_Subj(i,SubjCount1))))
%         end

    end
    
    Avg_DP_Subj
    figure()
    for dimCount=1:ndim_Global
        for SubjCount=1:size(IDs,2)
            subplot(ndim_Global, size(IDs,2), (dimCount-1)*size(IDs,2)+SubjCount)
            Syn_vector(1,:)=Subj_Avg_Synergy(SubjCount,dimCount,1:DOF);
            bar(Syn_vector, 'k')
            axis([0 DOF+1 0 1])
            text(0, 1.05, strcat('DP= ', num2str(floor(Avg_DP_Subj(dimCount,SubjCount)*100)/100)))
            if dimCount == ndim_Global
                if IDs(SubjCount) < 10
                    SubjID = strcat('0', num2str(IDs(SubjCount)));
                else
                    SubjID = num2str(IDs(SubjCount));
                end
                xlabel(strcat('H#', SubjID))
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compare subjects on VAF and DOF_VAF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DataFolder2 = strcat(CurrentDirectoryUp, '\1 identifying synergies\');    
    for SubjCount = 1:size(IDs,2)
        %load motion data of SubjCount
        if IDs(SubjCount) < 10
            SubjID = strcat('0', num2str(IDs(SubjCount)));
        else
            SubjID = num2str(IDs(SubjCount));
        end
        load(strcat(DataFolder2,'FullSet_', 'Y', SubjID, '.mat'))
        FullSet(:,2) = FullSet(:,2) + sign(mean(FullSet(:,2)))*mean(FullSet(:,2));
        FullSet(:,3) = FullSet(:,3) + sign(mean(FullSet(:,3)))*mean(FullSet(:,3));
        FullSet(:,4) = FullSet(:,4) + sign(max(FullSet(:,4)))*max(FullSet(:,4));
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12) ]; %RightFullSet = [FullSet(:,2:6) FullSet(:,8:10)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)];%LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:18)]; %t=0s is not included in FullSet (first row is t=1/30s)

        %make sure data is non-negative
        ProcessedRightSide = repmat([90 90 90 90 0 90 0 90 10 70],size(RightFullSet,1),1) + RightFullSet; %[90 90 90 90 0 90 0 90 10 70] from OpenSim model, abs(lower bound) of each DOF
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

        subjOne(:,:) = Subj_Avg_Synergy(SubjCount,1:ndim_Global,1:DOF);
        
        %reconstruct and calculate VAF & DOF_VAF
        if Hand_Dominance(SubjCount) == 1
            %Self
            Recons_Dom_Self = (ProcessedRightSide / subjOne) * subjOne;
            Recons_NonDom_Self = (ProcessedLeftSide / subjOne) * subjOne;
            VAF_Dom_Recons_Self(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_Dom_Self).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            VAF_NonDom_Recons_Self(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_NonDom_Self).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            DOF_VAF_Dom_Recons_Self(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_Dom_Self).^2,1) ./ sum((ProcessedRightSide).^2,1));
            DOF_VAF_NonDom_Recons_Self(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_NonDom_Self).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            Avg_DOF_VAF_Dom_Recons_Self(SubjCount, 1) = mean(DOF_VAF_Dom_Recons_Self(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons_Self(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons_Self(SubjCount, :));
            %Healthy template
            Recons_Dom_Template = (ProcessedRightSide / Healthy_Synergy_Template) * Healthy_Synergy_Template;
            Recons_NonDom_Template = (ProcessedLeftSide / Healthy_Synergy_Template) * Healthy_Synergy_Template;
            VAF_Dom_Recons_Template(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_Dom_Template).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            VAF_NonDom_Recons_Template(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_NonDom_Template).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            DOF_VAF_Dom_Recons_Template(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_Dom_Template).^2,1) ./ sum((ProcessedRightSide).^2,1));
            DOF_VAF_NonDom_Recons_Template(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_NonDom_Template).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            Avg_DOF_VAF_Dom_Recons_Template(SubjCount, 1) = mean(DOF_VAF_Dom_Recons_Template(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons_Template(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons_Template(SubjCount, :));
            %by others
            for others = 1:size(IDs,2)
                subjOther(:,:) = Subj_Avg_Synergy(others,1:ndim_Global,1:DOF);
                Recons_Dom_Other = (ProcessedRightSide / subjOther) * subjOther;
                Recons_NonDom_Other = (ProcessedLeftSide / subjOther) * subjOther;
                VAF_Dom_Recons_Other_temp(others, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_Dom_Other).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
                VAF_NonDom_Recons_Other_temp(others, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_NonDom_Other).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
                DOF_VAF_Dom_Recons_Other_temp(others, :) = 100*(1 - sum((ProcessedRightSide - Recons_Dom_Other).^2,1) ./ sum((ProcessedRightSide).^2,1));
                DOF_VAF_NonDom_Recons_Other_temp(others, :) = 100*(1 - sum((ProcessedLeftSide - Recons_NonDom_Other).^2,1) ./ sum((ProcessedLeftSide).^2,1));
                Avg_DOF_VAF_Dom_Recons_Other_temp(others, 1) = mean(DOF_VAF_Dom_Recons_Other_temp(others, :));
                Avg_DOF_VAF_NonDom_Recons_Other_temp(others, 1) = mean(DOF_VAF_NonDom_Recons_Other_temp(others, :));
            end
            VAF_Dom_Recons_Other(SubjCount, 1) = mean(VAF_Dom_Recons_Other_temp);
            VAF_NonDom_Recons_Other(SubjCount, 1) = mean(VAF_NonDom_Recons_Other_temp);
            DOF_VAF_Dom_Recons_Other(SubjCount, :) = mean(DOF_VAF_Dom_Recons_Other_temp);
            DOF_VAF_NonDom_Recons_Other(SubjCount, :) = mean(DOF_VAF_NonDom_Recons_Other_temp);
            Avg_DOF_VAF_Dom_Recons_Other(SubjCount, 1) = mean(Avg_DOF_VAF_Dom_Recons_Other_temp);
            Avg_DOF_VAF_NonDom_Recons_Other(SubjCount, 1) = mean(Avg_DOF_VAF_NonDom_Recons_Other_temp);
        else
            SubjCount
            %Self
            Recons_Dom_Self = (ProcessedLeftSide / subjOne) * subjOne;
            Recons_NonDom_Self = (ProcessedRightSide / subjOne) * subjOne;
            VAF_Dom_Recons_Self(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_Dom_Self).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            VAF_NonDom_Recons_Self(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_NonDom_Self).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            DOF_VAF_Dom_Recons_Self(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_Dom_Self).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            DOF_VAF_NonDom_Recons_Self(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_NonDom_Self).^2,1) ./ sum((ProcessedRightSide).^2,1));
            Avg_DOF_VAF_Dom_Recons_Self(SubjCount, 1) = mean(DOF_VAF_Dom_Recons_Self(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons_Self(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons_Self(SubjCount, :));
            %Healthy template
            Recons_Dom_Template = (ProcessedLeftSide / Healthy_Synergy_Template) * Healthy_Synergy_Template;
            Recons_NonDom_Template = (ProcessedRightSide / Healthy_Synergy_Template) * Healthy_Synergy_Template;
            VAF_Dom_Recons_Template(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_Dom_Template).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            VAF_NonDom_Recons_Template(SubjCount, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_NonDom_Template).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            DOF_VAF_Dom_Recons_Template(SubjCount, :) = 100*(1 - sum((ProcessedLeftSide - Recons_Dom_Template).^2,1) ./ sum((ProcessedLeftSide).^2,1));
            DOF_VAF_NonDom_Recons_Template(SubjCount, :) = 100*(1 - sum((ProcessedRightSide - Recons_NonDom_Template).^2,1) ./ sum((ProcessedRightSide).^2,1));
            Avg_DOF_VAF_Dom_Recons_Template(SubjCount, 1) = mean(DOF_VAF_Dom_Recons_Template(SubjCount, :));
            Avg_DOF_VAF_NonDom_Recons_Template(SubjCount, 1) = mean(DOF_VAF_NonDom_Recons_Template(SubjCount, :));
            %by others
            for others = 1:size(IDs,2)
                subjOther(:,:) = Subj_Avg_Synergy(others,1:ndim_Global,1:DOF);
                Recons_Dom_Other = (ProcessedLeftSide / subjOther) * subjOther;
                Recons_NonDom_Other = (ProcessedRightSide / subjOther) * subjOther;
                VAF_Dom_Recons_Other_temp(others, 1) = 100*(1 - (sum(sum((ProcessedLeftSide - Recons_Dom_Other).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
                VAF_NonDom_Recons_Other_temp(others, 1) = 100*(1 - (sum(sum((ProcessedRightSide - Recons_NonDom_Other).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
                DOF_VAF_Dom_Recons_Other_temp(others, :) = 100*(1 - sum((ProcessedLeftSide - Recons_Dom_Other).^2,1) ./ sum((ProcessedLeftSide).^2,1));
                DOF_VAF_NonDom_Recons_Other_temp(others, :) = 100*(1 - sum((ProcessedRightSide - Recons_NonDom_Other).^2,1) ./ sum((ProcessedRightSide).^2,1));
                Avg_DOF_VAF_Dom_Recons_Other_temp(others, 1) = mean(DOF_VAF_Dom_Recons_Other_temp(others, :));
                Avg_DOF_VAF_NonDom_Recons_Other_temp(others, 1) = mean(DOF_VAF_NonDom_Recons_Other_temp(others, :));
            end
            VAF_Dom_Recons_Other(SubjCount, 1) = mean(VAF_Dom_Recons_Other_temp);
            VAF_NonDom_Recons_Other(SubjCount, 1) = mean(VAF_NonDom_Recons_Other_temp);
            DOF_VAF_Dom_Recons_Other(SubjCount, :) = mean(DOF_VAF_Dom_Recons_Other_temp);
            DOF_VAF_NonDom_Recons_Other(SubjCount, :) = mean(DOF_VAF_NonDom_Recons_Other_temp);
            Avg_DOF_VAF_Dom_Recons_Other(SubjCount, 1) = mean(Avg_DOF_VAF_Dom_Recons_Other_temp);
            Avg_DOF_VAF_NonDom_Recons_Other(SubjCount, 1) = mean(Avg_DOF_VAF_NonDom_Recons_Other_temp);
        end
    end
    
    figure()
    
    subplot(3,2,1)
    boxplot([VAF_Dom_Recons_Self, VAF_Dom_Recons_Other, VAF_Dom_Recons_Template],'Colors','k')
    ylabel('VAF (%)', 'FontSize',11)
    hold on
    p1 = plot([0.5 3.5],[similarity_range_VAF similarity_range_VAF],'k-');
    Leg1 = 'VAF Similarity Limit = 89.7%';
    legend(p1, Leg1)
    axis([0.5 3.5 88 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Self' 'Others' 'Template'} ...
        ,'fontsize',11)
    title('Dominant Limb')

    subplot(3,2,2)
    boxplot([VAF_NonDom_Recons_Self, VAF_NonDom_Recons_Other, VAF_NonDom_Recons_Template],'Colors','k')
    hold on
    p2 = plot([0.5 3.5],[similarity_range_VAF similarity_range_VAF],'k-');
    Leg2 = 'VAF Similarity Limit = 89.7%';
    legend(p2, Leg2)
    axis([0.5 3.5 88 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Self' 'Others' 'Template'} ...
        ,'fontsize',11)
    title('Non-Dominant Limb')

    subplot(3,2,3)
    boxplot([Avg_DOF_VAF_Dom_Recons_Self, Avg_DOF_VAF_Dom_Recons_Other, Avg_DOF_VAF_Dom_Recons_Template],'Colors','k')
    ylabel('Average DOF VAF (%)', 'FontSize',11)
    hold on
    p3 = plot([0.5 3.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
    Leg3 = 'Average DOF VAF Similarity Limit = 85.4%';
    legend(p3, Leg3)
    axis([0.5 3.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Self' 'Others' 'Template'} ...
        ,'fontsize',11)

    subplot(3,2,4)
    boxplot([Avg_DOF_VAF_NonDom_Recons_Self, Avg_DOF_VAF_NonDom_Recons_Other, Avg_DOF_VAF_NonDom_Recons_Template],'Colors','k')
    hold on
    p4 = plot([0.5 3.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
    Leg4 = 'Average DOF VAF Similarity Limit = 85.4%';
    legend(p4, Leg4)
    axis([0.5 3.5 84 100])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Self' 'Others' 'Template'} ...
        ,'fontsize',11)

    subplot(3,2,5)
    bar([mean(DOF_VAF_Dom_Recons_Self); mean(DOF_VAF_Dom_Recons_Other); mean(DOF_VAF_Dom_Recons_Template)]', 'grouped')
    axis([0 DOF+1 0 149])
    ylabel('DOF VAF (%) by Muscles', 'FontSize',11)
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',10)
    legend('Self', 'Others', 'Template', 'fontsize',11)
    colormap(gray)
%     hold on
%     p5 = plot([0 9.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
%     Leg5 = 'DOF VAF Similarity Limit = 80.5%';
%     legend(p5, Leg5)

    subplot(3,2,6)
    bar([mean(DOF_VAF_NonDom_Recons_Self); mean(DOF_VAF_NonDom_Recons_Other); mean(DOF_VAF_NonDom_Recons_Template)]', 'grouped')
    axis([0 DOF+1 0 149])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',10)
    legend('Self', 'Others', 'Template', 'fontsize',11)
%     hold on
%     p6 = plot([0 9.5],[similarity_range_DOF_VAF similarity_range_DOF_VAF],'k-');
%     Leg6 = 'DOF VAF Similarity Limit = 80.5%';
%     legend(p6, Leg6)
    
    figure()
    boxplot(Avg_DP_Subj','Colors','k') %sort(Avg_DP_Subj, 2, 'descend')
    ylabel('Dot Product of Matched Synergy Vectors', 'FontSize',11)
    hold on
    p = plot([0.5 ndim_Global+0.5],[similarity_range_DP similarity_range_DP],'k-');
    Leg = 'Dot Product Similarity Limit = 0.91';
    legend(p, Leg)
    axis([0.5 ndim_Global+0.5 0.89 1])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'1st pair' '2nd pair' '3rd pair' } ,'fontsize',11)

    
    DataFolder2 = strcat(CurrentDirectoryUp, '\1 identifying synergies\');    
    for SubjCount = 1:size(IDs,2)
        %load motion data of SubjCount
        if IDs(SubjCount) < 10
            SubjID = strcat('0', num2str(IDs(SubjCount)));
        else
            SubjID = num2str(IDs(SubjCount));
        end
        load(strcat(DataFolder2,'FullSet_', 'Y', SubjID, '.mat'))
        FullSet(:,2) = FullSet(:,2) + sign(mean(FullSet(:,2)))*mean(FullSet(:,2));
        FullSet(:,3) = FullSet(:,3) + sign(mean(FullSet(:,3)))*mean(FullSet(:,3));
        FullSet(:,4) = FullSet(:,4) + sign(max(FullSet(:,4)))*max(FullSet(:,4));
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12) ]; %RightFullSet = [FullSet(:,2:6) FullSet(:,8:10)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)];%LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:18)]; %t=0s is not included in FullSet (first row is t=1/30s)

        %make sure data is non-negative
        ProcessedRightSide = repmat([90 90 90 90 0 90 0 90 10 70],size(RightFullSet,1),1) + RightFullSet; %[90 90 90 90 0 90 0 90 10 70] from OpenSim model, abs(lower bound) of each DOF
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
            for i = 1:ndim_Global
                recon = (ProcessedRightSide / Healthy_Synergy_Template(i,:)) * Healthy_Synergy_Template(i,:);
                VAF_Vector(i,SubjCount)=100*(1 - (sum(sum((ProcessedRightSide - recon).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));                
            end
            recon = (ProcessedRightSide / Healthy_Synergy_Template(:,:)) * Healthy_Synergy_Template(:,:);
            VAF_total=100*(1 - (sum(sum((ProcessedRightSide - recon).^2,2),1)) / (sum(sum((ProcessedRightSide).^2,2),1)));
            total = sum(VAF_Vector(:,SubjCount));
            for i = 1:ndim_Global
                VAF_Vector(i,SubjCount)=VAF_total*VAF_Vector(i,SubjCount)/total;                
            end
        else
            for i = 1:ndim_Global
                recon = (ProcessedLeftSide / Healthy_Synergy_Template(i,:)) * Healthy_Synergy_Template(i,:);
                VAF_Vector(i,SubjCount)=100*(1 - (sum(sum((ProcessedLeftSide - recon).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));                
            end
            recon = (ProcessedLeftSide / Healthy_Synergy_Template(:,:)) * Healthy_Synergy_Template(:,:);
            VAF_total=100*(1 - (sum(sum((ProcessedLeftSide - recon).^2,2),1)) / (sum(sum((ProcessedLeftSide).^2,2),1)));
            total = sum(VAF_Vector(:,SubjCount));
            for i = 1:ndim_Global
                VAF_Vector(i,SubjCount)=VAF_total*VAF_Vector(i,SubjCount)/total;                
            end
        end            
    end
    
    figure()
    for i=1:ndim_Global
        Healthy_Synergy_Template(i,1:DOF) = Healthy_Synergy_Template(i,:)/norm(Healthy_Synergy_Template(i,:));
        subplot(ndim_Global, 2, 2*i-1)
        
        bar(Healthy_Synergy_Template(i,1:DOF))
        ylabel(strcat('Synergy Vector #', num2str(i)))
        if i==1
            title('Healthy Synergy Template')
        end
        axis([0 DOF+1 0 1])
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit'  ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',10)
        colormap(gray)
        
        subplot(ndim_Global, 2, 2*i)
        boxplot(VAF_Vector(i,:),'Colors','k')
        ylabel(strcat('VAF of Synergy #', num2str(i)))
        if i==1
            title('VAF Contribution of Each Synergy Vector')
        end
        axis([0.85 1.15 28 38])
    end
    
    %%%%%%%%%%%%%
    %experimental stuff
    %%%%%%%%%%%%%
    
%     figure()
%     for i=1:ndim_Global
%        subplot(ndim_Global, 4, (i-1)*4+1 )
%        bar(Healthy_Synergy_Template(i,1:DOF))
%        axis([0 DOF+1 0 1])
%     end
%     for SubjCount=16:18 
%         if IDs(SubjCount) < 10
%             SubjID = strcat('0', num2str(IDs(SubjCount)));
%         else
%             SubjID = num2str(IDs(SubjCount));
%         end
%         load(strcat(OrigSynergiesFolder,'Y', SubjID, '_Synergies_ndimGlobal' )); 
%         %match vectors and calculate the average
%         DP_All_NonImp = zeros(ndim_Global, ndim_Global); %DP for dot product
%         for i = 1:ndim_Global 
%             for j = 1:ndim_Global
%                 DP_All_NonImp(i,j) = Healthy_Synergy_Template(i,:)*Synergy_Dom(j,:)'; %Synergies comes from the loaded data
%             end
%         end   
%         NumMatched = 0;
%         while NumMatched < ndim_Global
%             for RowNum = 1:ndim_Global
%                 %find the max DP in the row
%                 [MaxDPRow, IndRow] = max(DP_All_NonImp(RowNum,:));
%                 %check to see if the max in row is also max in column
%                 %if so it is a match
%                 [~, IndCol] = max(DP_All_NonImp(:,IndRow));
%                 if MaxDPRow ~= 0 %this row is not previously matched
%                     if IndCol == RowNum
%                         %this is a match
%                         NumMatched = NumMatched+1;
%                         subplot(ndim_Global, 4, (NumMatched-1)*4+(SubjCount-14) )
%                         bar(Synergy_Dom(IndRow,1:DOF))
%                         title(num2str(MaxDPRow))
%                         axis([0 DOF+1 0 1])
%                         %normalize it
%                         
%                         DP_All_NonImp(RowNum,:) = zeros(size(DP_All_NonImp(RowNum,:)));
%                         DP_All_NonImp(:,IndRow) = zeros(size(DP_All_NonImp(:,IndRow)));
%                         %RowNum
%                         %IndRow                                               
%                     end
%                     if NumMatched == ndim_Global
%                         break
%                     end
%                 end
%                 
%             end    
%         end    
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    for i = 1:ndim_Global 
        for j = 1:DOF
            normalize_Healthy_Synergy_Template(i,j) = Healthy_Synergy_Template(i,j)*RangeofMotion(j)/220; %
        end
    end
    
    figure()
    for i=1:ndim_Global
        normalize_Healthy_Synergy_Template(i,1:DOF) = normalize_Healthy_Synergy_Template(i,:)/norm(normalize_Healthy_Synergy_Template(i,:));
        subplot(ndim_Global, 2, 2*i-1)
        bar(normalize_Healthy_Synergy_Template(i,1:DOF))
        ylabel(strcat('Synergy Vector #', num2str(i)))
        if i==1
            title('Healthy Synergy Template')
        end
        axis([0 DOF+1 0 1])
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', {'TrRol'  'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',10)
        colormap(gray)
        
        subplot(ndim_Global, 2, 2*i)
        boxplot(VAF_Vector(i,:),'Colors','k')
        ylabel(strcat('VAF of Synergy #', num2str(i)))
        if i==1
            title('VAF Contribution of Each Synergy Vector')
        end
        axis([0.85 1.15 17 29])
    end
    
    for i = 1:ndim_Global 
        constant = 0.5/Healthy_Synergy_Template(i,1);
        for j = 1:DOF
            poses(i,j) = constant*Healthy_Synergy_Template(i,j)*RangeofMotion(j); %
        end
        poses(i,:) = poses(i,:) - [90 90 90 90 0 90 0 90 10 70];
        
        AllGood = poses(i,:)- [90 90 90 130 180 20 130 90 25 70];
        if max(AllGood) > 0
            [~,colnum] = max(AllGood);
            Healthy_Synergy_Template(i,colnum)
            constant = 1/Healthy_Synergy_Template(i,colnum)-0.1
            constant*Healthy_Synergy_Template(i,:)
            for j = 1:DOF
                poses(i,j) = constant*Healthy_Synergy_Template(i,j)*RangeofMotion(j); %
            end
            poses(i,:) = poses(i,:) - [90 90 90 90 0 90 0 90 10 70];
        end
        
    end
    poses
    Healthy_Synergy_Template
end %function
    
    
    
  