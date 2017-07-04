%Operation_normalize([1,2,5,6,9,10,11,12,14,15,17,18,19,20,21])
function Operation_normalize(SubjectIDs)
    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end

        DOF = 8;
        %load the data 
        load(strcat('Processed_Subj_', SubjID, '_Left.mat'))
        load(strcat('Processed_Subj_', SubjID, '_Right.mat'))
        ProcessedRightSide(:,2:9) = ProcessedRightSide(:,2:9) + 4.5*rand(size(ProcessedRightSide(:,2:9)));
        ProcessedLeftSide(:,2:9) = ProcessedLeftSide(:,2:9) + 4.5*rand(size(ProcessedLeftSide(:,2:9)));

        for i = 2:DOF+1 %1st column is time
            ProcessedRightSide(:,i) = 100*ProcessedRightSide(:,i) / max(ProcessedRightSide(:,i)); 
            ProcessedLeftSide(:,i) = 100*ProcessedLeftSide(:,i) / max(ProcessedLeftSide(:,i)); 
        end
        
        save(strcat('Processed_Subj_', SubjID, '_Right.mat'), 'ProcessedRightSide')
        save(strcat('Processed_Subj_', SubjID, '_Left.mat'), 'ProcessedLeftSide')
        
%         save(strcat('Processed_Subj_02_Right.mat'), 'ProcessedRightSide')
%         save(strcat('Processed_Subj_02_Left.mat'), 'ProcessedLeftSide')
        
    end