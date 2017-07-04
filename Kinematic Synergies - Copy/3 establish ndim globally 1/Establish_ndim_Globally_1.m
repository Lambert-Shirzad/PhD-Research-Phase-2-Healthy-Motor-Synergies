% Phase 1-B: Analysis step 3

% In step 1, I factorized the muscle/kinematic data of healthy participants
% and generated a set of muscle synergy vectors and numerical values for
% goodness of fit. In step 2, the similarity lower bound for dot products
% of two vectors was established (95th percentile of dot product values =
% 0.8240). 
% In this step, I will look into how many synergy vectors (ndim_Global) is 
% enough to achieve the training criteria when averaging the performance of
% ndim_Globaldim vectors over the population. This will be used to extract 
% the same number of synergies for all subjects, making it easy to compare 
% between subjects.
% If Dom and NonDom ndim_Globals are the same, then it makes it easier to
% compare within subjects, between the two sides of the body. 

% 20160624 Written by Navid Lambert-Shirzad

function Establish_ndim_Globally_1
    
    DOF = 10;
    IDs = [1:15];
    NumParticipants = size(IDs,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the stats of subject-individual extraction of synergies and
    %plot them
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CurrentDirectory = cd;
    CurrentDirectoryUp = strrep(CurrentDirectory,'3 establish ndim globally 1',''); 
    OrigSynergiesFolder = strcat(CurrentDirectoryUp, '\1 identifying synergies\');
    NumDomR = 0; NumDomL = 0;
    for i=1:NumParticipants
        if IDs(i) < 10
            SubjID = strcat('0', num2str(IDs(i)));
        else
            SubjID = num2str(IDs(i));
        end
        load(strcat(OrigSynergiesFolder,'Y', SubjID, '_Dom_Synergies.mat' )); 
        load(strcat(OrigSynergiesFolder,'Y', SubjID, '_NonDom_Synergies.mat' )); 
        
        if SubjDomHand == 'R'
            ndim_Dom(i,1) = ndim_R;
            VAF_Dom(i,1) = VAF_R;
            DOF_VAF_Dom(i,:) = DOF_VAF_R;
            ndim_NonDom(i,1) = ndim_L;
            VAF_NonDom(i,1) = VAF_L;
            DOF_VAF_NonDom(i,:) = DOF_VAF_L;
            NumDomR = NumDomR+1;
        else
            ndim_Dom(i,1) = ndim_L;
            VAF_Dom(i,1) = VAF_L;
            DOF_VAF_Dom(i,:) = DOF_VAF_L;
            ndim_NonDom(i,1) = ndim_R;
            VAF_NonDom(i,1) = VAF_R;
            DOF_VAF_NonDom(i,:) = DOF_VAF_R;
            NumDomL = NumDomL+1;
        end
    end
    
    ndim_Dom_Categ = zeros(1,DOF);
    ndim_NonDom_Categ = zeros(1,DOF);
    for i=1:DOF
        for j=1:NumParticipants
           if ndim_Dom(j,1) == i
               ndim_Dom_Categ(1,i) = ndim_Dom_Categ(1,i)+1;
           end
           if ndim_NonDom(j,1) == i
               ndim_NonDom_Categ(1,i) = ndim_NonDom_Categ(1,i)+1;
           end
        end
    end
    
   figure()
    subplot(2,2,1)
    bar([100*ndim_Dom_Categ/NumParticipants; 100*ndim_NonDom_Categ/NumParticipants]', 1, 'grouped' );
    axis([0 DOF+1 0 105])
    ylabel('Percentage of Participants', 'FontSize',12)
    xlabel('Number of Synergies', 'FontSize',12)
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'},'fontsize',12)
    legend('Dominant Limb','Non-Dominant Limb', 'fontsize',12)
    colormap(gray)

    subplot(2,2,2)
    bStacked = bar([100*ndim_Dom_Categ/NumParticipants; 100*ndim_NonDom_Categ/NumParticipants], 0.5, 'stacked' );
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Dominant Limb' 'Non-Dominant Limb'},'fontsize',12)
    numVector=num2str([1:DOF].','ndim = %d');
	legend(numVector,'location','northeast')
    ylabel('Percentage of Participants', 'FontSize',12)
    axis([0 3 0 105])
    
    subplot(2,2,3)
    boxplot(DOF_VAF_Dom,'Colors','k')
    hold on
    p1 = plot([0 DOF+1],[65 65],'k-');
    LegDom = 'Min required DOF VAF = 65%';
    legend(p1, LegDom)
    text(4, 50,strcat('Mean Global VAF = ', num2str(ceil(10*mean(VAF_Dom))/10)))
    ylabel('DOF VAF: Dominant Limb (%)', 'FontSize',12)
    axis([0 DOF+1 25 105])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx' },'fontsize',10)


    subplot(2,2,4)
    boxplot(DOF_VAF_NonDom,'Colors','k')
    hold on
    p2 = plot([0 DOF+1],[65 65],'k-');
    LegNonDom = 'Min required DOF VAF = 65%';
    legend(p2, LegNonDom)
    text(4, 50,strcat('Mean Global VAF = ', num2str(floor(10*mean(VAF_NonDom))/10)))
    ylabel('DOF VAF: Non-Dominant Limb (%)', 'FontSize',12)
    axis([0 DOF+1 25 105])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'TrRol' 'TrYaw' 'TrPit' ...
        'ShFlEx' 'ShAbAd' 'ShRot' 'ElFlEx' 'ElPrSu' 'WrDev' 'WrFlEx'},'fontsize',10)
 
    
 
    
end %function
    
    
    
  