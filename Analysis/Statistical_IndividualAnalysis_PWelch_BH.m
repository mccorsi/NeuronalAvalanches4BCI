function Statistical_IndividualAnalysis_PWelch_BH(subject_ID,Data_MEG_MI_DK_indiv, Data_MEG_Baseline_DK_indiv, db_path, numperm, nb_ROIs_DK, p_thresh, freqband)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform individual-level statistical analysis - PWelch
    % Authors: MCC
    % Date: 25/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; %clc;

Node_MI=Data_MEG_MI_DK_indiv;
Node_Base=Data_MEG_Baseline_DK_indiv; % nb_ROIs x nb_trials

nb_trials_V2=min(size(Node_MI, 2), size(Node_Base,2));

output=[];
mask=logical(triu(ones(nb_ROIs_DK),1));
for kk_trials=1:nb_trials_V2
    %nodes
    temp_node_200=Node_Base(:,kk_trials);
    temp_node_100=Node_MI(:,kk_trials);
    output_node(:,kk_trials,1)=temp_node_100;
    output_node(:,kk_trials,2)=temp_node_200;
end

% % nodes
output_node_MI=sum(output_node(:,:,1),2)./sum(output_node(:,:,1)~=0,2);
output_node_MI(isnan(output_node_MI))=0;
output_node_Base=sum(output_node(:,:,2),2)./sum(output_node(:,:,2)~=0,2);
output_node_Base(isnan(output_node_Base))=0;   

obs_node_diff= output_node_MI-output_node_Base;

% function to randomize
group_generator = @(x)x(randperm(numel(x)));

% nodes
rep_obs_node_diff=[];
rep_obs_node_diff=repmat(obs_node_diff,1,numperm);

output_node_cat=cat(2,output_node(:,:,1),output_node(:,:,2));

for kk_nperm=1:numperm
    % nodes
    temp_node=output_node_cat(:,group_generator(1:size(output_node_cat,2)));
    MI_temp_node=sum(temp_node(:,1:nb_trials_V2),2)./sum(temp_node(:,1:nb_trials_V2)~=0,2);
    MI_temp_node(isnan(MI_temp_node))=0;
    Base_temp_node=sum(temp_node(:,nb_trials_V2+1:end),2)./sum(temp_node(:,nb_trials_V2+1:end)~=0,2);
    Base_temp_node(isnan(Base_temp_node))=0;
    perm_node_diff(:,kk_nperm)=MI_temp_node- Base_temp_node;
end

% nodes 
nb_times=zeros(size(perm_node_diff,1),numperm);
nb_times(perm_node_diff>rep_obs_node_diff)=1; 
pval_node_diff_perm_uncorrected=sum(nb_times,2)./numperm;
pval_node_diff_perm_corrected_BH=pval_adjust(pval_node_diff_perm_uncorrected, 'BH')<p_thresh ;  

clearvars nb_times sum_times
  
save(strcat(db_path,'Stat_Analysis_Indiv_PWelch_',subject_ID,'freq','_',num2str(freqband(1)), '_', num2str(freqband(2)),'_BH_p_thresh_',num2str(p_thresh),'.mat'),...
                'pval_node_diff_perm_uncorrected',...
                'pval_node_diff_perm_corrected_BH',...
                '-v7.3');
       
end