function Statistical_IndividualAnalysis_PLV_BH(PLV_MI, PLV_Baseline, subject_ID, fig_path, db_path, numperm, nb_ROIs_DK, p_thresh, freqband)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform individual-level statistical analysis - PLV
    % Authors: MCC
    % Date: 26/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; %clc;

% concat. x avere nb_ROIs x nb_ROIs x nb_trials

data_concat=cat(3,PLV_MI, PLV_Baseline);

MI=PLV_MI;
Node_MI=squeeze(sum(MI,2)); % nb_ROIs x nb_trials
Base=PLV_Baseline;
Node_Base=squeeze(sum(Base,2)); % nb_ROIs x nb_trials

nb_trials_V2=min(size(MI, 3), size(Base,3));

output=[];
mask=logical(triu(ones(nb_ROIs_DK),1));
for kk_trials=1:nb_trials_V2
    % edges
    temp_200=Base(:,:,kk_trials);
    temp_100=MI(:,:,kk_trials);
    output(:,kk_trials,1)=temp_100(mask);
    output(:,kk_trials,2)=temp_200(mask);
    %nodes
    temp_node_200=Node_Base(:,kk_trials);
    temp_node_100=Node_MI(:,kk_trials);
    output_node(:,kk_trials,1)=temp_node_100;
    output_node(:,kk_trials,2)=temp_node_200;
end

% edges
output_MI=sum(output(:,:,1),2)./sum(output(:,:,1)~=0,2);
output_MI(isnan(output_MI))=0;
output_Base=sum(output(:,:,2),2)./sum(output(:,:,2)~=0,2);
output_Base(isnan(output_Base))=0;

obs_diff= output_MI-output_Base;
obs_abs_diff=abs(obs_diff);

% % nodes
output_node_MI=sum(output_node(:,:,1),2)./sum(output_node(:,:,1)~=0,2);
output_node_MI(isnan(output_node_MI))=0;
output_node_Base=sum(output_node(:,:,2),2)./sum(output_node(:,:,2)~=0,2);
output_node_Base(isnan(output_node_Base))=0;

obs_node_diff= output_node_MI-output_node_Base;
obs_node_abs_diff=abs(obs_node_diff);

% function to randomize
group_generator = @(x)x(randperm(numel(x)));

% edges
rep_obs_abs_diff=[];
rep_obs_abs_diff=repmat(obs_abs_diff,1,numperm);

rep_obs_diff=[];
rep_obs_diff=repmat(obs_abs_diff,1,numperm);

output_cat=cat(2,output(:,:,1),output(:,:,2));%output_MI, output_Base);

% nodes
rep_obs_node_abs_diff=[];
rep_obs_node_abs_diff=repmat(obs_node_abs_diff,1,numperm);

rep_obs_node_diff=[];
rep_obs_node_diff=repmat(obs_node_abs_diff,1,numperm);

output_node_cat=cat(2,output_node(:,:,1),output_node(:,:,2));

for kk_nperm=1:numperm
    % edges
    temp=output_cat(:,group_generator(1:size(output_cat,2)));
    MI_temp=sum(temp(:,1:nb_trials_V2),2)./sum(temp(:,1:nb_trials_V2)~=0,2);
    MI_temp(isnan(MI_temp))=0;
    Base_temp=sum(temp(:,nb_trials_V2+1:end),2)./sum(temp(:,nb_trials_V2+1:end)~=0,2);
    Base_temp(isnan(Base_temp))=0;
    perm_abs_diff(:,kk_nperm)=abs(MI_temp- Base_temp);
    perm_diff(:,kk_nperm)=MI_temp- Base_temp;

    % nodes
    temp_node=output_node_cat(:,group_generator(1:size(output_node_cat,2)));
    MI_temp_node=sum(temp_node(:,1:nb_trials_V2),2)./sum(temp_node(:,1:nb_trials_V2)~=0,2);
    MI_temp_node(isnan(MI_temp_node))=0;
    Base_temp_node=sum(temp_node(:,nb_trials_V2+1:end),2)./sum(temp_node(:,nb_trials_V2+1:end)~=0,2);
    Base_temp_node(isnan(Base_temp_node))=0;
    perm_node_abs_diff(:,kk_nperm)=abs(MI_temp_node- Base_temp_node);
    perm_node_diff(:,kk_nperm)=MI_temp_node- Base_temp_node;
end

% edges
nb_times=zeros(size(perm_abs_diff,1),numperm);
nb_times(perm_abs_diff>rep_obs_abs_diff)=1;
sum_times=sum(nb_times,2)./numperm;

mask=logical(triu(ones(nb_ROIs_DK),1));
temp=pval_adjust(sum_times, 'BH')<0.05;
pval_abs_diff_perm_corrected_BH=zeros(nb_ROIs_DK);
pval_abs_diff_perm_corrected_BH(mask)=temp;
pval_abs_diff_perm_corrected_BH=(pval_abs_diff_perm_corrected_BH+pval_abs_diff_perm_corrected_BH')./(eye(nb_ROIs_DK)+1);

mask=logical(triu(ones(nb_ROIs_DK),1));
temp = pval_adjust(sum_times, 'BH')<0.05;
pval_diff_perm_corrected_BH=zeros(nb_ROIs_DK);
pval_diff_perm_corrected_BH(mask)=temp;
pval_diff_perm_corrected_BH=(pval_diff_perm_corrected_BH+pval_diff_perm_corrected_BH')./(eye(nb_ROIs_DK)+1);


% Plots to have a rough idea of the results at the individual level
    % edges
figure(1);
subplot(1,2,1);
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,pval_diff_perm_corrected_BH);%DoMyMatrix_FromUpper(pval_diff_perm_corrected_BH,nb_ROIs_DK))
title('Significant MI vs Rest - diff - BH corrected')
axis square
colorbar
caxis([0 0.05])
subplot(1,2,2);
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,pval_abs_diff_perm_corrected_BH);%DoMyMatrix_FromUpper(pval_abs_diff_perm_corrected_BH,nb_ROIs_DK))
title('Significant MI vs Rest - abs(diff) - BH corrected')
axis square
colorbar
caxis([0 0.05])
filename=strcat(fig_path,'MI_vs_Rest_Edges_PLV_diff_BH_corrected_p_thresh_', num2str(p_thresh), '_freq_',num2str(freqband(1)), '_', num2str(freqband(2)),'_',subject_ID);
saveas(gcf,strcat(filename,'.pdf'));


save(strcat(db_path,'Stat_Analysis_Indiv_PLV_',subject_ID,'freq','_',num2str(freqband(1)), '_', num2str(freqband(2)),'_BH_p_thresh_',num2str(p_thresh),'.mat'),...
                'pval_abs_diff_perm_corrected_BH',... % BH
                'pval_diff_perm_corrected_BH',... % BH
                'obs_diff',...
                'obs_abs_diff',...
                'obs_node_diff',...
                'obs_node_abs_diff',...
                '-v7.3');

end
