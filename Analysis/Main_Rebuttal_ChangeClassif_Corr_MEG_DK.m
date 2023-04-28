%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN - Analysi fatta sul dataset a disposizione - 20 soggetti - Sess. 4 - MEG - Desikan-Killiany
    % Authors: MCC
    % Date: 21/03/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% to be updated

root_path='/Users/marie-constance.corsi/Documents/GitHub/NeuronalAvalanches4BCI/';

%%
db_path=strcat(root_path,'Database/');
fig_path=strcat(root_path,'Figures/');

cd(strcat(root_path,'/Scripts'))
Start_SessionMatlab

subject_IDs={'',... % list of subjects_IDs
    };
nb_subj=size(subject_IDs,2);


load(strcat(root_path,'/Scripts/Visualization_Cortex/cortex_15002V_MNI.mat'));
idx_DK=3;
nb_ROIs_DK=size(cortex_15002V_MNI.Atlas(idx_DK).Scouts,2);
labels_DK={cortex_15002V_MNI.Atlas(idx_DK).Scouts(1:nb_ROIs_DK).Label};

modality='MEG';
atlas='DK';
fs=250;


%% Avalanches computed from python code


%% Statistical Analysis - Individual
clc
numperm=10000;
list_p_thresh=[0.05;0.025;0.05];%0.025;%0.05; % reliable analysis
p_thresh_stat=0.05; % individual stat
nb_perm_rand=10000;%200
        
for kk_subj=1:nb_subj % ogni soggetto
    filename_Avalanches_subj=strcat(db_path,'/Classification/3_Dataset-netbci-meg-sess4-DK/Avalanches_Analysis_LDA_ds_subj_',num2str(kk_subj),'.mat');
    filename_Avalanches_subj_file= strcat('Avalanches_Analysis_LDA_ds_subj_',subject_IDs{kk_subj});

    disp(strcat('%%%%%%%%%%%%%%%%%% Individual stat analysis - MEG - DK - LDA -',subject_IDs{kk_subj}, ' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
    % individual-level
    Statistical_IndividualAnalysis_LDA_BH(filename_Avalanches_subj, filename_Avalanches_subj_file, fig_path, db_path, numperm, nb_ROIs_DK, p_thresh_stat, cortex_15002V_MNI, idx_DK);
    % prepare materials to compute reliability:
    load(filename_Avalanches_subj);
    load(strcat(db_path,'Stat_Analysis_Indiv_LDA_',filename_Avalanches_subj_file,'_BH_p_thresh_',num2str(p_thresh_stat),'.mat'));
    matrix_signif_diff_LDA(:,:,kk_subj)=pval_diff_perm_corrected_BH; % analisi BH
    matrix_signif_abs_diff_LDA(:,:,kk_subj)=pval_abs_diff_perm_corrected_BH; % analisi BH
end

    group_generator = @(x)x(randperm(numel(x)));
    temp_rr=[];
    temp_rr_abs=[];
    % reliablity over the subjects - group level
    for kk1=1:nb_perm_rand
        for kk_subj=1:nb_subj
            rand_pos=group_generator(1:nb_ROIs_DK);
            newrand_abs(:,:,kk_subj)=matrix_signif_abs_diff_LDA(rand_pos,rand_pos,kk_subj);
            newrand(:,:,kk_subj)=matrix_signif_diff_LDA(rand_pos,rand_pos,kk_subj);

        end
        rr_abs(:,:,kk1)=sum(newrand_abs,3);
        rr(:,:,kk1)=sum(newrand,3);

        temp_rr(:,kk1)=sum(rr(:,:,kk1)<0.05);
        temp_rr_abs(:,kk1)=sum(rr_abs(:,:,kk1)<0.05);
    end
    rr_thresh=mean(mean(temp_rr));
    rr_abs_thresh=mean(mean(temp_rr_abs));

    for kk_roi=1:nb_ROIs_DK
        temp_node_abs=sum(sum(matrix_signif_abs_diff_LDA,3),1);
        pval_node_abs_reliab_uncorrect(kk_roi)=length(find(temp_rr_abs(kk_roi,:)>temp_node_abs(kk_roi)))/nb_perm_rand;
        temp_node=sum(sum(matrix_signif_diff_LDA,3),1);
        pval_node_reliab_uncorrect(kk_roi)=length(find(temp_rr(kk_roi,:)>temp_node(kk_roi)))/nb_perm_rand;
    end

    % correction - BH
    pval_node_abs_reliab_correct=pval_adjust(pval_node_abs_reliab_uncorrect,'BH');
    pval_node_reliab_correct=pval_adjust(pval_node_reliab_uncorrect,'BH');


% edges
test_abs=sum(matrix_signif_abs_diff_LDA,3);
test2_abs=repmat(test_abs, 1, 1, nb_perm_rand);
abs_cmp=rr_abs>test2_abs; 
signif_abs_cmp=pval_adjust(sum(abs_cmp,3)./numperm,'BH')<0.05; % BH


test=sum(matrix_signif_diff_LDA,3);
test2=repmat(test, 1, 1, nb_perm_rand);
cmp=rr>test2; 
signif_cmp=pval_adjust(sum(cmp,3)./numperm,'BH')<0.05;


% nodes - reliability based on edges results
node_cmp=sum(signif_cmp,2);
node_abs_cmp=sum(signif_abs_cmp,2);


save(strcat(db_path,'Rebuttal_Corr_LDA_Reliable_Analysis_Edge_Based_BH_p_thresh_edges_', num2str(0.05),'_p_thresh_nodes_',num2str(0.05),'_Group_''.mat'),...
        'signif_abs_cmp',... % post-correction bh
        'signif_cmp',... % post-correction bh
        'node_cmp',... % node signif
        'node_abs_cmp',... % node signif
        'rr_thresh',... % thresh for the nodes selection
        'rr_abs_thresh',...
        'pval_node_abs_reliab_uncorrect',...
        'pval_node_reliab_uncorrect',...
        'pval_node_abs_reliab_correct',...
        'pval_node_reliab_correct',...
        '-v7.3');

close all;
figure(1);
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,signif_cmp)
title('Significant diff - BH corrected')
axis square
filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_diff_LDA_Sess4_MEG_DK_p_thresh_', num2str(0.05),'.png');
print('-f1', filename, '-dpng','-r600');
figure(2)
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,signif_abs_cmp)
title('Significant abs(diff) - BH corrected')
axis square
filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_abs_diff_LDA_Sess4_MEG_DK_p_thresh_', num2str(0.05), '.png');
print('-f2', filename, '-dpng','-r600');

DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, node_abs_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_abs_diff_LDA_Sess4_MEG_DK_p_thresh_', num2str(0.05)), 'Sess4_MEG_DK')
DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, node_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_diff_LDA_Sess4_MEG_DK_p_thresh_', num2str(0.05)), 'Sess4_MEG_DK')


%% with SVM this time...
for kk_subj=1:nb_subj % ogni soggetto
    filename_Avalanches_subj=strcat(db_path,'/Classification/3_Dataset-netbci-meg-sess4-DK/Avalanches_Analysis_SVM_ds_subj_',num2str(kk_subj),'.mat');
    filename_Avalanches_subj_file= strcat('Avalanches_Analysis_SVM_ds_subj_',subject_IDs{kk_subj});

    disp(strcat('%%%%%%%%%%%%%%%%%% Individual stat analysis - MEG - DK - SVM -',subject_IDs{kk_subj}, ' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
    % individual-level
    Statistical_IndividualAnalysis_SVM_BH(filename_Avalanches_subj, filename_Avalanches_subj_file, fig_path, db_path, numperm, nb_ROIs_DK, p_thresh_stat, cortex_15002V_MNI, idx_DK);
    % prepare materials to compute reliability:
    load(filename_Avalanches_subj);
    load(strcat(db_path,'Stat_Analysis_Indiv_SVM_',filename_Avalanches_subj_file,'_BH_p_thresh_',num2str(p_thresh_stat),'.mat'));
    matrix_signif_diff_SVM(:,:,kk_subj)=pval_diff_perm_corrected_BH; % analisi BH
    matrix_signif_abs_diff_SVM(:,:,kk_subj)=pval_abs_diff_perm_corrected_BH; % analisi BH
end

    group_generator = @(x)x(randperm(numel(x)));
    temp_rr=[];
    temp_rr_abs=[];
    % reliablity over the subjects - group level
    for kk1=1:nb_perm_rand
        for kk_subj=1:nb_subj
            rand_pos=group_generator(1:nb_ROIs_DK);
            newrand_abs(:,:,kk_subj)=matrix_signif_abs_diff_SVM(rand_pos,rand_pos,kk_subj);
            newrand(:,:,kk_subj)=matrix_signif_diff_SVM(rand_pos,rand_pos,kk_subj);

        end
        rr_abs(:,:,kk1)=sum(newrand_abs,3);
        rr(:,:,kk1)=sum(newrand,3);

        temp_rr(:,kk1)=sum(rr(:,:,kk1)<0.05);
        temp_rr_abs(:,kk1)=sum(rr_abs(:,:,kk1)<0.05);
    end
    rr_thresh=mean(mean(temp_rr));
    rr_abs_thresh=mean(mean(temp_rr_abs));

    for kk_roi=1:nb_ROIs_DK
        temp_node_abs=sum(sum(matrix_signif_abs_diff_SVM,3),1);
        pval_node_abs_reliab_uncorrect(kk_roi)=length(find(temp_rr_abs(kk_roi,:)>temp_node_abs(kk_roi)))/nb_perm_rand;
        temp_node=sum(sum(matrix_signif_diff_SVM,3),1);
        pval_node_reliab_uncorrect(kk_roi)=length(find(temp_rr(kk_roi,:)>temp_node(kk_roi)))/nb_perm_rand;
    end

    % correction - BH
    pval_node_abs_reliab_correct=pval_adjust(pval_node_abs_reliab_uncorrect,'BH');
    pval_node_reliab_correct=pval_adjust(pval_node_reliab_uncorrect,'BH');


% edges
test_abs=sum(matrix_signif_abs_diff_SVM,3);
test2_abs=repmat(test_abs, 1, 1, nb_perm_rand);
abs_cmp=rr_abs>test2_abs; 
signif_abs_cmp=pval_adjust(sum(abs_cmp,3)./numperm,'BH')<0.05; % BH


test=sum(matrix_signif_diff_SVM,3); % quanto spesso ho una stima di un dato arco
test2=repmat(test, 1, 1, nb_perm_rand);
cmp=rr>test2; 
signif_cmp=pval_adjust(sum(cmp,3)./numperm,'BH')<0.05; % BH


% nodes - reliability based on edges results
node_cmp=sum(signif_cmp,2);
node_abs_cmp=sum(signif_abs_cmp,2);



save(strcat(db_path,'Rebuttal_Corr_SVM_Reliable_Analysis_Edge_Based_BH_p_thresh_edges_', num2str(0.05),'_p_thresh_nodes_',num2str(0.05),'_Group_''.mat'),...
        'signif_abs_cmp',... % post-correction bh
        'signif_cmp',... % post-correction bh
        'node_cmp',... % node signif
        'node_abs_cmp',... % node signif
        'rr_thresh',... % thresh for the nodes selection
        'rr_abs_thresh',...
        'pval_node_abs_reliab_uncorrect',...
        'pval_node_reliab_uncorrect',...
        'pval_node_abs_reliab_correct',...
        'pval_node_reliab_correct',...
        '-v7.3');

close all;
figure(1);
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,signif_cmp)
title('Significant diff - BH corrected')
axis square
filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_diff_SVM_Sess4_MEG_DK_p_thresh_', num2str(0.05),'.png');
print('-f1', filename, '-dpng','-r600');
figure(2)
imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,signif_abs_cmp)
title('Significant abs(diff) - BH corrected')
axis square
filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_abs_diff_SVM_Sess4_MEG_DK_p_thresh_', num2str(0.05), '.png');
print('-f2', filename, '-dpng','-r600');

DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, node_abs_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_abs_diff_SVM_Sess4_MEG_DK_p_thresh_', num2str(0.05)), 'Sess4_MEG_DK')
DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, node_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_diff_SVM_Sess4_MEG_DK_p_thresh_', num2str(0.05)), 'Sess4_MEG_DK')
