%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN - PWelch analysis - 20 subjects - MEG - Desikan-Killiany
    % Authors: MCC
    % Date: 25/10/2022
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
freq_band_def={[5, 7],[8 12],[13,30],[3 40]};

%% Load data
path_file_MI=strcat(db_path,'Data_PWelch_MEG_MI_DK.mat');
path_file_Rest=strcat(db_path,'Data_PWelch_MEG_Baseline_DK.mat');

load(path_file_MI)
load(path_file_Rest)

%% Statistical Analysis - Individual
clc
numperm=10000;
p_thresh_stat=0.05; % individual stat
nb_perm_rand=10000;%200

for kk_subj=1:nb_subj % ogni soggetto
    subject_ID=subject_IDs{kk_subj};
    PWelch_Baseline_Indiv={Data_MEG_Baseline_DK{kk_subj,:}};
    PWelch_MI_Indiv={Data_MEG_MI_DK{kk_subj,:}};
    nb_trials_MI=size(find(~cellfun(@isempty,PWelch_MI_Indiv)),2);
    nb_trials_baseline=size(find(~cellfun(@isempty,PWelch_Baseline_Indiv)),2);
    nb_trials=min(nb_trials_MI,nb_trials_baseline);
    disp(strcat('%%%%%%%%%%%%%%%%%% Individual stat analysis - MEG - DK - ',subject_IDs{kk_subj}, '%%%%%%%%%%%%%%%%%%'));


    for kk_freq=1:size(freq_band_def,2) % ogni banda di frequenza considerata

        disp(strcat('%%%%%%%%%%%%%%%%%%', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
        % individual-level

        for kk_trials=1:nb_trials
            PWelch_Baseline_temp{kk_trials}=squeeze(PWelch_Baseline_Indiv{1,kk_trials}(:,1,kk_freq));
            PWelch_MI_temp{kk_trials}=squeeze(PWelch_MI_Indiv{1,kk_trials}(:,1,kk_freq));
        end
            Data_MEG_MI_DK_indiv=reshape([PWelch_MI_temp{:}], nb_ROIs_DK,  nb_trials); % nb_rois x nb_trials
            Data_MEG_Baseline_DK_indiv=reshape([PWelch_Baseline_temp{:}], nb_ROIs_DK,nb_trials); % nb_rois x nb_trials

        Data_MEG_MI_DK_indiv=cell2mat(Data_MEG_MI_DK(kk_subj,:));% nb_rois x nb_trials
        Data_MEG_Baseline_DK_indiv=cell2mat(Data_MEG_Baseline_DK(kk_subj,:));% nb_rois x nb_trials

        Statistical_IndividualAnalysis_PWelch_BH(subject_ID,Data_MEG_MI_DK_indiv, Data_MEG_Baseline_DK_indiv, db_path, numperm, nb_ROIs_DK, p_thresh_stat, freq_band_def{kk_freq});
    end
    clearvars subject_ID Data_MEG_Baseline_DK_indiv Data_MEG_MI_DK_indiv nb_trials_MI nb_trials_baseline nb_trials PWelch_Baseline_temp PWelch_MI_temp PWelch_MI PWelch_Baseline
end
%%
clc
for kk_freq=1:size(freq_band_def,2)

    for kk_subj=1:nb_subj
        subject_ID=subject_IDs{kk_subj};

        % prepare materials to compute reliability:
        load(strcat(db_path,'Stat_Analysis_Indiv_PWelch_',subject_ID,'freq','_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_BH_p_thresh_',num2str(p_thresh_stat),'.mat'));
        matrix_signif_diff(:,kk_subj)=pval_node_diff_perm_corrected_BH; % analisi BH
        clearvars pval_node_diff_perm_corrected_BH
    end
    close all
    % plot signif matrices - indiv
    figure()
    imagesc(matrix_signif_diff); axis square
    xlabel('Subjects')
    ylabel('ROIs')
    colorbar
    filename=strcat(fig_path,'SignifMatrix_PWlech_BH_diff_Sess4_MEG_DK_p_thresh_', num2str(p_thresh_stat),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas, '.png');
    print('-f1', filename, '-dpng','-r600');

     disp(strcat('%%%%%%%%%%%%%%%%%% Group stat analysis - MEG - DK - ', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));

    group_generator = @(x)x(randperm(numel(x)));
    temp_rr=[];
    % reliablity over the subjects - group level
    for kk1=1:nb_perm_rand

        rand_pos=group_generator(1:nb_ROIs_DK);
        rand_sog=group_generator(1:nb_subj);
        newrand=matrix_signif_diff(rand_pos,rand_sog);
        rr(:,kk1)=sum(newrand,2);

    end

    test=sum(matrix_signif_diff,2); % quanto spesso ho una stima di un dato arco
    test2=repmat(test, 1, nb_perm_rand);
    cmp=rr>=test2;
    pval_pwelch_reliab_uncorrect=sum(cmp,2)./numperm;
    pval_pwelch_reliab_correct=pval_adjust(pval_pwelch_reliab_uncorrect,'BH')<0.05; % BH


    save(strcat(db_path,'Reliable_Analysis_PWelch_Based_BH_p_thresh_edges_', num2str(0.05),'_p_thresh_nodes_',num2str(0.05),'_Group_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas,'.mat'),...
            'pval_pwelch_reliab_uncorrect',...
            'pval_pwelch_reliab_correct',...
            '-v7.3');

   DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, pval_pwelch_reliab_correct, strcat(fig_path,'GroupAnalysis_Reliable_PWelch_BH_diff_Sess4_MEG_DK_p_thresh_', num2str(0.05),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas), 'Sess4_MEG_DK')

end
