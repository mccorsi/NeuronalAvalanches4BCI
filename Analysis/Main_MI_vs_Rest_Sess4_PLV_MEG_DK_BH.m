%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN - Analysi fatta sul dataset a disposizione - PLV - 20 soggetti - Sess. 4 - MEG - Desikan-Killiany
    % Authors: MCC
    % Date: 26/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% Da modificare x ogni persona

% root_path='/Users/marie-constance.corsi/Documents/GitHub/Fenicotteri-equilibristi/';
root_path='/Users/marieconstance.corsi/Documents/GitHub/Fenicotteri-equilibristi/';

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

%% Load PLV matrices
load(strcat(db_path,'PLV_Timeseries_Baseline_20Subj_Sess4_MEG_DK.mat'));

%% Statistical Analysis - Individual
clc
numperm=10000;%200
p_thresh_stat=0.05; % individual stat
nb_perm_rand=10000;%200

for kk_subj=1:nb_subj
    subject_ID=subject_IDs{kk_subj};
    PLV_Baseline_Indiv={PLV_Baseline_matrix{kk_subj,:}};
    PLV_MI_Indiv={PLV_MI_matrix{kk_subj,:}};
    nb_trials_MI=size(find(~cellfun(@isempty,PLV_MI_Indiv)),2);
    nb_trials_baseline=size(find(~cellfun(@isempty,PLV_Baseline_Indiv)),2);
    nb_trials=min(nb_trials_MI,nb_trials_baseline);
    disp(strcat('%%%%%%%%%%%%%%%%%% Individual stat analysis - MEG - DK - ',subject_IDs{kk_subj}, '%%%%%%%%%%%%%%%%%%'));
    for kk_freq=1:size(freq_band_def,2)
        for kk_trials=1:nb_trials
            PLV_Baseline_temp{kk_trials}=PLV_Baseline_Indiv{1,kk_trials}{1,kk_freq};
            PLV_MI_temp{kk_trials}=PLV_MI_Indiv{1,kk_trials}{1,kk_freq};
        end
            PLV_MI=reshape([PLV_MI_temp{:}], nb_ROIs_DK, nb_ROIs_DK, nb_trials);
            PLV_Baseline=reshape([PLV_Baseline_temp{:}], nb_ROIs_DK, nb_ROIs_DK, nb_trials);

        disp(strcat('%%%%%%%%%%%%%%%%%%', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
        % individual-level
        Statistical_IndividualAnalysis_PLV_BH(PLV_MI, PLV_Baseline, subject_ID, fig_path, db_path, numperm, nb_ROIs_DK, p_thresh_stat,freq_band_def{kk_freq});
    end
    clearvars subject_ID PLV_Baseline_Indiv PLV_MI_Indiv nb_trials_MI nb_trials_baseline nb_trials PLV_Baseline_temp PLV_MI_temp PLV_MI PLV_Baseline
end
%%
clc
for kk_freq=1:size(freq_band_def,2)

    for kk_subj=1:nb_subj
        subject_ID=subject_IDs{kk_subj};
        % prepare materials to compute reliability:
        load(strcat(db_path,'Stat_Analysis_Indiv_PLV_',subject_ID,'freq','_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_BH_p_thresh_',num2str(p_thresh_stat),'.mat'));
        matrix_signif_diff(:,:,kk_subj)=pval_diff_perm_corrected_BH;


    end
     disp(strcat('%%%%%%%%%%%%%%%%%% Group stat analysis - MEG - DK - ', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - p_thresh:' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
            group_generator = @(x)x(randperm(numel(x)));
            temp_rr=[];
            temp_rr_abs=[];
            % reliablity over the subjects - group level
            for kk1=1:nb_perm_rand
                rand_pos=group_generator(1:nb_ROIs_DK);
                rand_subj=group_generator(1:nb_subj);
                newrand=matrix_signif_diff(rand_pos,rand_pos,rand_subj);
                rr(:,:,kk1)=sum(newrand,3);

                temp_rr(:,kk1)=sum(rr(:,:,kk1)<0.05);
            end
            rr_thresh=mean(mean(temp_rr));

            for kk_roi=1:nb_ROIs_DK
                temp_node=sum(sum(matrix_signif_diff,3),1);
                pval_node_reliab_uncorrect(kk_roi)=length(find(temp_rr(kk_roi,:)>temp_node(kk_roi)))/nb_perm_rand;
            end

            % correction - BH
            pval_node_reliab_correct=pval_adjust(pval_node_reliab_uncorrect,'BH');


        % edges
        test=sum(matrix_signif_diff,3); % quanto spesso ho una stima di un dato arco
        test2=repmat(test, 1, 1, nb_perm_rand);
        cmp=rr>=test2;
        signif_cmp=pval_adjust(sum(cmp,3)./numperm,'BH')<0.05; % BH

        % nodes - reliability based on edges results
        node_cmp=sum(signif_cmp,2); % nodi dov'? signif


        save(strcat(db_path,'Reliable_Analysis_PLV_Edge_Based_BH_p_thresh_edges_', num2str(p_thresh_stat),'_p_thresh_nodes_',num2str(p_thresh_stat),'_Group_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas,'.mat'),...
                'signif_cmp',... % post-correction bh
                'node_cmp',... % node signif
                'rr_thresh',... % thresh for the nodes selection
                'pval_node_reliab_uncorrect',...
                'pval_node_reliab_correct',...
                '-v7.3');

        close all;
        figure();
        imagesc(1:nb_ROIs_DK,1:nb_ROIs_DK,signif_cmp)
        title('Significant diff - BH corrected')
        axis square
        filename=strcat(fig_path,'GroupAnalysis_PLV_Edge_Based_Reliable_Edges_BH_diff_Sess4_MEG_DK_p_thresh_', num2str(p_thresh_stat),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas, '.png');
        print('-f1', filename, '-dpng','-r600');

        DoMyViz_node_rel(cortex_15002V_MNI, idx_DK, node_cmp, strcat(fig_path,'GroupAnalysis_Reliable_PLV_Edge_Based_Nodes_BH_diff_Sess4_MEG_DK_p_thresh_', num2str(p_thresh_stat),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_Sess4_', modality, '_',atlas), 'Sess4_MEG_DK')

end
