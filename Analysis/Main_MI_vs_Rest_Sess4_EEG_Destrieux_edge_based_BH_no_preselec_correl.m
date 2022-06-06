%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN - EEG - Destrieux
    % Authors: MCC
    % Date: 21/03/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% To be updated

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
idx_Destrieux=2;
nb_ROIs_Destrieux=size(cortex_15002V_MNI.Atlas(idx_Destrieux).Scouts,2);
labels_Destrieux={cortex_15002V_MNI.Atlas(idx_Destrieux).Scouts(1:nb_ROIs_Destrieux).Label};

modality='EEG';
atlas='Destrieux';
fs=250;


%% Parameters to be tested

% per default
freq_band_default={[3 40]};
z_thresh_default=[3]; 
min_size_aval_default=[5]; 
n_binning_default=[2]; %=1

% for the robustness test
freq_band_def=freq_band_default;
z_thresh=[2.5, 2.7, 3, 3.2]; 
min_size_aval=[1,5,10,15]; 
n_binning=1:3;

%% ATMs and avalanches computation - to be updated with proper filenames
path_file_MI=strcat(db_path,'Data_EEG_MI_Destrieux_V2.mat');
path_file_Rest=strcat(db_path,'Data_EEG_Baseline_Destrieux_V2.mat');

Compute_Avalanches_EEG_Destrieux(path_file_MI, path_file_Rest, subject_IDs, db_path, labels_Destrieux, freq_band_def, z_thresh, fs, min_size_aval)


%% Statistical Analysis - Individual
clc
numperm=10000;
list_p_thresh=[0.05;0.025;0.05];
p_thresh_stat=0.05; % individual stat
nb_perm_rand=10000;

for kk_freq=1:size(freq_band_def,1) % ogni banda di frequenza considerata
    for kk_zthresh=1:length(z_thresh) % ogni soglia x z-score 
        for kk_binning=1:length(n_binning) % loop sul binning
            for kk_min_size_aval=1:length(min_size_aval)
        %% TODO : add loop sul binning
        
            for kk_subj=1:nb_subj % ogni soggetto
                filename_Avalanches_subj=strcat(db_path,'Avalanches_Analysis_ds_',subject_IDs{kk_subj} ,'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh), '_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_EEG_Destrieux.mat');
                filename_Avalanches_subj_file= strcat('Avalanches_Analysis_ds_',subject_IDs{kk_subj} ,'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh), '_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_', modality, '_',atlas);

                disp(strcat('%%%%%%%%%%%%%%%%%% Individual stat analysis - EEG - Destrieux - ',subject_IDs{kk_subj}, ' - ', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - z_thresh:', num2str(z_thresh), ' - p_thresh:,' , num2str(p_thresh_stat),'%%%%%%%%%%%%%%%%%%'));
                % individual-level
                Statistical_IndividualAnalysis_BH(filename_Avalanches_subj, filename_Avalanches_subj_file, fig_path, db_path, numperm, nb_ROIs_Destrieux, p_thresh_stat, cortex_15002V_MNI, idx_Destrieux, n_binning(kk_binning));
                % prepare materials to compute reliability:
                load(filename_Avalanches_subj);
                load(strcat(db_path,'Stat_Analysis_Indiv_',filename_Avalanches_subj_file,'_BH_p_thresh_',num2str(p_thresh_stat),'.mat'));
                matrix_signif_diff(:,:,kk_subj)=pval_diff_perm_corrected_BH; % analisi BH
                matrix_signif_abs_diff(:,:,kk_subj)=pval_abs_diff_perm_corrected_BH; % analisi BH
            end

            group_generator = @(x)x(randperm(numel(x)));
            temp_rr=[];
            temp_rr_abs=[];
            % reliablity over the subjects - group level
            for kk1=1:nb_perm_rand
                for kk_subj=1:nb_subj
                    rand_pos=group_generator(1:nb_ROIs_Destrieux);
                    newrand_abs(:,:,kk_subj)=matrix_signif_abs_diff(rand_pos,rand_pos,kk_subj);
                    newrand(:,:,kk_subj)=matrix_signif_diff(rand_pos,rand_pos,kk_subj);
                    
                end
                rr_abs(:,:,kk1)=sum(newrand_abs,3);
                rr(:,:,kk1)=sum(newrand,3);
                
                temp_rr(:,kk1)=sum(rr(:,:,kk1)<0.05);
                temp_rr_abs(:,kk1)=sum(rr_abs(:,:,kk1)<0.05);
            end
            rr_thresh=mean(mean(temp_rr));
            rr_abs_thresh=mean(mean(temp_rr_abs));
            
            for kk_roi=1:nb_ROIs_Destrieux
                temp_node_abs=sum(sum(matrix_signif_abs_diff,3),1);
                pval_node_abs_reliab_uncorrect(kk_roi)=length(find(temp_rr_abs(kk_roi,:)>temp_node_abs(kk_roi)))/nb_perm_rand;
                temp_node=sum(sum(matrix_signif_diff,3),1);
                pval_node_reliab_uncorrect(kk_roi)=length(find(temp_rr(kk_roi,:)>temp_node(kk_roi)))/nb_perm_rand;
            end
            
            
            % correction - BH
            pval_node_abs_reliab_correct=pval_adjust(pval_node_abs_reliab_uncorrect,'BH');
            pval_node_reliab_correct=pval_adjust(pval_node_reliab_uncorrect,'BH');
            
%             % find nodes where test_corrected<0.05
%             list_nodes_abs_reliab=find(pval_node_abs_reliab_correct<0.05);
%             list_nodes_reliab=find(pval_node_reliab_correct<0.05);
            
            
        % edges
        test_abs=sum(matrix_signif_abs_diff,3);
        test2_abs=repmat(test_abs, 1, 1, nb_perm_rand);
        abs_cmp=rr_abs>test2_abs; 
        signif_abs_cmp=pval_adjust(sum(abs_cmp,3)./numperm,'BH')<0.05; % BH
%         signif_abs_cmp = signif_abs_cmp - diag(diag(signif_abs_cmp));
%         signif_abs_cmp_bhfdr=reshape(mafdr(reshape(sum(abs_cmp,3)./numperm,nb_ROIs_Destrieux*nb_ROIs_Destrieux,1),'BHFDR','true'),nb_ROIs_Destrieux,nb_ROIs_Destrieux);

%         
%         mask=tril(ones(nb_ROIs_Destrieux),-1);
%         temp=sum(abs_cmp,3)./numperm;
%         p_vect=temp(logical(mask));
%         signif_abs_cmp_bhfdr=mafdr(p_vect,'BHFDR','true');
        
        
        
        
        test=sum(matrix_signif_diff,3); % quanto spesso ho una stima di un dato arco
        test2=repmat(test, 1, 1, nb_perm_rand);
        cmp=rr>test2; 
        signif_cmp=pval_adjust(sum(cmp,3)./numperm,'BH')<0.05; % BH

        % nodes - reliability based on edges results
        node_cmp=sum(signif_cmp,2); 
        node_abs_cmp=sum(signif_abs_cmp,2); 

        
        
        
        save(strcat(db_path,'Reliable_Analysis_Edge_Based_BH_p_thresh_edges_', num2str(0.05),'_p_thresh_nodes_',num2str(0.05),'_Group_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh), '_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_', modality, '_',atlas,'.mat'),...
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
        imagesc(1:nb_ROIs_Destrieux,1:nb_ROIs_Destrieux,signif_cmp)
        title('Significant diff - BH corrected')
        axis square
        filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_diff_Sess4_EEG_Destrieux_p_thresh_', num2str(0.05),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh),'_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)), '_Sess4_', modality, '_',atlas, '.png');
        print('-f1', filename, '-dpng','-r600');
        figure(2)
        imagesc(1:nb_ROIs_Destrieux,1:nb_ROIs_Destrieux,signif_abs_cmp)
        title('Significant abs(diff) - BH corrected')
        axis square
        filename=strcat(fig_path,'GroupAnalysis_Edge_Based_Reliable_Edges_BH_abs_diff_Sess4_EEG_Destrieux_p_thresh_', num2str(0.05), '_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh),'_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)), '_Sess4_', modality, '_',atlas, '.png');
        print('-f2', filename, '-dpng','-r600');

        DoMyViz_node_rel(cortex_15002V_MNI, idx_Destrieux, node_abs_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_abs_diff_Sess4_EEG_Destrieux_p_thresh_', num2str(0.05),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh), '_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_', modality, '_',atlas), 'Sess4_EEG_Destrieux')
        DoMyViz_node_rel(cortex_15002V_MNI, idx_Destrieux, node_cmp, strcat(fig_path,'GroupAnalysis_Reliable_Edge_Based_Nodes_BH_diff_Sess4_EEG_Destrieux_p_thresh_', num2str(0.05),'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh), '_binning_',num2str(n_binning(kk_binning)),'_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_', modality, '_',atlas), 'Sess4_EEG_Destrieux')
            end
        end
   end
end
