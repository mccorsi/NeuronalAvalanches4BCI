function Compute_Avalanches_EEG_Destrieux(path_file_MI, path_file_Rest, subject_IDs, db_path, labels_Destrieux, freq_band_def, z_thresh, fs, min_size_aval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Avalanches x dataset - 20 subj. - Sess. 4 - EEG - Desikan-Killiany
    % Authors: MCC
    % Date: 27/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

%% load data - can require time (>1Go for each file)...
load(path_file_MI);
load(path_file_Rest);

%% parameters definition
nb_subj=size(Data_EEG_Baseline_Destrieux, 1); % number of subjects
nb_ROIs=size(Data_EEG_Baseline_Destrieux{1,1}, 1); % to check the number of ROIs
nbsamples =size(Data_EEG_Baseline_Destrieux{1,1}, 2); % number of samples
ROIs = [1:nb_ROIs]; 
thresh = z_thresh;
temp = zeros(nb_ROIs,nbsamples);

%% Avalanches
for kk_subj=1:nb_subj % each subject
    %% define dataset
    close all; clc;
    notempty_MI_idx=find(~cellfun(@isempty,{Data_EEG_MI_Destrieux{kk_subj,:}}));
    notempty_Baseline_idx=find(~cellfun(@isempty,{Data_EEG_Baseline_Destrieux{kk_subj,:}}));
    nb_trials_tot=size(notempty_MI_idx,2) + size(notempty_Baseline_idx,2);

    % 0 x Baseline
    Dataset_EEG_Baseline_Destrieux.data={Data_EEG_Baseline_Destrieux{kk_subj,1:size(notempty_Baseline_idx,2)}};
    Dataset_EEG_Baseline_Destrieux.condition(1:size({Data_EEG_Baseline_Destrieux{kk_subj,1:size(notempty_Baseline_idx,2)}},2))=zeros(1,size({Data_EEG_Baseline_Destrieux{kk_subj,1:size(notempty_Baseline_idx,2)}},2));

    % 1 x MI
    Dataset_EEG_MI_Destrieux.data={Data_EEG_MI_Destrieux{kk_subj,1:size(notempty_MI_idx,2) }};
    Dataset_EEG_MI_Destrieux.condition(1:size({Data_EEG_MI_Destrieux{kk_subj,1:size(notempty_MI_idx,2) }},2))=ones(1,size({Data_EEG_MI_Destrieux{kk_subj,1:size(notempty_MI_idx,2) }},2)); 

    Dataset_EEG_Destrieux.data=[Dataset_EEG_MI_Destrieux.data,Dataset_EEG_Baseline_Destrieux.data];
    Dataset_EEG_Destrieux.condition=[Dataset_EEG_MI_Destrieux.condition,Dataset_EEG_Baseline_Destrieux.condition];

    %% "rasterplot" for each trial
    % concat MI & Baseline to compute z-scores
    data=[Dataset_EEG_Destrieux.data{:}];
    c=[1:nbsamples:(nb_trials_tot+1)*nbsamples];
    
    for kk_freq=1:size(freq_band_def,1) % each frequency band
        for kk_zthresh=1:length(z_thresh) % each z-score threshold
            for kk_min_size_aval=1:size(min_size_aval)
                
                disp(strcat('%%%%%%%%%%%%%%%%%% EEG - Destrieux - ',subject_IDs{kk_subj}, ' - ', num2str(freq_band_def{kk_freq}(1)), '-', num2str(freq_band_def{kk_freq}(2)),'Hz ',' - z_thresh:', num2str(z_thresh(kk_zthresh)), ' %%%%%%%%%%%%%%%%%%'));

                h_freq=freq_band_def{kk_freq}(2); % 40Hz
                l_freq=freq_band_def{kk_freq}(1); % 3 Hz
                [b,a]=butter(2,[l_freq h_freq]/(fs/2)); 

                filtered_series=cell(1,size(data,2));
                w=cell(1,size(data,2));
                filtered_series=filter(b,a,data,[],2); 
                w=zscore(filtered_series,0,2);
                activations=abs(w)>z_thresh(kk_zthresh);

                activations_trials=cell(1,nb_trials_tot);%nb_trials
                for  kk_trials=2:nb_trials_tot+1 %1:size(c,2)-1
                    activations_trials{kk_trials-1}=activations(:,c(kk_trials-1) : c(kk_trials)-1);

                 %% binning
                    % varies the binning of the avalanches 
                    n_binning=2:3;
    %                 aval_binned=cell(1,size(n_binning,2));
                    aval_nobin=activations_trials{kk_trials-1};

                    % binns the data
                    for kk_binning= 2:3 %size(n_binning,2) %loops across binnings
                       div=mod(size(aval_nobin,2),n_binning(kk_binning-1)); %is the size divisible by the binning
                       if div~=0
                           aval_nobin= aval_nobin(:,1:end-div); 
                       end
                       for kk5=1:size(aval_nobin,1) %scrolls across raws of the avalanche
                             temp=buffer(aval_nobin(kk5,:),n_binning(kk_binning-1));
                             aval_binned{kk_binning-1}{kk_trials-1}(kk5,:)=any(temp,1); 
                       end
                    end


                end
                    aval_binned_2(2:3) = aval_binned;
                    aval_binned_2{1} = activations_trials;


                    for kk_binning = 1:size(aval_binned_2,2)

                         for kk_trials = 1:nb_trials_tot
                            [avalanches{kk_binning}{kk_trials},patterns_aval{kk_trials}]= avalanches_global_pattern(aval_binned_2{kk_binning}{kk_trials},nb_ROIs);

                         end

                         for kk_trials=1:nb_trials_tot 
                            kk1=1;
                            out=[];
                            for kk2=1:size(avalanches{kk_binning}{kk_trials},2)
                                   if size(avalanches{kk_binning}{kk_trials}{kk2},2)>min_size_aval(kk_min_size_aval)
                                       out(1:nb_ROIs,1:nb_ROIs,kk1)=func_transition_matrix(avalanches{kk_binning}{kk_trials}{kk2});
                                       kk1=kk1+1;
                                   end    
                            end
                        mask=sum(out~=0,3)>(0.20*size(out,3)); 
                        out_temp2=zeros(size(out));
                        out_temp2(repmat(mask,1,1,size(out,3)))=out(repmat(mask,1,1,size(out,3)));

                        out_temp=sum(out_temp2,3)./sum(out_temp2~=0,3);
                        % to make the matrix symmetrical
                        out_temp(isnan(out_temp))=0;
                        TM_out{kk_binning}{kk_trials}=(out_temp+out_temp')./2;

                        clearvars out

                end
                val_empty.values=cellfun(@isempty,TM_out{kk_binning});
                val_empty.conditions=Dataset_EEG_Destrieux.condition;

                nb_trials_tot=size(notempty_MI_idx,2) + size(notempty_Baseline_idx,2);

                out_MI={TM_out{kk_binning}{1:size(notempty_MI_idx,2)}};
                TM_MI{kk_binning}={mean(cat(3,out_MI{:}),3)};

                out_Baseline={TM_out{kk_binning}{size(notempty_MI_idx,2)+1:nb_trials_tot}};
                TM_Baseline{kk_binning}={mean(cat(3,out_Baseline{:}),3)};

                TM_Diff{kk_binning}={(TM_MI{kk_binning}{:}-TM_Baseline{kk_binning}{:})};
                TM_abs_Diff{kk_binning}={abs(TM_Diff{kk_binning}{:})};
      end
           % save results in a single .mat file

            save(strcat(db_path,'Avalanches_Analysis_ds_',subject_IDs{kk_subj} ,'_',num2str(freq_band_def{kk_freq}(1)), '_', num2str(freq_band_def{kk_freq}(2)),'_z_', num2str(z_thresh(kk_zthresh)), '_min_size_aval_',num2str(min_size_aval(kk_min_size_aval)),'_Sess4_EEG_Destrieux.mat'),...
                'Dataset_EEG_Destrieux',...
                'TM_out',...
                'TM_abs_Diff',...
                'TM_Diff',...
                'val_empty',...
                'labels_Destrieux',...
                'notempty_MI_idx',...
                'notempty_Baseline_idx',...
                'avalanches',...
                '-v7.3');

        %     clearvars Dataset_EEG_Destrieux Dataset_EEG_Baseline_Destrieux Dataset_EEG_MI_Destrieux avalanches TM_out
            end
end
            

            
            
    end
 end
end
