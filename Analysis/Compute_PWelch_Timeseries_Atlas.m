%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % "Batch" code to compute PWelch & PLV after spatial downsampling - 20 subjects - MEG & DK only
    % Authors: MCC
    % Date: 25/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;
%% to be adapted to each user:
root_path='/Users/marieconstance.corsi/Documents/GitHub/NeuronalAvalanches4BCI/';
db_path=strcat(root_path,'Database/');
fig_path=strcat(root_path,'Figures/');

cd(strcat(root_path,'/Analysis'))
Start_SessionMatlab(root_path)
%% lista degli progetti Brainstorm da usare
nb_subj=20;
nb_runs=6;
nb_trials=16; % per condizione & per run
conditions_IDs={'MI','Baseline'};
atlas_IDs={'Desikan-Killiany'};
nb_atlas=size(atlas_IDs,2);
freq_band_def={[5, 7],[8 12],[13,30],[3 40]};
fs=250;

info_sess='';
subject_IDs={'',... % list of subjects_IDs
    };


root_bst_db='';
for kk_subj=1:nb_subj
    Database_bst.Protocol_IDs{kk_subj}=strcat(subject_IDs{kk_subj},info_sess);
end

% load la lista degli files da usare per ogni soggetto --> sFilesGroup
List_Files_Power

%%
% downsampling
Data_MEG_MI_DK=cell(nb_subj,96); % nb subj x nb max of trials
Data_MEG_Baseline_DK=cell(nb_subj,96);
Timeseries_MEG_MI_DK=cell(nb_subj,96); % nb subj x nb max of trials
Timeseries_MEG_Baseline_DK=cell(nb_subj,96);

for kk_subj=1:nb_subj
    ProtocolName = Database_bst.Protocol_IDs{kk_subj};
    iProtocol = bst_get('Protocol', ProtocolName);
    if isempty(iProtocol)
        error(['Unknown protocol: ' ProtocolName]);
    end

    % seleziona il protocollo
    gui_brainstorm('SetCurrentProtocol', iProtocol);

        idx_MEG_MI_DK=1;
        idx_MEG_Baseline_DK=1;

        % Input files
        sFiles = sFilesGroup{kk_subj};


%         %% extract timeseries used afterwards to compute PLV
%             % Start a new report
%             bst_report('Start', sFiles);
%
%             % Process: Extract values: [0.000s,6.996s]
%             sFiles_out_extract = bst_process('CallProcess', 'process_extract_time', sFiles, [], ...
%                 'timewindow', [0, 6.996], ...
%                 'overwrite',  0);
%
%             % Save and display report
%             ReportFile = bst_report('Save', sFiles_out_extract);
%             bst_report('Open', ReportFile);
%
%             path_project=strcat(root_bst_db,subject_IDs{kk_subj},'_RebootMRI/',ProtocolName,'/data/');
%             for kk_files= 1:size(sFiles_out_extract,2)
%                 path_file=strcat(path_project,sFiles_out_extract(kk_files).FileName);
%                 if contains(path_file,'_time.mat')
%                     temp=load(strcat(path_project,sFiles_out_extract(kk_files).FileName));
%
%                     if contains(path_file, 'MI') && contains(temp.Comment, 'MEG') && strcmp(temp.Atlas.Name,'Desikan-Killiany')
%                         Timeseries_MEG_MI_DK{kk_subj,idx_MEG_MI_DK}=temp.ImageGridAmp;
%                         PLV_MI_matrix{kk_subj,idx_MEG_MI_DK}=DoMyPLV_MEG_DK(freq_band_def,Timeseries_MEG_MI_DK{kk_subj,idx_MEG_MI_DK},fs);
%                         idx_MEG_MI_DK=idx_MEG_MI_DK+1;
%                     elseif contains(path_file, 'Baseline') && contains(temp.Comment, 'MEG') && strcmp(temp.Atlas.Name,'Desikan-Killiany')
%                         Timeseries_Baseline_DK{kk_subj,idx_MEG_Baseline_DK}=temp.ImageGridAmp;
%                         PLV_Baseline_matrix{kk_subj,idx_MEG_Baseline_DK}=DoMyPLV_MEG_DK(freq_band_def,Timeseries_Baseline_DK{kk_subj,idx_MEG_Baseline_DK},fs);
%                         idx_MEG_Baseline_DK=idx_MEG_Baseline_DK+1;
%
%                     end
%                 end
%             end

        %% compute pwelch
        bst_report('Start', sFiles);

        % Process: Power spectrum density (Welch)
        sFiles_out = bst_process('CallProcess', 'process_psd', sFiles, [], ...
            'timewindow',  [0, 6.996], ...
            'win_length',  1, ...
            'win_overlap', 50, ...
            'units',       'physical', ...  % Physical: U2/Hz
            'clusters',    {}, ...
            'scoutfunc',   1, ...  % Mean
            'win_std',     0, ...
            'edit',        struct(...
                 'Comment',         'Power,FreqBands', ...
                 'TimeBands',       [], ...
                 'Freqs',           {{'theta', '5, 7', 'mean'; 'alpha', '8, 12', 'mean'; 'beta', '13, 29', 'mean'; 'default', '3, 40', 'mean'}}, ... % sufficiente per noi qui
                 'ClusterFuncTime', 'none', ...
                 'Measure',         'power', ...
                 'Output',          'all', ...
                 'SaveKernel',      0));

        % Save and display report
        ReportFile = bst_report('Save', sFiles_out);
        bst_report('Open', ReportFile);

        % load results & gather them in a separate file


    path_project=strcat(root_bst_db,subject_IDs{kk_subj},'_RebootMRI/',ProtocolName,'/data/');
    for kk_files= 1:size(sFiles_out,2)
        temp=load(strcat(path_project,sFiles_out(kk_files).FileName));
        if contains(temp.DataFile, 'MI') && contains(temp.Comment, 'Power,FreqBands') && strcmp(temp.Atlas.Name,'Desikan-Killiany')
            Data_MEG_MI_DK{kk_subj,idx_MEG_MI_DK,:}=temp.TF;
            idx_MEG_MI_DK=idx_MEG_MI_DK+1;
        elseif contains(temp.DataFile, 'Baseline') && contains(temp.Comment, 'Power,FreqBands') && strcmp(temp.Atlas.Name,'Desikan-Killiany')
            Data_MEG_Baseline_DK{kk_subj,idx_MEG_Baseline_DK,:}=temp.TF;
            idx_MEG_Baseline_DK=idx_MEG_Baseline_DK+1;

        end
    end
end


%% save results
save(strcat(db_path,'Data_PWelch_MEG_MI_DK.mat'),...
    'Data_MEG_MI_DK', '-V7.3');
save(strcat(db_path,'Data_PWelch_MEG_Baseline_DK.mat'),...
    'Data_MEG_Baseline_DK', '-V7.3');

save(strcat(db_path, 'PLV_Timeseries_Baseline_20Subj_Sess4_MEG_DK.mat'),...
    'Timeseries_Baseline_DK',...
    'Timeseries_MI_DK',...
    'PLV_Baseline_matrix',...
    'PLV_MI_matrix',...
    'subject_IDs', ...
    '-v7.3');
