%% To be modified

root_path='/Users/marie-constance.corsi/Documents/GitHub/NeuronalAvalanches4BCI';
%% packages to be downloaded and added to the path
addpath(strcat(root_path,'Scripts/0_packages/'));
cd(strcat(root_path,'Scripts/0_packages/'))
addpath(genpath('geostd'));
addpath(genpath('BrewerMap-master'));
addpath(genpath('matplotlib'));
addpath(genpath('EEGLAB'));
addpath(genpath('redblue'));
addpath(genpath('pval_adjust-master'));
%% additional scripts for plots
addpath(strcat(root_path,'Scripts/Visualization_Cortex/'))
cd(strcat(root_path,'Scripts'));
