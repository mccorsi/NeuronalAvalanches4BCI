
# create a script from a matrix and a vector of labels - plot connectome

import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
import mne
import numpy as np
#%%
## functions
def isSquare (m): return all (len (row) == len (m) for row in m)

def regions_bst_to_lobes( argument ):
    '''
        gather brainstorm labels in 5 lobes/areas
    '''
    switcher = {
        # gathers Left and Right in one single lobe
        # Frontal: includes both Frontal and PreFrontal in this way
        'LF': 'F' ,
        'RF': 'F' ,
        'LPF': 'F' ,
        'RPF': 'F' ,
        # Central: C --> Motor: M, except postcentral that is going to be placed in parietal...
        'LC': 'M' ,
        'RC': 'M' ,
        'M': 'M', # in case that already done
        # Temporal: includes both Temporal and Lateral
        'LT': 'T' ,
        'RT': 'T' ,
        'LL': 'T' ,
        'RL': 'T' ,
        # Parietal
        'LP': 'P' ,
        'RP': 'P' ,
        'P': 'P',  # in case that already done
        # Occipital
        'LO': 'O' ,
        'RO': 'O'
    }
    return switcher.get ( argument )

def Lobes_Partition( MatrixToBeDivided , idx_F , idx_M , idx_P , idx_T , idx_O ):
    SubMatrix_F_F = MatrixToBeDivided[ np.ix_ ( idx_F , idx_F ) ]
    SubMatrix_F_M = MatrixToBeDivided[ np.ix_ ( idx_F , idx_M ) ]
    SubMatrix_F_T = MatrixToBeDivided[ np.ix_ ( idx_F , idx_T ) ]
    SubMatrix_F_P = MatrixToBeDivided[ np.ix_ ( idx_F , idx_P ) ]
    SubMatrix_F_O = MatrixToBeDivided[ np.ix_ ( idx_F , idx_O ) ]

    SubMatrix_M_M = MatrixToBeDivided[ np.ix_ ( idx_M , idx_M ) ]
    SubMatrix_M_T = MatrixToBeDivided[ np.ix_ ( idx_M , idx_T ) ]
    SubMatrix_M_P = MatrixToBeDivided[ np.ix_ ( idx_M , idx_P ) ]
    SubMatrix_M_O = MatrixToBeDivided[ np.ix_ ( idx_M , idx_O ) ]

    SubMatrix_T_T = MatrixToBeDivided[ np.ix_ ( idx_T , idx_T ) ]
    SubMatrix_T_P = MatrixToBeDivided[ np.ix_ ( idx_T , idx_P ) ]
    SubMatrix_T_O = MatrixToBeDivided[ np.ix_ ( idx_T , idx_O ) ]

    SubMatrix_P_P = MatrixToBeDivided[ np.ix_ ( idx_P , idx_P ) ]
    SubMatrix_P_O = MatrixToBeDivided[ np.ix_ ( idx_P , idx_O ) ]

    SubMatrix_O_O = MatrixToBeDivided[ np.ix_ ( idx_O , idx_O ) ]

    return SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O


#%%
#TODO: update the root_path!
root_path = '/Users/marie-constance.corsi/Documents/GitHub/NeuronalAvalanches4BCI'
db_path = root_path + '/Database/'
fig_path = root_path + '/Figures/'

# Desikan-Kiliany
ROI_DK_list = [ 'bankssts L' , 'bankssts R' , 'caudalanteriorcingulate L' , 'caudalanteriorcingulate R' ,
                'caudalmiddlefrontal L' , 'caudalmiddlefrontal R' , 'cuneus L' , 'cuneus R' , 'entorhinal L' ,
                'entorhinal R' , 'frontalpole L' , 'frontalpole R' , 'fusiform L' , 'fusiform R' ,
                'inferiorparietal L' , 'inferiorparietal R' , 'inferiortemporal L' , 'inferiortemporal R' , 'insula L' ,
                'insula R' , 'isthmuscingulate L' , 'isthmuscingulate R' , 'lateraloccipital L' , 'lateraloccipital R' ,
                'lateralorbitofrontal L' , 'lateralorbitofrontal R' , 'lingual L' , 'lingual R' ,
                'medialorbitofrontal L' , 'medialorbitofrontal R' , 'middletemporal L' , 'middletemporal R' ,
                'paracentral L' , 'paracentral R' , 'parahippocampal L' , 'parahippocampal R' , 'parsopercularis L' ,
                'parsopercularis R' , 'parsorbitalis L' , 'parsorbitalis R' , 'parstriangularis L' ,
                'parstriangularis R' , 'pericalcarine L' , 'pericalcarine R' , 'postcentral L' , 'postcentral R' ,
                'posteriorcingulate L' , 'posteriorcingulate R' , 'precentral L' , 'precentral R' , 'precuneus L' ,
                'precuneus R' , 'rostralanteriorcingulate L' , 'rostralanteriorcingulate R' , 'rostralmiddlefrontal L' ,
                'rostralmiddlefrontal R' , 'superiorfrontal L' , 'superiorfrontal R' , 'superiorparietal L' ,
                'superiorparietal R' , 'superiortemporal L' , 'superiortemporal R' , 'supramarginal L' ,
                'supramarginal R' , 'temporalpole L' , 'temporalpole R' , 'transversetemporal L' ,
                'transversetemporal R' ]
ROIs = '%0d'.join ( ROI_DK_list )
Region_DK = [ 'LT' , 'RT' , 'LL' , 'RL' , 'LF' , 'RF' , 'LO' , 'RO' , 'LT' , 'RT' , 'LPF' , 'RPF' , 'LT' , 'RT' , 'LP' ,
              'RP' , 'LT' , 'RT' , 'LT' , 'RT' , 'LL' , 'RL' , 'LO' , 'RO' , 'LPF' , 'RPF' , 'LO' , 'RO' , 'LPF' ,
              'RPF' , 'LT' , 'RT' , 'LC' , 'RC' , 'LT' , 'RT' , 'LF' , 'RF' , 'LPF' , 'RPF' , 'LF' , 'RF' , 'LO' , 'RO' , 'LC' ,
              'RC' , 'LL' , 'RL' , 'LC' , 'RC' , 'LP' , 'RP' , 'LL' , 'RL' , 'LF' , 'RF' , 'LF' , 'RF' , 'LP' , 'RP' ,
              'LT' , 'RT' , 'LP' , 'RP' , 'LT' , 'RT' , 'LT' , 'RT' ] # partition provided by Brainstorm


idx_DK_postcent=[ i for i , j in enumerate ( ROI_DK_list ) if j == 'postcentral L' or  j == 'postcentral R']
idx_DK_Maudalmiddlefront=[ i for i , j in enumerate ( ROI_DK_list ) if j == 'caudalmiddlefrontal L' or  j == 'caudalmiddlefrontal R']
for i in idx_DK_postcent:
    Region_DK[i]='P'
for i in idx_DK_Maudalmiddlefront:
    Region_DK[i]='M'

# partition into lobes - DK
nb_ROIs_DK=len(Region_DK)
Lobes_DK = [ "" for x in range ( len ( Region_DK ) ) ]
for kk_roi in range ( len ( Region_DK ) ):
    Lobes_DK[ kk_roi ] = regions_bst_to_lobes ( Region_DK[ kk_roi ] )

idx_F_DK = [ i for i , j in enumerate ( Lobes_DK ) if j == 'F' ]
idx_M_DK = [ i for i , j in enumerate ( Lobes_DK ) if j == 'M' ]
idx_P_DK = [ i for i , j in enumerate ( Lobes_DK ) if j == 'P' ]
idx_O_DK = [ i for i , j in enumerate ( Lobes_DK ) if j == 'O' ]
idx_T_DK = [ i for i , j in enumerate ( Lobes_DK ) if j == 'T' ]
idx_lobes_DK = [ idx_F_DK , idx_M_DK , idx_T_DK , idx_P_DK , idx_O_DK ]

# Destrieux
ROI_Destrieux_list = ['G_Ins_lg_and_S_Ment_ins L','G_Ins_lg_and_S_Ment_ins R','G_and_S_Mingul-Ant L','G_and_S_Mingul-Ant R','G_and_S_Mingul-Mid-Ant L','G_and_S_Mingul-Mid-Ant R','G_and_S_Mingul-Mid-Post L','G_and_S_Mingul-Mid-Post R','G_and_S_frontomargin L','G_and_S_frontomargin R','G_and_S_occipital_inf L','G_and_S_occipital_inf R','G_and_S_paracentral L','G_and_S_paracentral R','G_and_S_subcentral L','G_and_S_subcentral R','G_and_S_transv_frontopol L','G_and_S_transv_frontopol R','G_Mingul-Post-dorsal L','G_Mingul-Post-dorsal R','G_Mingul-Post-ventral L','G_Mingul-Post-ventral R','G_Muneus L','G_Muneus R','G_front_inf-Opercular L','G_front_inf-Opercular R','G_front_inf-Orbital L','G_front_inf-Orbital R','G_front_inf-Triangul L','G_front_inf-Triangul R','G_front_middle L','G_front_middle R','G_front_sup L','G_front_sup R','G_insular_short L','G_insular_short R','G_oc-temp_lat-fusifor L','G_oc-temp_lat-fusifor R','G_oc-temp_med-Lingual L','G_oc-temp_med-Lingual R','G_oc-temp_med-Parahip L','G_oc-temp_med-Parahip R','G_occipital_middle L','G_occipital_middle R','G_occipital_sup L','G_occipital_sup R','G_orbital L','G_orbital R','G_pariet_inf-Angular L','G_pariet_inf-Angular R','G_pariet_inf-Supramar L','G_pariet_inf-Supramar R','G_parietal_sup L','G_parietal_sup R','G_postcentral L','G_postcentral R','G_precentral L','G_precentral R','G_precuneus L','G_precuneus R','G_rectus L','G_rectus R','G_subcallosal L','G_subcallosal R','G_temp_sup-G_T_transv L','G_temp_sup-G_T_transv R','G_temp_sup-Lateral L','G_temp_sup-Lateral R','G_temp_sup-Plan_polar L','G_temp_sup-Plan_polar R','G_temp_sup-Plan_tempo L','G_temp_sup-Plan_tempo R','G_temporal_inf L','G_temporal_inf R','G_temporal_middle L','G_temporal_middle R','Lat_Fis-ant-Horizont L','Lat_Fis-ant-Horizont R','Lat_Fis-ant-Vertical L','Lat_Fis-ant-Vertical R','Lat_Fis-post L','Lat_Fis-post R','Pole_occipital L','Pole_occipital R','Pole_temporal L','Pole_temporal R','S_Malcarine L','S_Malcarine R','S_Mentral L','S_Mentral R','S_Mingul-Marginalis L','S_Mingul-Marginalis R','S_Mircular_insula_ant L','S_Mircular_insula_ant R','S_Mircular_insula_inf L','S_Mircular_insula_inf R','S_Mircular_insula_sup L','S_Mircular_insula_sup R','S_Mollat_transv_ant L','S_Mollat_transv_ant R','S_Mollat_transv_post L','S_Mollat_transv_post R','S_front_inf L','S_front_inf R','S_front_middle L','S_front_middle R','S_front_sup L','S_front_sup R','S_interm_prim-Jensen L','S_interm_prim-Jensen R','S_intrapariet_and_P_trans L','S_intrapariet_and_P_trans R','S_oc-temp_lat L','S_oc-temp_lat R','S_oc-temp_med_and_Lingual L','S_oc-temp_med_and_Lingual R','S_oc_middle_and_Lunatus L','S_oc_middle_and_Lunatus R','S_oc_sup_and_transversal L','S_oc_sup_and_transversal R','S_occipital_ant L','S_occipital_ant R','S_orbital-H_Shaped L','S_orbital-H_Shaped R','S_orbital_lateral L','S_orbital_lateral R','S_orbital_med-olfact L','S_orbital_med-olfact R','S_parieto_occipital L','S_parieto_occipital R','S_pericallosal L','S_pericallosal R','S_postcentral L','S_postcentral R','S_precentral-inf-part L','S_precentral-inf-part R','S_precentral-sup-part L','S_precentral-sup-part R','S_suborbital L','S_suborbital R','S_subparietal L','S_subparietal R','S_temporal_inf L','S_temporal_inf R','S_temporal_sup L','S_temporal_sup R','S_temporal_transverse L','S_temporal_transverse R']
Region_Destrieux = ['LT','RT','LL','RL','LL','RL','LL','RL','LPF','RPF','LO','RO','LC','RC','LC','RC','LPF','RPF','LL','RL','LL','RL','LO','RO','LF','RF','LF','RF','LF','RF','LF','RF','LF','RF','LT','RT','LT','RT','LO','RO','LT','RT','LO','RO','LO','RO','LPF','RPF','LP','RP','LP','RP','LP','RP','LC','RC','LC','RC','LP','RP','LPF','RPF','LL','RL','LT','RT','LT','RT','LT','RT','LT','RT','LT','RT','LT','RT','LF','RF','LF','RF','LT','RT','LO','RO','LT','RT','LO','RO','LC','RC','LC','RC','LT','RT','LT','RT','LT','RT','LT','RT','LO','RO','LF','RF','LF','RF','LF','RF','LP','RP','LP','RP','LT','RT','LT','RT','LO','RO','LO','RO','LO','RO','LPF','RPF','LPF','RPF','LPF','RPF','LP','RP','LL','RL','LC','RC','LC','RC','LC','RC','LPF','RPF','LP','RP','LT','RT','LT','RT','LT','RT']

nb_ROIs_Destrieux=len(Region_Destrieux)
Lobes_Destrieux=[1,1,1,1,1,1,3,3,1,1,4,4,5,5,1,1,1,1,3,3,3,3,4,4,1,1,1,1,1,1,5,5,1,1,2,2,2,2,2,2,2,2,4,4,4,4,1,1,3,3,3,3,3,3,3,3,5,5,3,3,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,4,4,2,2,4,4,5,5,3,3,2,2,2,2,2,2,2,2,2,2,1,1,5,5,1,1,3,3,3,3,4,4,4,4,4,4,4,4,4,4,1,1,1,1,1,1,3,3,1,1,3,3,5,5,5,5,1,1,3,3,2,2,2,2,2,2]

# cf subdivision made by PS
idx_F_Destrieux = [ i for i , j in enumerate ( Lobes_Destrieux ) if j == 1 ]
idx_M_Destrieux = [ i for i , j in enumerate ( Lobes_Destrieux ) if j == 5 ]
idx_P_Destrieux = [ i for i , j in enumerate ( Lobes_Destrieux ) if j == 3 ]
idx_O_Destrieux = [ i for i , j in enumerate ( Lobes_Destrieux ) if j == 4 ]
idx_T_Destrieux = [ i for i , j in enumerate ( Lobes_Destrieux ) if j == 2 ]
idx_lobes_Destrieux = [ idx_F_Destrieux , idx_M_Destrieux , idx_T_Destrieux , idx_P_Destrieux , idx_O_Destrieux ]

#%% Figure 1 - Methods
plt.close('all')
fig_corr = plt.figure()
plt.imshow(Test, cmap='RdBu', vmin=-max, vmax=max)
plt.colorbar()
plt.show()
fig_corr.savefig(fig_path+'Connect_Matrix_'+ f_name + '_p_thresh_0,001.png', dpi=600, facecolor='white')

plt.close('all')
fig_corr = plt.figure()
plt.imshow(sig_r_diff_val_thresh, cmap='RdBu_r', vmin=-max, vmax=max)
plt.colorbar()
plt.show()
fig_corr.savefig(fig_path+'Connect_Matrix_'+ f_name + '.png', dpi=600, facecolor='white')

# random shuffle columns matrix
rand_sig_r_diff_val_thresh = sig_r_diff_val_thresh[:, np.random.permutation(sig_r_diff_val_thresh.shape[1])]
rand_sig_r_diff_val_thresh = rand_sig_r_diff_val_thresh[np.random.permutation(sig_r_diff_val_thresh.shape[1]),:]
plt.close('all')
fig_corr = plt.figure()
plt.imshow(rand_sig_r_diff_val_thresh, cmap='RdBu_r', vmin=-max, vmax=max)
plt.colorbar()
plt.show()
fig_corr.savefig(fig_path+'Connect_Matrix_Rand_Shuffle'+ f_name + '.png', dpi=600, facecolor='white')



#%% Figura 2A - plots connectome - reliability
import mat73
p_thresh=0.05
#f_name = "Reliable_Analysis_Edge_Based_BH_p_thresh_edges_0.05_p_thresh_nodes_0.05_Group.mat"
f_name = "Reliable_Analysis_Edge_Based_BH_p_thresh_edges_0.05_p_thresh_nodes_0.05_Group_3_40_z_3_binning_2_min_size_aval_5_Sess4_MEG_DK"
data_dict = mat73.loadmat (db_path+f_name+'.mat')

sig_diff_val = data_dict[ 'signif_cmp' ]
sig_diff_val = sig_diff_val.astype ( int )

[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_diff_val, idx_F_DK , idx_M_DK , idx_P_DK , idx_T_DK , idx_O_DK )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_DK_list_F = [ROI_DK_list[i] for i in idx_F_DK]
ROI_DK_list_M = [ROI_DK_list[i] for i in idx_M_DK]
ROI_DK_list_T = [ROI_DK_list[i] for i in idx_T_DK]
ROI_DK_list_P = [ROI_DK_list[i] for i in idx_P_DK]
ROI_DK_list_O = [ROI_DK_list[i] for i in idx_O_DK]
ROI_DK_list_by_lobes=list(np.concatenate((ROI_DK_list_F,ROI_DK_list_M,ROI_DK_list_T,ROI_DK_list_P,ROI_DK_list_O)))

Region_DK_colors=Region_DK
for i in idx_F_DK:
    Region_DK_colors[i]='firebrick'
for i in idx_M_DK:
    Region_DK_colors[i]='darkorange'
for i in idx_T_DK:
    Region_DK_colors[i]='darkolivegreen'
for i in idx_P_DK:
    Region_DK_colors[i]='cadetblue'
for i in idx_O_DK:
    Region_DK_colors[i]='mediumpurple'

Region_DK_color_F = [Region_DK_colors[i] for i in idx_F_DK]
Region_DK_color_M = [Region_DK_colors[i] for i in idx_M_DK]
Region_DK_color_T = [Region_DK_colors[i] for i in idx_T_DK]
Region_DK_color_P = [Region_DK_colors[i] for i in idx_P_DK]
Region_DK_color_O = [Region_DK_colors[i] for i in idx_O_DK]
Region_DK_color_by_lobes=list(np.concatenate((Region_DK_color_F,Region_DK_color_M,Region_DK_color_T,Region_DK_color_P,Region_DK_color_O)))


sig_diff_val_plot=Test

fig_reliab = plt.figure ( figsize=(15 ,15) , facecolor='white' )
plot_connectivity_circle ( sig_diff_val_plot , node_names=ROI_DK_list_by_lobes , node_colors=Region_DK_color_by_lobes,
                                       colorbar=False,
                                       #colorbar_pos=(0.3 , -0.1) ,
                                       fontsize_names=21 , #fontsize_title=9 ,
                                       colormap='Greys', facecolor='white', textcolor='black',
                                       node_edgecolor='white',
                                       #title='p-values < 0.05, BH corrected, reliability' ,
                                       fig=fig_reliab )
fig_reliab.savefig ( fig_path + 'Connect_Reliable_Analysis_diff_Edge_Based_BH_p_thresh_edges_0.05_p_thresh_nodes_0.05_Group' + '.png' , dpi=600 ,
                   facecolor='white' )

#%% Figura 2B - matrix w/ all the edges - 0 preselezione
from mne_connectivity.viz import plot_connectivity_circle


f_name='Correlation_Analysis_BCI_Based_Nodi_Group_OptPerm_Paper_3_40_z_3_binning_2_min_size_aval_5_Sess4_MEG_DK_no_preselec_correl'
data_dict = scipy.io.loadmat(db_path+f_name+'.mat')
sig_r_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_r']
sig_p_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_p']
#sig_r_diff_val_thresh=sig_r_diff_val_thresh.astype ( int )

p_thresh=0.05#01
idx_signif=np.where(sig_p_diff_val_thresh<p_thresh)
sig_r_diff_val_thresh_test=np.zeros(np.shape(sig_r_diff_val_thresh))
sig_r_diff_val_thresh_test[idx_signif]=sig_r_diff_val_thresh[idx_signif]


[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_r_diff_val_thresh_test, idx_F_DK , idx_M_DK , idx_P_DK , idx_T_DK , idx_O_DK )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_DK_list_F = [ROI_DK_list[i] for i in idx_F_DK]
ROI_DK_list_M = [ROI_DK_list[i] for i in idx_M_DK]
ROI_DK_list_T = [ROI_DK_list[i] for i in idx_T_DK]
ROI_DK_list_P = [ROI_DK_list[i] for i in idx_P_DK]
ROI_DK_list_O = [ROI_DK_list[i] for i in idx_O_DK]
ROI_DK_list_by_lobes=list(np.concatenate((ROI_DK_list_F,ROI_DK_list_M,ROI_DK_list_T,ROI_DK_list_P,ROI_DK_list_O)))

Region_DK_colors=Region_DK
for i in idx_F_DK:
    Region_DK_colors[i]='firebrick'
for i in idx_M_DK:
    Region_DK_colors[i]='darkorange'
for i in idx_T_DK:
    Region_DK_colors[i]='darkolivegreen'
for i in idx_P_DK:
    Region_DK_colors[i]='cadetblue'
for i in idx_O_DK:
    Region_DK_colors[i]='mediumpurple'

Region_DK_color_F = [Region_DK_colors[i] for i in idx_F_DK]
Region_DK_color_M = [Region_DK_colors[i] for i in idx_M_DK]
Region_DK_color_T = [Region_DK_colors[i] for i in idx_T_DK]
Region_DK_color_P = [Region_DK_colors[i] for i in idx_P_DK]
Region_DK_color_O = [Region_DK_colors[i] for i in idx_O_DK]
Region_DK_color_by_lobes=list(np.concatenate((Region_DK_color_F,Region_DK_color_M,Region_DK_color_T,Region_DK_color_P,Region_DK_color_O)))

max=np.max(np.abs((Test)))
fig_corr = plt.figure(figsize=(15,15), facecolor='white')
plot_connectivity_circle ( con=Test , node_names=ROI_DK_list_by_lobes , node_colors=Region_DK_color_by_lobes,
                                    colorbar=True,
                                    colorbar_pos=(0.75 , -0.1) ,
                                    colormap='RdBu_r',
                                    facecolor='white', textcolor='black',
                                    node_edgecolor='white',
                                    #title='r-values, correlation (p<0.05)' ,
                                    fontsize_names=21 ,
                                    vmin=-max, vmax=max,
                                    fig=fig_corr )
fig_corr.savefig(fig_path+'Connect_'+ f_name + '_p_thresh_0,05_colorbar.png', dpi=600, facecolor='white')



#%% Supplementary - MEEG - DK/Destrieux - matrix w/ all the edges
from mne_connectivity.viz import plot_connectivity_circle

f_name='Correlation_Analysis_BCI_Based_Nodi_Group_OptPerm_Paper_3_40_z_3_binning_2_min_size_aval_5_Sess4_MEG_DK_no_preselec_correl'
data_dict = scipy.io.loadmat(db_path+f_name+'.mat')
sig_r_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_r']
sig_p_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_p']
[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_r_diff_val_thresh, idx_F_DK , idx_M_DK , idx_P_DK , idx_T_DK , idx_O_DK )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_DK_list_F = [ROI_DK_list[i] for i in idx_F_DK]
ROI_DK_list_M = [ROI_DK_list[i] for i in idx_M_DK]
ROI_DK_list_T = [ROI_DK_list[i] for i in idx_T_DK]
ROI_DK_list_P = [ROI_DK_list[i] for i in idx_P_DK]
ROI_DK_list_O = [ROI_DK_list[i] for i in idx_O_DK]
ROI_DK_list_by_lobes=list(np.concatenate((ROI_DK_list_F,ROI_DK_list_M,ROI_DK_list_T,ROI_DK_list_P,ROI_DK_list_O)))

Region_DK_colors=Region_DK
for i in idx_F_DK:
    Region_DK_colors[i]='firebrick'
for i in idx_M_DK:
    Region_DK_colors[i]='darkorange'
for i in idx_T_DK:
    Region_DK_colors[i]='darkolivegreen'
for i in idx_P_DK:
    Region_DK_colors[i]='cadetblue'
for i in idx_O_DK:
    Region_DK_colors[i]='mediumpurple'

Region_DK_color_F = [Region_DK_colors[i] for i in idx_F_DK]
Region_DK_color_M = [Region_DK_colors[i] for i in idx_M_DK]
Region_DK_color_T = [Region_DK_colors[i] for i in idx_T_DK]
Region_DK_color_P = [Region_DK_colors[i] for i in idx_P_DK]
Region_DK_color_O = [Region_DK_colors[i] for i in idx_O_DK]
Region_DK_color_by_lobes=list(np.concatenate((Region_DK_color_F,Region_DK_color_M,Region_DK_color_T,Region_DK_color_P,Region_DK_color_O)))

max=np.max(np.abs(Test[np.isnan(Test)==0]))
Test_r_thresh=Test
Test_r_thresh[np.where(abs(Test)<0.6)] = 0

fig_corr = plt.figure(figsize=(15,15), facecolor='white')
plot_connectivity_circle ( con=Test , node_names=ROI_DK_list_by_lobes , node_colors=Region_DK_color_by_lobes,
                                    colorbar=True,
                                    colorbar_pos=(0.75 , -0.1) ,
                                    colormap='RdBu_r',
                                    facecolor='white', textcolor='black',
                                    node_edgecolor='white',
                                    #title='r-values, correlation (p<0.05)' ,
                                    fontsize_names=21 ,
                                    vmin=-max, vmax=max,
                                    fig=fig_corr )
fig_corr.savefig(fig_path+'Connect_'+ f_name + '_r_thresh_0,6_colorbar.png', dpi=600, facecolor='white')


f_name='Correlation_Analysis_BCI_Based_Nodi_Group_OptPerm_Paper_3_40_z_3_binning_2_min_size_aval_5_Sess4_EEG_DK_no_preselec_correl'
data_dict = scipy.io.loadmat(db_path+f_name+'.mat')
sig_r_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_r']
sig_p_diff_val_thresh=data_dict['matr_edges_diff_corr_BCI_p']

[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_r_diff_val_thresh, idx_F_DK , idx_M_DK , idx_P_DK , idx_T_DK , idx_O_DK )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_DK_list_F = [ROI_DK_list[i] for i in idx_F_DK]
ROI_DK_list_M = [ROI_DK_list[i] for i in idx_M_DK]
ROI_DK_list_T = [ROI_DK_list[i] for i in idx_T_DK]
ROI_DK_list_P = [ROI_DK_list[i] for i in idx_P_DK]
ROI_DK_list_O = [ROI_DK_list[i] for i in idx_O_DK]
ROI_DK_list_by_lobes=list(np.concatenate((ROI_DK_list_F,ROI_DK_list_M,ROI_DK_list_T,ROI_DK_list_P,ROI_DK_list_O)))

Region_DK_colors=Region_DK
for i in idx_F_DK:
    Region_DK_colors[i]='firebrick'
for i in idx_M_DK:
    Region_DK_colors[i]='darkorange'
for i in idx_T_DK:
    Region_DK_colors[i]='darkolivegreen'
for i in idx_P_DK:
    Region_DK_colors[i]='cadetblue'
for i in idx_O_DK:
    Region_DK_colors[i]='mediumpurple'

Region_DK_color_F = [Region_DK_colors[i] for i in idx_F_DK]
Region_DK_color_M = [Region_DK_colors[i] for i in idx_M_DK]
Region_DK_color_T = [Region_DK_colors[i] for i in idx_T_DK]
Region_DK_color_P = [Region_DK_colors[i] for i in idx_P_DK]
Region_DK_color_O = [Region_DK_colors[i] for i in idx_O_DK]
Region_DK_color_by_lobes=list(np.concatenate((Region_DK_color_F,Region_DK_color_M,Region_DK_color_T,Region_DK_color_P,Region_DK_color_O)))

max=np.max(np.abs(Test[np.isnan(Test)==0]))
Test_r_thresh=Test
Test_r_thresh[np.where(abs(Test)<0.6)] = 0

fig_corr = plt.figure(figsize=(15,15), facecolor='white')
plot_connectivity_circle ( con=Test , node_names=ROI_DK_list_by_lobes , node_colors=Region_DK_color_by_lobes,
                                    colorbar=True,
                                    colorbar_pos=(0.75 , -0.1) ,
                                    colormap='RdBu_r',
                                    facecolor='white', textcolor='black',
                                    node_edgecolor='white',
                                    fontsize_names=21 ,
                                    vmin=-max, vmax=max,
                                    fig=fig_corr )
fig_corr.savefig(fig_path+'Connect_'+ f_name + '_r_thresh_0,6_colorbar.png', dpi=600, facecolor='white')

#%% Destrieux

f_name='Correlation_Analysis_BCI_Based_Nodi_Group_OptPerm_Paper_3_40_z_3_binning_2_min_size_aval_5_Sess4_MEG_Destrieux_no_preselec_correl'
data_dict = scipy.io.loadmat(db_path+f_name+'.mat')
sig_r_diff_val=data_dict['matr_edges_diff_corr_BCI_r']
sig_p_diff_val=data_dict['matr_edges_diff_corr_BCI_p']


# cf ways to denominate lobes differ from PS's parcellation and what has been done with DK
idx_F_Destrieux_plot=idx_F_Destrieux
idx_M_Destrieux_plot=idx_O_Destrieux
idx_T_Destrieux_plot=idx_M_Destrieux
idx_P_Destrieux_plot=idx_T_Destrieux
idx_O_Destrieux_plot=idx_P_Destrieux

[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_r_diff_val, idx_F_Destrieux, idx_M_Destrieux , idx_P_Destrieux , idx_T_Destrieux , idx_O_Destrieux )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_Destrieux_list_F = [ROI_Destrieux_list[i] for i in idx_F_Destrieux_plot]
ROI_Destrieux_list_M = [ROI_Destrieux_list[i] for i in idx_M_Destrieux_plot]
ROI_Destrieux_list_T = [ROI_Destrieux_list[i] for i in idx_T_Destrieux_plot]
ROI_Destrieux_list_P = [ROI_Destrieux_list[i] for i in idx_P_Destrieux_plot]
ROI_Destrieux_list_O = [ROI_Destrieux_list[i] for i in idx_O_Destrieux_plot]
ROI_Destrieux_list_by_lobes=list(np.concatenate((ROI_Destrieux_list_F,ROI_Destrieux_list_M,ROI_Destrieux_list_T,ROI_Destrieux_list_P,ROI_Destrieux_list_O)))

Region_Destrieux_colors=Region_Destrieux
for i in idx_F_Destrieux:
    Region_Destrieux_colors[i]='firebrick'
for i in idx_M_Destrieux:
    Region_Destrieux_colors[i]='mediumpurple'
for i in idx_T_Destrieux:
    Region_Destrieux_colors[i]='darkorange'
for i in idx_P_Destrieux:
    Region_Destrieux_colors[i]='darkolivegreen'
for i in idx_O_Destrieux:
    Region_Destrieux_colors[i]='cadetblue'




ROI_Destrieux_color_F = [Region_Destrieux_colors[i] for i in idx_F_Destrieux]
ROI_Destrieux_color_M = [Region_Destrieux_colors[i] for i in idx_M_Destrieux]
ROI_Destrieux_color_T = [Region_Destrieux_colors[i] for i in idx_T_Destrieux]
ROI_Destrieux_color_P = [Region_Destrieux_colors[i] for i in idx_P_Destrieux]
ROI_Destrieux_color_O = [Region_Destrieux_colors[i] for i in idx_O_Destrieux]
Region_Destrieux_color_by_lobes=list(np.concatenate((ROI_Destrieux_color_F,ROI_Destrieux_color_M,ROI_Destrieux_color_T,ROI_Destrieux_color_P,ROI_Destrieux_color_O)))


max=np.max(np.abs(Test[np.isnan(Test)==0]))
Test_r_thresh=Test
Test_r_thresh[np.where(abs(Test)<0.6)] = 0

fig_corr = plt.figure(figsize=(21,21), facecolor='white')
plot_connectivity_circle ( con=Test_r_thresh , node_names=ROI_Destrieux_list_by_lobes , node_colors=Region_Destrieux_color_by_lobes,
                                    colorbar=True,
                                    colorbar_pos=(0.55, -0.1),
                                    colormap='RdBu_r',
                                    facecolor='white', textcolor='black',
                                    node_edgecolor='white',
                                    #title='r-values, correlation (p<0.05)' ,
                                    fontsize_names=24 ,
                                    vmin=-max, vmax=max,
                                    fig=fig_corr )
fig_corr.savefig(fig_path+'Connect_'+ f_name + '_r_thresh_0,6_colorbar.png', dpi=600, facecolor='white')


f_name='Correlation_Analysis_BCI_Based_Nodi_Group_OptPerm_Paper_3_40_z_3_binning_2_min_size_aval_5_Sess4_EEG_Destrieux_no_preselec_correl'
data_dict = scipy.io.loadmat(db_path+f_name+'.mat')
sig_r_diff_val=data_dict['matr_edges_diff_corr_BCI_r']
sig_p_diff_val=data_dict['matr_edges_diff_corr_BCI_p']


# cf ways to denominate lobes differ from PS's parcellation and what has been done with DK
idx_F_Destrieux_plot=idx_F_Destrieux
idx_M_Destrieux_plot=idx_O_Destrieux
idx_T_Destrieux_plot=idx_M_Destrieux
idx_P_Destrieux_plot=idx_T_Destrieux
idx_O_Destrieux_plot=idx_P_Destrieux

[SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O , \
           SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O , \
           SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O , \
           SubMatrix_P_P , SubMatrix_P_O , \
           SubMatrix_O_O] = Lobes_Partition(sig_r_diff_val, idx_F_Destrieux, idx_M_Destrieux , idx_P_Destrieux , idx_T_Destrieux , idx_O_Destrieux )

temp1=np.hstack((SubMatrix_F_F , SubMatrix_F_M , SubMatrix_F_T , SubMatrix_F_P , SubMatrix_F_O))
temp2= np.hstack((np.transpose(SubMatrix_F_M), SubMatrix_M_M , SubMatrix_M_T , SubMatrix_M_P , SubMatrix_M_O))
temp1b=np.vstack((temp1,temp2))
temp3=np.hstack((np.transpose(SubMatrix_F_T) , np.transpose(SubMatrix_M_T), SubMatrix_T_T , SubMatrix_T_P , SubMatrix_T_O))
temp4=np.hstack((np.transpose(SubMatrix_F_P) , np.transpose(SubMatrix_M_P), np.transpose(SubMatrix_T_P) , SubMatrix_P_P , SubMatrix_P_O))
temp3b=np.vstack((temp3,temp4))
temp4b=np.vstack((temp1b,temp3b))
temp5=np.hstack((np.transpose(SubMatrix_F_O), np.transpose(SubMatrix_M_O), np.transpose(SubMatrix_T_O), np.transpose(SubMatrix_P_O), SubMatrix_O_O))
temp5b=np.vstack((temp4b,temp5))

Test=temp5b
ROI_Destrieux_list_F = [ROI_Destrieux_list[i] for i in idx_F_Destrieux_plot]
ROI_Destrieux_list_M = [ROI_Destrieux_list[i] for i in idx_M_Destrieux_plot]
ROI_Destrieux_list_T = [ROI_Destrieux_list[i] for i in idx_T_Destrieux_plot]
ROI_Destrieux_list_P = [ROI_Destrieux_list[i] for i in idx_P_Destrieux_plot]
ROI_Destrieux_list_O = [ROI_Destrieux_list[i] for i in idx_O_Destrieux_plot]
ROI_Destrieux_list_by_lobes=list(np.concatenate((ROI_Destrieux_list_F,ROI_Destrieux_list_M,ROI_Destrieux_list_T,ROI_Destrieux_list_P,ROI_Destrieux_list_O)))

Region_Destrieux_colors=Region_Destrieux
for i in idx_F_Destrieux:
    Region_Destrieux_colors[i]='firebrick'
for i in idx_M_Destrieux:
    Region_Destrieux_colors[i]='mediumpurple'
for i in idx_T_Destrieux:
    Region_Destrieux_colors[i]='darkorange'
for i in idx_P_Destrieux:
    Region_Destrieux_colors[i]='darkolivegreen'
for i in idx_O_Destrieux:
    Region_Destrieux_colors[i]='cadetblue'




ROI_Destrieux_color_F = [Region_Destrieux_colors[i] for i in idx_F_Destrieux]
ROI_Destrieux_color_M = [Region_Destrieux_colors[i] for i in idx_M_Destrieux]
ROI_Destrieux_color_T = [Region_Destrieux_colors[i] for i in idx_T_Destrieux]
ROI_Destrieux_color_P = [Region_Destrieux_colors[i] for i in idx_P_Destrieux]
ROI_Destrieux_color_O = [Region_Destrieux_colors[i] for i in idx_O_Destrieux]
Region_Destrieux_color_by_lobes=list(np.concatenate((ROI_Destrieux_color_F,ROI_Destrieux_color_M,ROI_Destrieux_color_T,ROI_Destrieux_color_P,ROI_Destrieux_color_O)))


max=np.max(np.abs(Test[np.isnan(Test)==0]))
Test_r_thresh=Test
Test_r_thresh[np.where(abs(Test)<0.6)] = 0

fig_corr = plt.figure(figsize=(21,21), facecolor='white')
plot_connectivity_circle ( con=Test_r_thresh , node_names=ROI_Destrieux_list_by_lobes , node_colors=Region_Destrieux_color_by_lobes,
                                    colorbar=True,
                                    colorbar_pos=(0.55, -0.1),
                                    colormap='RdBu_r',
                                    facecolor='white', textcolor='black',
                                    node_edgecolor='white',
                                    #title='r-values, correlation (p<0.05)' ,
                                    fontsize_names=24 ,
                                    vmin=-max, vmax=max,
                                    fig=fig_corr )
fig_corr.savefig(fig_path+'Connect_'+ f_name + '_r_thresh_0,6_colorbar.png', dpi=600, facecolor='white')
