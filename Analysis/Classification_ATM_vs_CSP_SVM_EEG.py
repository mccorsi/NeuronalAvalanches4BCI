#%%
"""
==============================================================
Attempt to classify EEG data in the source space - neuronal avalanches vs classical approaches - classification on longer trials w/ SVM
===============================================================

"""
# Authors:
#   Marie-Constance Corsi
#   Pierpaolo Sorrentino
# Date: 17/04/2023
# License: BSD (3-clause)

import matplotlib.pyplot as plt
import ptitprince as pt
import scipy.stats
import seaborn as sns

import pandas as pd
import mat73

from tqdm import tqdm
import pickle
import mne
from mne import create_info, EpochsArray
from mne.decoding import CSP as CSP_MNE

from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import ShuffleSplit, cross_val_score

import numpy as np
from scipy.stats import zscore
from moabb.paradigms import MotorImagery


#%%
# to be updated
root_path='/Users/marie-constance.corsi/Documents/GitHub/NeuronalAvalanches4BCI/';

db_path=root_path + 'Database/'
path_figures_root=root_path + 'Figures/'
path_csv_root=db_path
path_data_root=db_path

#%% functions

def transprob(aval,nregions): # (t,r)
    mat = np.zeros((nregions, nregions))
    norm = np.sum(aval, axis=0)
    for t in range(len(aval) - 1):
        ini = np.where(aval[t] == 1)
        mat[ini] += aval[t + 1]
    mat[norm != 0] = mat[norm != 0] / norm[norm != 0][:, None]
    return mat

def Transprob(ZBIN,nregions, val_duration): # (t,r)
    mat = np.zeros((nregions, nregions))
    A = np.sum(ZBIN, axis=1)
    a = np.arange(len(ZBIN))
    idx = np.where(A != 0)[0]
    aout = np.split(a[idx], np.where(np.diff(idx) != 1)[0] + 1)
    ifi = 0
    for iaut in range(len(aout)):
        if len(aout[iaut]) > val_duration:
            mat += transprob(ZBIN[aout[iaut]],nregions)
            ifi += 1
    mat = mat / ifi
    return mat,aout

def threshold_mat(data,thresh=3):
    current_data=data
    binarized_data=np.where(np.abs(current_data)>thresh,1,0)
    return (binarized_data)

def find_avalanches(data,thresh=3, val_duration=2):
    binarized_data=threshold_mat(data,thresh=thresh)
    N=binarized_data.shape[0]
    mat, aout = Transprob(binarized_data.T, N, val_duration)
    aout=np.array(aout,dtype=object)
    list_length=[len(i) for i in aout]
    unique_sizes=set(list_length)
    min_size,max_size=min(list_length),max(list_length)
    list_avalanches_bysize={i:[] for i in unique_sizes}
    for s in aout:
        n=len(s)
        list_avalanches_bysize[n].append(s)
    return(aout,min_size,max_size,list_avalanches_bysize, mat)
#%% load data from matlab
from moabb.datasets.base import BaseDataset
class EEG_DK_Dataset(BaseDataset):
    """
    Dataset from the NETBCI protocol, source space, EEG, DK, Sess4
    """

    def __init__(self):
        super().__init__(
            subjects=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],
            sessions_per_subject=1, # for the moment...
            events={"right_hand": 1, "rest": 2},
            code="NETBCI_EEG_DK_Dataset",
            interval=[1, 6], # starts when the target is displayed, ends before displaying the result -- longer trials
            paradigm="imagery",
            doi="", # not available yet...
        )

    def _get_single_subject_data(self, subject):
        """return data for a single subject"""
        file_path = self.data_path(subject)

        temp = mat73.loadmat(file_path[0]) # huge files...
        x = temp['Data_concat_moabb'][subject][0] # 2D: nb rois x (nb_trials x nb_samples)
        fs = temp["fs"]
        ch_names = temp['labels_DK'] # actually, roi names here...
        #subj_ID = temp['subject_IDs'][subject]
        events=temp["Events_moabb"][subject] #1: right_hand, 2:rest

        ch_types = ["eeg" for i in range(np.shape(ch_names)[0])]
        info = mne.create_info(ch_names, fs, ch_types)
        raw = mne.io.RawArray(data=np.array(x), info=info)

        mapping = {1: 'right_hand', 2: 'rest'}
        annot_from_events = mne.annotations_from_events(
            events=events, event_desc=mapping, sfreq=raw.info['sfreq'])
        raw.set_annotations(annot_from_events)

        return {"session_0": {"run_0": raw}}

    def data_path(
        self, path=None, force_update=False, update_path=None, verbose=None
    ):
        """Download the data from all the group"""
        path = '/Users/marieconstance.corsi/Documents/GitHub/Fenicotteri-equilibristi/Database/0_BCI/Classification/3_Dataset-netbci-EEG-sess4-DK/EEG_DK.mat'
        return [path]

#%%

from sklearn.pipeline import Pipeline
from sklearn.model_selection import ShuffleSplit, cross_val_score

dataset_EEG_DK = EEG_DK_Dataset()
subjects = dataset_EEG_DK.subject_list
events = ["right_hand", "rest"]

freqbands= {'paper': [3, 40]}

opt_zthresh = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
opt_val_duration = [2, 3, 4, 5, 6, 7, 8]

kk_components = 8

#%%
for subject in subjects:

    for f in freqbands:
        fmin = freqbands[f][0]
        fmax = freqbands[f][1]

        results = pd.DataFrame()

        paradigm_EEG = MotorImagery(
            events=events, n_classes=len(events), fmin=fmin, fmax=fmax
        )
        ep_EEG, labels, meta = paradigm_EEG.get_data(
            dataset=dataset_EEG_DK, subjects=[subject], return_epochs=True
        )
        epochs_data=ep_EEG.get_data()
        nb_ROIs=np.shape(epochs_data)[1]
        nb_trials=np.shape(epochs_data)[0]

        class_balance = np.mean(labels == labels[0])
        class_balance = max(class_balance, 1. - class_balance)

        # Preparation of the next steps - Computing CSP with MNE & applying it to time series to obtain new epochs
        csp_mne = CSP_MNE(n_components=kk_components, transform_into='csp_space').fit(epochs_data, labels)
        epochs_data_csp_m = csp_mne.transform(epochs_data)
        info = create_info([f'CSP{i}' for i in range(kk_components)], sfreq=ep_EEG.info['sfreq'], ch_types='eeg')
        ep_csp_m = EpochsArray(epochs_data_csp_m, info, ep_EEG.events, event_id=ep_EEG.event_id)

        # parameters for the default classifier & cross-validation
        svm = GridSearchCV(SVC(), {"kernel": ("linear", "rbf"), "C": [0.1, 1, 10]}, cv=3)
        csp = CSP_MNE(n_components=kk_components, reg=None, log=True, norm_trace=False) # average power
        nbSplit=50#75#50#10
        cv = ShuffleSplit(nbSplit,test_size=0.2, random_state=21)

    #%% CSP + SVM: classical
        clf_0 = Pipeline([('CSP', csp), ('SVM', svm)])
        score_CSP_SVM = cross_val_score(clf_0, epochs_data, labels, cv=cv, n_jobs=None)
        print("Classification accuracy CSP+SVM: %f / Chance level: %f" % (np.mean(score_CSP_SVM),
                                                                  class_balance))

        for kk_zthresh in opt_zthresh:
            for kk_val_duration in opt_val_duration:
                #%% ATM + SVM
                Nep_EEG, labels, meta = paradigm_EEG.get_data(
                    dataset=dataset_EEG_DK, subjects=[subject], return_epochs=False
                )
                temp = np.transpose(Nep_EEG, (1, 0, 2))
                temp_nc = np.reshape(temp, (np.shape(temp)[0], np.shape(temp)[1] * np.shape(temp)[2]))
                zscored_data = zscore(temp_nc, axis=1)
                # epoching here before computing the avalanches
                temp_zscored_data_ep = np.reshape(zscored_data, (np.shape(temp)[0], np.shape(temp)[1], np.shape(temp)[2]))
                zscored_data_ep = np.transpose(temp_zscored_data_ep, (1, 0, 2))

                ATM = np.empty((nb_trials, nb_ROIs, nb_ROIs))
                for kk_trial in range(nb_trials):
                    list_avalanches, min_size_avalanches, max_size_avalanches, list_avalanches_bysize, temp_ATM = find_avalanches(
                        zscored_data_ep[kk_trial, :, :], thresh=kk_zthresh, val_duration=kk_val_duration)
                    # ATM: nb_trials x nb_ROIs x nb_ROIs matrix
                    ATM[kk_trial, :, :] = temp_ATM

                clf_2 = Pipeline([('SVM', svm)])
                reshape_ATM = np.reshape(ATM, (np.shape(ATM)[0], np.shape(ATM)[1] * np.shape(ATM)[2]))
                score_ATM_SVM = cross_val_score(clf_2, reshape_ATM, labels, cv=cv, n_jobs=None)
                print("Classification accuracy ATM+SVM: %f / Chance level: %f" % (np.mean(score_ATM_SVM),
                                                                                  class_balance))


        #%% concatenate results in a single dataframe
                pd_CSP_SVM=pd.DataFrame(score_CSP_SVM, columns=["CSP+SVM"])
                pd_ATM_SVM = pd.DataFrame(score_ATM_SVM, columns=["ATM+SVM"])

                results = pd.concat([pd_CSP_SVM, pd_ATM_SVM,
                                     ],axis=1)

                c_scores = np.concatenate((score_CSP_SVM, score_ATM_SVM))
                ppl = np.concatenate((["CSP+SVM"]*len(score_CSP_SVM),
                                      ["ATM+SVM"]*len(score_ATM_SVM)))
                c_subject = [subject]*len(ppl)

                split = np.concatenate((np.arange(nbSplit),
                                        np.arange(nbSplit)))
                zthresh = [kk_zthresh]*len(ppl)
                val_duration = [kk_val_duration] * len(ppl)
                n_csp_comp = [kk_components]*len(ppl)
                freq = [str(fmin)+'-'+str(fmax)]*len(ppl)

                data2=np.transpose(np.vstack((c_scores,ppl, split, n_csp_comp, zthresh, val_duration, freq, c_subject)))
                results=pd.DataFrame(data2,columns=["score","pipeline","split",
                                                    "n_csp_comp", "zthresh", "val_duration", "freq",
                                                    "subject"])
                results.to_csv(
                    path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(subject) +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                    "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                         "_val_duration_" +str(kk_val_duration) + "_EEG.csv"
                )
                print(
                    "saved " +
                    path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(subject) +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                    "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                         "_val_duration_" +str(kk_val_duration) + "_EEG.csv"
                )

#%% study Diff
import statsmodels.stats.multitest
opt_zthresh = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
opt_val_duration = [2, 3, 4, 5, 6, 7, 8]
kk_components = 8
freqbands= {'paper': [3, 40]}

nbSplit=50
res_csv=pd.DataFrame()
res_diff_csv=pd.DataFrame()
res_csv_global = pd.DataFrame()
results_ttest=dict()
nb_split=list(range(0,nbSplit))

for f in freqbands:
    fmin = freqbands[f][0]
    fmax = freqbands[f][1]
    for kk_val_duration in opt_val_duration:
        for kk_zthresh in opt_zthresh:
            res_csv = pd.DataFrame()
            res_diff_csv = pd.DataFrame()
            x_Diff=[]
            vect_sign=[]
            pval_wilcoxon=[]
            pval_ttest=[]
            Diff = dict()
            Diff_median = []
            Diff_mean = []
            for subj in tqdm(subjects, desc="subject"):
                print(str(subj))
                temp = pd.read_csv(
                    path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(subj) +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                    "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                         "_val_duration_" +str(kk_val_duration) + "_EEG.csv")

                sc_pipeline=temp.loc[temp['pipeline'] == 'ATM+SVM', 'score']
                Baseline=temp.loc[temp['pipeline'] == 'CSP+SVM', 'score']

                Diff[subj]  = list(Baseline.values - sc_pipeline.values)

                if min(Diff[subj]) == 0 :
                    x_Diff = np.sign(max(Diff[subj]))
                elif max(Diff[subj]) == 0:
                    x_Diff = np.sign(min(Diff[subj]))
                else:
                    x_Diff = np.sign(min(Diff[subj]))*np.sign(max(Diff[subj]))
                vect_sign.append(x_Diff)
                Diff_median.append(np.median(Diff[subj]))
                Diff_mean.append(np.mean(Diff[subj]))
                results_wilcoxon = scipy.stats.wilcoxon(Baseline.values,sc_pipeline.values)
                results_ttest = scipy.stats.ttest_rel(Baseline.values,sc_pipeline.values)
                pval_wilcoxon.append(results_wilcoxon.pvalue)
                pval_ttest.append(results_ttest.pvalue)

            [reject_wilcoxon, pval_wilcoxon_corrected, alphacSidak, alphacBonf] = statsmodels.stats.multitest.multipletests(pvals=pval_wilcoxon,
                                                     alpha=0.05, method='fdr_bh')#method='bonferroni')
            [reject_ttest, pval_ttest_corrected, alphacSidak, alphacBonf] = statsmodels.stats.multitest.multipletests(pvals=pval_ttest,
                                                     alpha=0.05, method='fdr_bh')#method='bonferroni')


            filename = path_csv_root + "/SVM/IndivOpt_ComparisonDiff_SVM_Classification-allnode-2class-right_hand-rest-" + "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) + "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) + "_val_duration_" + str(kk_val_duration) + "_EEG"

            with open(filename + '.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                    pickle.dump([Diff, Diff_median, Diff_mean, vect_sign, pval_wilcoxon, pval_wilcoxon_corrected, pval_ttest, pval_ttest_corrected, nbSplit], f)

                # # Getting back the objects:
                # with open((filename + '.pkl') as f:  # Python 3: open(..., 'rb')
                #     Diff, vect_sign, results_wilcoxon,  results_ttest, nbSplit= pickle.load(f)


#%% study diff
opt_zthresh = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
opt_val_duration = [2, 3, 4, 5, 6, 7, 8]
kk_components = 8
freqbands= {'paper': [3, 40]}

for kk_zthresh in opt_zthresh:
    for f in freqbands:
        fmin = freqbands[f][0]
        fmax = freqbands[f][1]
        res_csv = pd.DataFrame()
        for kk_val_duration in opt_val_duration:
            filename = path_csv_root + "/SVM/IndivOpt_ComparisonDiff_SVM_Classification-allnode-2class-right_hand-rest-" + "n_csp_cmp-" + str(
                kk_components) + "zthresh-" + str(kk_zthresh) + "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(
                nbSplit) + "_val_duration_" + str(kk_val_duration) + "_EEG"
            # # Getting back the objects:
            with open(filename + '.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
                 Diff, Diff_median, Diff_mean, vect_sign, pval_wilcoxon, pval_wilcoxon_corrected, pval_ttest, pval_ttest_corrected, nbSplit= pickle.load(f)

            vect_sign=list(vect_sign)
            pval_wilcoxon=list(pval_wilcoxon)
            pval_wilcoxon_corrected = list(pval_wilcoxon_corrected)
            pval_ttest=list(pval_ttest)
            pval_ttest_corrected = list(pval_ttest_corrected)
            Diff_median= list(Diff_median)
            Diff_mean = list(Diff_mean)
            results_Diff=pd.DataFrame.from_dict(Diff)
            results_stats=pd.DataFrame(list(zip(vect_sign, Diff_median, Diff_mean, pval_wilcoxon, pval_wilcoxon_corrected, pval_ttest, pval_ttest_corrected)), columns=["vect_sign", "Diff_median", "Diff_mean", "pval_wilcoxon", "pval_wilcoxon_corrected", "pval_ttest", "pval_ttest_corrected"])

            results_Diff.to_csv(
                path_csv_root + "/SVM/Exhaustive_DiffComparison_SVM_Classification-allnode-2class-right_hand-rest-" +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                     "_val_duration_" +str(kk_val_duration) + "_EEG.csv"
            )
            results_stats.to_csv(
                path_csv_root + "/SVM/Exhaustive_StatsComparison_SVM_Classification-allnode-2class-right_hand-rest-" +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                     "_val_duration_" +str(kk_val_duration) + "_EEG.csv"
            )

#%% Load results for statistical analysis

opt_zthresh = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
opt_val_duration = [2, 3, 4, 5, 6, 7, 8]
kk_components = 8
freqbands= {'paper': [3, 40]}
fmin=3
fmax=40
res_opt = pd.DataFrame()
res_stats_opt = pd.DataFrame()
OptConfig_sc_ATM = np.empty((20, 4))
Diff_opt_cfg = np.empty((nbSplit, 20))
OptConfig_sc_ATM_std = np.empty((20, 4))
for subj in tqdm(subjects, desc="subject"):
    res = pd.DataFrame()
    res_stats = pd.DataFrame()
    res_diff = pd.DataFrame()
    res_score=[]
    kk=0
    cfg=np.empty((len(opt_zthresh)*len(opt_val_duration),2))

    for i_zth, kk_zthresh in enumerate(opt_zthresh):
        for i_val, kk_val_duration in enumerate(opt_val_duration):
            print(str(subj))
            temp_res = pd.read_csv(
                path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(
                    subj) + "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) +
                "_val_duration_" + str(kk_val_duration) + "_EEG.csv")

            res=pd.concat((res,temp_res))

            temp_stat = pd.read_csv(
                path_csv_root + "/SVM/Exhaustive_StatsComparison_SVM_Classification-allnode-2class-right_hand-rest-" + "n_csp_cmp-" + str(
                    kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) +
                "_val_duration_" + str(kk_val_duration) + "_EEG.csv") # 1 file for all subjects here

            res_stats=pd.concat((res_stats,temp_stat.iloc[[subj]]))

            temp_diff = pd.read_csv(
                path_csv_root + "/SVM/Exhaustive_DiffComparison_SVM_Classification-allnode-2class-right_hand-rest-" +  "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) +'-'+ str(fmax) +'-nbSplit' + str(nbSplit) +
                "_val_duration_" +str(kk_val_duration) + "_EEG.csv") # 1 file for all subjects here

            temp2=temp_diff[str(subj)]
            res_diff= pd.concat((res_diff,temp2))

            # scores
            temp_scores=pd.read_csv(
                path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(
                    subj) + "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(kk_zthresh) +
                "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) +
                "_val_duration_" + str(kk_val_duration) + "_EEG.csv"
            )
            sc_ATM=temp_scores[temp_scores["pipeline"]=="ATM+SVM"]
            res_score=np.concatenate((res_score,[sc_ATM["score"].mean()]))#mean()]))

            cfg[kk, :] = np.array([kk_zthresh, kk_val_duration])
            kk = kk + 1

    cfg_opt = np.where(res_score == max(res_score))

    optimal_zthresh = cfg[cfg_opt[0][0], 0]
    if optimal_zthresh==2.0:
        optimal_zthresh=2
    elif optimal_zthresh==3.0:
        optimal_zthresh=3
    optimal_val_duration = cfg[cfg_opt[0][0], 1]
    if optimal_val_duration==2.0:
        optimal_val_duration=2
    elif optimal_val_duration==3.0:
        optimal_val_duration=3
    elif optimal_val_duration == 4.0:
        optimal_val_duration = 4
    elif optimal_val_duration == 5.0:
        optimal_val_duration = 5
    elif optimal_val_duration == 6.0:
        optimal_val_duration = 6
    elif optimal_val_duration == 7.0:
        optimal_val_duration = 7
    elif optimal_val_duration == 8.0:
        optimal_val_duration = 8
    print("Optimal config - zthresh:" + str(optimal_zthresh) + "& val_duration:" +
          str(optimal_val_duration))

    temp_res_opt = pd.read_csv(
        path_csv_root + "/SVM/IndivOpt_Comparison_SVM_Classification-allnode-2class-right_hand-rest-subject-" + str(
            subj) + "n_csp_cmp-" + str(kk_components) + "zthresh-" + str(optimal_zthresh) +
        "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) +
        "_val_duration_" + str(int(optimal_val_duration)) + "_EEG.csv"
    )
    res_opt = pd.concat((res_opt, temp_res_opt))

    sc_pipeline = temp_res_opt.loc[temp_res_opt['pipeline'] == 'ATM+SVM', 'score']
    Baseline = temp_res_opt.loc[temp_res_opt['pipeline'] == 'CSP+SVM', 'score']

    temp_Diff_opt_cfg = Baseline.values - sc_pipeline.values
    Diff_opt_cfg[:,subj]=temp_Diff_opt_cfg

    temp_stat_opt = pd.read_csv(
        path_csv_root + "/SVM/Exhaustive_StatsComparison_SVM_Classification-allnode-2class-right_hand-rest-" + "n_csp_cmp-" + str(
            kk_components) + "zthresh-" + str(optimal_zthresh) +
        "-freq-" + str(fmin) + '-' + str(fmax) + '-nbSplit' + str(nbSplit) +
        "_val_duration_" + str(optimal_val_duration) + "_EEG.csv")  # 1 file for all subjects here

    res_stats_opt = pd.concat((res_stats_opt, temp_stat_opt.iloc[[subj]]))


    OptConfig_sc_ATM[subj,:] = [optimal_zthresh, optimal_val_duration, Baseline.values.mean(),sc_pipeline.values.mean()]
    OptConfig_sc_ATM_std[subj, :] = [optimal_zthresh, optimal_val_duration, Baseline.values.std(),
                                 sc_pipeline.values.std()]

df_optcfg=pd.DataFrame(OptConfig_sc_ATM, columns=["zthresh", "val_duration", "CSP+SVM", "ATM+SVM"])
df_optcfg_std=pd.DataFrame(OptConfig_sc_ATM_std, columns=["zthresh", "val_duration", "CSP+SVM", "ATM+SVM"])

#%% plot
path_figures_root = "/Users/marieconstance.corsi/Documents/GitHub/Fenicotteri-equilibristi/Figures/Classification/"

plt.style.use("dark_background")
g = sns.catplot(x="score",
                y='pipeline',
                hue="pipeline",
                kind='swarm',
                height=6, aspect=2,
                data=res_opt)
plt.savefig(path_figures_root + "IndivOpt_SVM_Classification_"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_EEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.pdf", dpi=300)

plt.style.use("dark_background")
g = sns.catplot(x="score",
                y='pipeline',
                col="subject",
                col_wrap=4,
                hue="pipeline",
                kind='box', # swarm
                height=6, aspect=2,
                data=res_opt)
plt.savefig(path_figures_root + "IndivOpt_SVM_Classification_IndivPlot"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_EEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.pdf", dpi=300)

#%% mean scores ATM - opt cfg
res_opt_ATM=res_opt[res_opt["pipeline"]=="ATM+SVM"]

mean_sc_ATM=[]
for subj in subjects:
    temp_subj=res_opt_ATM[res_opt_ATM["subject"]==subj]
    mean_sc_ATM.append(temp_subj["score"].mean())

mean_sc_ATM=np.array(mean_sc_ATM)

#%% plot best config vs mean(score_ATM)
plt.close('all')
sns.scatterplot(data=df_optcfg, x="zthresh", y="val_duration", hue="ATM+SVM", size="ATM+SVM", palette="magma", legend="brief")
plt.savefig(path_figures_root + "BestConfig_Study_IndivOpt_SVM_Classification_IndivPlot"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_EEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.pdf", dpi=300)


test_sc_pipeline = df_optcfg['ATM+SVM']
test_Baseline = df_optcfg['CSP+SVM']
test_Diff= test_Baseline - test_sc_pipeline

pd_Diff=pd.DataFrame()
pd_ATM=pd.DataFrame()
pd_CSP=pd.DataFrame()

pd_Diff=pd.DataFrame(test_Diff)
pd_Diff["pipeline"]=["Diff"]*len(pd_Diff)
pd_Diff.columns=["Avgscore", "pipeline"]

pd_ATM=pd.DataFrame(test_sc_pipeline)
pd_ATM["pipeline"]=["ATM+SVM"]*len(pd_ATM)
pd_ATM.columns=["Avgscore", "pipeline"]

pd_CSP=pd.DataFrame(test_Baseline)
pd_CSP["pipeline"]=["CSP+SVM"]*len(pd_CSP)
pd_CSP.columns=["Avgscore", "pipeline"]

pd_plot=pd.concat((pd_CSP,pd_ATM,pd_Diff))

# swarm diff
pd_plot_Diff = pd_plot[pd_plot["pipeline"]=="Diff"]
plt.style.use("dark_background")
g = sns.catplot(x="Avgscore",
                y='pipeline',
                hue="pipeline",
                kind='swarm', # swarm
                height=6, aspect=2,
                data=pd_plot_Diff)
plt.savefig(path_figures_root + "SwarmMeanDiff_IndivOpt_SVM_Classification_IndivPlot"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_EEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.pdf", dpi=300)

# swarm avg scores x pipelines
list_ppl=["ATM+SVM","CSP+SVM"]
pd_plot_avgsc= pd_plot[pd_plot["pipeline"].isin(list_ppl)]
plt.style.use("dark_background")
g = sns.catplot(x="Avgscore",
                y='pipeline',
                hue="pipeline",
                kind='swarm', # swarm
                height=6, aspect=2,
                data=pd_plot_avgsc)
plt.savefig(path_figures_root + "SwarmMeanScores_IndivOpt_SVM_Classification_IndivPlot"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_EEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.pdf", dpi=300)
#%% plot raincloud CSP+SVM vs ATM+ SVM
filename = "OptConfig_AvgSplits_Raincloud_"
plt.style.use("dark_background")
ort = "h"
pal = "Set2"
sigma = 0.2
dx = "pipeline"
dy = "Avgscore"
dhue = "pipeline"
f, ax = plt.subplots(figsize=(12, 9))
ax = pt.RainCloud(
    x=dx,
    y=dy,
    hue=dhue,
    data=pd_plot_avgsc,
    palette=pal,
    bw=sigma,
    #width_viol=0.7,
    ax=ax,
    orient=ort,
    alpha=0.65,
    dodge=True,
    pointplot=True,
    move=0.2,
)
ax.get_legend().remove()
ax.xaxis.label.set_size(9)
ax.yaxis.label.set_size(9)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.xlim((0, 1))  # to have always the same scale
# plt.yticks(rotation=45)
plt.savefig(path_figures_root + filename + "_FE_EEG.pdf", dpi=300)

#%% paired plot - CSP+SVM vs ATM+SVM - best configurations
import matplotlib.pyplot as plt
plt.style.use('default')
plt.style.use('seaborn-colorblind')

test=pd.DataFrame(list(np.concatenate((subjects,subjects))), columns=["subject"])
pd_plot_avgsc_concat= pd.concat([pd_plot_avgsc.reset_index(drop=True), test.reset_index(drop=True)],axis=1)

signif=dict()
for kk in list(range(len(res_stats_opt))):
    if res_stats_opt["pval_ttest_corrected"][kk]<0.05:
        if np.sign(res_stats_opt["Diff_median"][kk])==-1:
            signif[kk] ="ATM+SVM better"
        elif np.sign(res_stats_opt["Diff_median"][kk])==1:
            signif[kk] = "CSP+SVM better"
    else:
        signif[kk] = "equals"

temp_stat_plot=pd.DataFrame.from_dict(data=signif, orient='index', columns=["signif"])
temp_stat_plot2=pd.concat((temp_stat_plot,temp_stat_plot))
pd_plot_avgsc_concat2=pd.concat([pd_plot_avgsc_concat.reset_index(drop=True),temp_stat_plot2.reset_index(drop=True)], axis=1)

fig, axes = plt.subplots(figsize=[15,15])
paired = pd_plot_avgsc_concat2.pivot_table(
    values="Avgscore", columns="pipeline", index=["subject", "signif"]
)
paired = paired.reset_index()

fig, axes = plt.subplots(figsize=[15,15])
test_magma=sns.color_palette("magma")
colors = [test_magma[4], test_magma[1], [0.5,0.5,0.5]]# Set your custom color palette
palette_signif=sns.color_palette(colors)
g = sns.JointGrid(data=paired, x=paired["ATM+SVM"], y=paired["CSP+SVM"], hue=paired["signif"], palette=palette_signif)
g.plot_joint(sns.scatterplot, s=200, alpha=.5)
g.ax_joint.plot([0.5, 1], [0.5, 1], ls="--", c="k", linewidth=3)
g.ax_marg_x.hist(paired["ATM+SVM"], color=palette_signif[0], alpha=.6)
g.ax_marg_y.hist(paired["CSP+SVM"], color=palette_signif[1], alpha=.6, orientation="horizontal")
g.fig.axes[0].spines["bottom"].set_linewidth(6)
g.fig.axes[0].spines["left"].set_linewidth(6)
g.set_axis_labels('ATM+SVM', 'CSP+SVM', fontsize=24, weight='bold')
g.ax_joint.legend_.remove()
g.ax_joint.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
g.ax_joint.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
g.ax_joint.tick_params(labelsize=21)
for xlabel, ylabel in zip(g.ax_joint.get_xticklabels(), g.ax_joint.get_yticklabels()):
    ylabel.set_weight('bold')
    xlabel.set_weight('bold')
    xlabel.set_rotation(45)
filename = "OptConfig_AvgSplits_Comparison_"
g.savefig(path_figures_root + filename + "__EEG.png", dpi=300)

plt.close('all')
plt.style.use('default')
plt.style.use('seaborn-colorblind')
fig, axes = plt.subplots(figsize=[15,15])
g=sns.catplot(data=res_opt, x="pipeline", y="score", kind="swarm", palette='magma',
              aspect=1.2, height=3.6, s=1.8)
g.ax.set(xlabel=None, ylabel=None)
g.ax.tick_params(axis='both', which='major', labelsize=12)
for xlabel in g.ax.get_xticklabels():
    xlabel.set_weight('bold')
for ylabel in g.ax.get_yticklabels():
    ylabel.set_weight('bold')
g.ax.spines["bottom"].set_linewidth(6)
g.ax.spines["top"].set_visible(False)
g.ax.spines["left"].set_linewidth(6)
g.ax.spines["right"].set_visible(False)
plt.savefig(path_figures_root + "IndivOpt_SVM_Classification_"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_MEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.png", dpi=300)

#%% plot variance
test_sc_pipeline = df_optcfg_std['ATM+SVM']
test_Baseline = df_optcfg_std['CSP+SVM']

pd_ATM_std=pd.DataFrame()
pd_CSP_std=pd.DataFrame()

pd_ATM_std=pd.DataFrame(test_sc_pipeline)
pd_ATM_std["pipeline"]=["ATM+SVM"]*len(pd_ATM_std)
pd_ATM_std.columns=["STD_score", "pipeline"]

pd_CSP_std=pd.DataFrame(test_Baseline)
pd_CSP_std["pipeline"]=["CSP+SVM"]*len(pd_CSP_std)
pd_CSP_std.columns=["STD_score", "pipeline"]

pd_plot_std=pd.concat((pd_CSP_std,pd_ATM_std))

import matplotlib.pyplot as plt
plt.style.use('default')
plt.style.use('seaborn-colorblind')
fig, axes = plt.subplots()
sns.stripplot(
    data=pd_plot_std,
    y="STD_score",
    x="pipeline",
    ax=axes,
    jitter=True,
    alpha=0.5,
    zorder=1,
    palette="magma",
)
sns.pointplot(data=pd_plot_std, y="STD_score", x="pipeline", ax=axes, zorder=1, palette="magma")

axes.set_ylabel("Intra-subject variance")
axes.set_ylim(0, 0.1)
plt.savefig(path_figures_root + "IndivOpt_SVM_Classification_STD_"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_MEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.png", dpi=300)

#%%
fig, axes = plt.subplots(figsize=[15,15])
test_magma=sns.color_palette("magma")
colors = [test_magma[4], test_magma[1], [0.5,0.5,0.5]]# Set your custom color palette
palette_signif=sns.color_palette(colors)
x_std=pd_plot_std.STD_score[pd_plot_std["pipeline"]=="ATM+SVM"]
y_std=pd_plot_std.STD_score[pd_plot_std["pipeline"]=="CSP+SVM"]
g = sns.JointGrid(data=pd_plot_std, x=x_std, y=y_std, hue=paired["signif"], palette=palette_signif)
g.plot_joint(sns.scatterplot, s=200, alpha=.5)
g.ax_joint.plot([0, 0.1], [0, 0.1], ls="--", c="k", linewidth=3)
g.ax_marg_x.hist(x_std, color=palette_signif[0], alpha=.6)
g.ax_marg_y.hist(y_std, color=palette_signif[1], alpha=.6, orientation="horizontal")
g.fig.axes[0].spines["bottom"].set_linewidth(6)
g.fig.axes[0].spines["left"].set_linewidth(6)
g.set_axis_labels('ATM+SVM', 'CSP+SVM', fontsize=24, weight='bold')
g.ax_joint.legend_.remove()
g.ax_joint.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10], )
g.ax_joint.set_yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10],weight ="bold")
g.ax_joint.tick_params(labelsize=21)
for xlabel, ylabel in zip(g.ax_joint.get_xticklabels(), g.ax_joint.get_yticklabels()):
    ylabel.set_weight('bold')
    xlabel.set_weight('bold')
    xlabel.set_rotation(45)
g.savefig(path_figures_root + "IndivOpt_SVM_Classification_STD_joint_"+str(fmin)+'-'+str(fmax)+"NETBCI_GlobalResults_MEG_DK_Sess4-2class_rest_rh_nbSplits"+str(nbSplit)+"_EEG.png", dpi=300)


#%% compute inter-subject variance over all the splits
res_opt_CSP_SVM=res_opt[res_opt["pipeline"]=="CSP+SVM"]
res_opt_CSP_SVM.score.std()
res_opt_ATM.score.std()

res_opt_CSP_SVM.score.min()
res_opt_CSP_SVM.score.max()
res_opt_ATM.score.min()
res_opt_ATM.score.max()

#%% compute intra-subject variability
pd_plot_std[pd_plot_std["pipeline"]=="CSP+SVM"].STD_score.median()
pd_plot_std[pd_plot_std["pipeline"]=="ATM+SVM"].STD_score.median()