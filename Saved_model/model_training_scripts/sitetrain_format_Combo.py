import os
import gzip
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from numpy.lib.stride_tricks import sliding_window_view
from sitetrain_format_Combo import *
#
def get_args():
    parser = argparse.ArgumentParser()
    # ZMW prediction files
    parser.add_argument("-hifizmw", "--HiFiZMWpre", required=True, help="Molecule level 5mCpG prediction using HiFi reads. File format: txt.gz")
    parser.add_argument("-subreadszmw", "--SubreadsZMWpre", required=True, help="Molecule level 5mCpG prediction using subreads. File format: txt.gz")
    parser.add_argument("-minZMWdep", "--minZMWdepth", required=False, default='4', type=int, help="min ZMW_molecule/HiFi_reads/CCS_reads depth")
    parser.add_argument("-r","--windowRadius", required=False, default='5', type=int, help="how many CpG sites in a sliding window? (2 * r + 1)")
    # Ground-truth filles
    parser.add_argument("-trueMe", "--TrueMeFile", required=True, help="true DNAm dataset, in bed format, and indexed by Cytosine position in Forward Strand ( 1-based coordinate system )")
    parser.add_argument("-dep", "--covdep",  required=True, default='10', help="true DNAm dataset reads depth of coverage, e.g., wgbs reads >= 5X: -dep 5")
    # output files
    parser.add_argument("-o", "--outputfile", required=True, help="outputfile prefix")
    parser.add_argument("-maxCpGnum", "--maxCpGnumPerFile", type=int, default='1000',help="Max CpG sites for each output file, default is 1000 CpG sites per file")
    # include or exclude sites
    parser.add_argument("-include", "--includecglist", default="NA", type=str, help="Include features from a list of the CpGs, in the format of chr3_323223")
    parser.add_argument("-exclude", "--excludecglist", default="NA", type=str, help="Exclude features from a list of the CpGs, in the format of chr3_323223")
    return parser.parse_args()
#
if __name__ == "__main__":
    # --------
    args = get_args()
    #
    ZMWprefile_HiFi = args.HiFiZMWpre
    ZMWprefile_ZMW = args.SubreadsZMWpre
    minZMWdepth = int(args.minZMWdepth)
    R = int(args.windowRadius)
    trueMefile=args.TrueMeFile
    cutoffdep = int(args.covdep)
    outputfile = args.outputfile
    maxCpGlimit = int(args.maxCpGnumPerFile)
    #
    include_list_file = args.includecglist
    exclude_list_file = args.excludecglist
    # ------------------------------------------------------------------------------------------
    if include_list_file != "NA":
        target_cg_list = get_cglist(include_list_file)
    else:
        target_cg_list = "NA"
    # ------------------------------------------------------------------------------------------
    if exclude_list_file != "NA":
        target_ex_cg_list = get_cglist(exclude_list_file)
    else:
        target_ex_cg_list = "NA"
    # ------------------------------------------------------------------------------------------
    # read in true wgbs labelling file, key: "Chr_Pos":wgbs_value
    trueMeDict = {}
    trueMeDict = getMethyInfo(trueMefile, cutoffdep)
    # ------------------------------------------------------------------------------------------
    # read in files and fomratting
    # key: "Chr_Pos", value:[zmw1prob,zmw2prob, ... ]
    #
    ZMWprobDict_norm = {}
    ZMWprobDict_norm = normZMWprediction_HCO(ZMWprefile_HiFi, ZMWprefile_ZMW, minZMWdepth)

    # ------------------------------------------------------------------------------------------
    # intersect the two dicts
    update_ZMWprobDict_norm, WGBS_label = intersect_function(ZMWprobDict_norm, trueMeDict)


    # ------------------------------------------------------------------------------------------
    # keyof_X_train is chr_pos
    X_train, keyof_X_train = CreateDataMatrix_distance(update_ZMWprobDict_norm, R)

    #
    # ------------------------------------------------------------------------------------------
    # 2024-10-15: include and exclude list
    if target_cg_list != "NA" and target_ex_cg_list == "NA":
        shared_target_cg_list = [x for x in target_cg_list if x in keyof_X_train]
        subset_cg_index = [keyof_X_train.index(cg) for cg in shared_target_cg_list]
        y_train = [WGBS_label[x] for x in shared_target_cg_list]
        y_train = np.array(y_train)/100
        X_train_subset = X_train[subset_cg_index]
        del(X_train)
        n_samples = len(subset_cg_index)

    if target_ex_cg_list != "NA" and target_cg_list == "NA":
        shared_target_cg_list = [x for x in keyof_X_train if x not in set(target_ex_cg_list)]
        trimed_keyofX_dict = {cg: i for i, cg in enumerate(keyof_X_train)}
        subset_cg_index = [trimed_keyofX_dict[cg] for cg in shared_target_cg_list if cg in trimed_keyofX_dict]
        y_train = [WGBS_label[x] for x in shared_target_cg_list]
        y_train = np.array(y_train)/100
        X_train_subset = X_train[subset_cg_index]
        del(X_train)
        n_samples = len(subset_cg_index)

    if target_cg_list == "NA" and target_ex_cg_list == "NA":
        y_train = [WGBS_label[x] for x in keyof_X_train]
        y_train = np.array(y_train)/100
        X_train_subset = X_train
        n_samples = len(keyof_X_train)

    # ------------------------------------------------------------------------------------------
    #
    #
    # n_samples = len(subset_cg_index)

    batch_size = maxCpGlimit
    n_batches = (n_samples + batch_size - 1) // batch_size

    # ----------------------------------------------------------------------------------------
    processed_0_20_cpg_num = 0
    processed_20_40_cpg_num = 0
    processed_40_60_cpg_num = 0
    processed_60_80_cpg_num = 0
    processed_80_100_cpg_num = 0
    #
    FileIndex_0_20 = 0
    FileIndex_20_40 = 0
    FileIndex_40_60 = 0
    FileIndex_60_80 = 0
    FileIndex_80_100 = 0
    #
    #
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, n_samples)
        # X_batch = X_train[start_idx:end_idx]
        X_batch = X_train_subset[start_idx:end_idx]
        Y_batch = y_train[start_idx:end_idx]
        #
        # 2023-08-29: split into 0-20, 20-40, 40-60, 60-80, 80-100 groups and output to files
        # ----------------------------------------------
        # wgbs: 0-20
        #
        Y_batch_0_20_indices = np.where( (Y_batch >= 0) & (Y_batch <= 0.2) )
        if len(Y_batch_0_20_indices) > 0:
            if ('X_batch_0_20' not in globals()) and ('Y_batch_0_20' not in globals()):
                X_batch_0_20 = X_batch[Y_batch_0_20_indices]
                Y_batch_0_20 = Y_batch[Y_batch_0_20_indices]
                processed_0_20_cpg_num = processed_0_20_cpg_num + len(Y_batch_0_20)
            else:
                x = X_batch[Y_batch_0_20_indices]
                y = Y_batch[Y_batch_0_20_indices]
                X_batch_0_20 = np.concatenate((X_batch_0_20, x), axis = 0)
                Y_batch_0_20 = np.concatenate((Y_batch_0_20, y), axis = 0)
                processed_0_20_cpg_num = processed_0_20_cpg_num + len(y)
            if processed_0_20_cpg_num >= maxCpGlimit:
                FileIndex_0_20 = FileIndex_0_20 + 1
                npz_file_name_0_20 = outputfile + "." + str(FileIndex_0_20) + ".wgbs_0_20.chunk.npz"
                np.savez_compressed(npz_file_name_0_20, X = X_batch_0_20, Y = Y_batch_0_20)
                processed_0_20_cpg_num = 0
                del(X_batch_0_20)
                del(Y_batch_0_20)
        # ----------------------------------------------
        # wgbs: 20-40
        #
        Y_batch_20_40_indices = np.where( (Y_batch > 0.2) & (Y_batch <= 0.4) )
        if len(Y_batch_20_40_indices) > 0:
            if ('X_batch_20_40' not in globals()) and ('Y_batch_20_40' not in globals()):
                X_batch_20_40 = X_batch[Y_batch_20_40_indices]
                Y_batch_20_40 = Y_batch[Y_batch_20_40_indices]
                processed_20_40_cpg_num = processed_20_40_cpg_num + len(Y_batch_20_40)
            else:
                x = X_batch[Y_batch_20_40_indices]
                y = Y_batch[Y_batch_20_40_indices]
                X_batch_20_40 = np.concatenate((X_batch_20_40, x), axis = 0)
                Y_batch_20_40 = np.concatenate((Y_batch_20_40, y), axis = 0)
                processed_20_40_cpg_num = processed_20_40_cpg_num + len(y)
            if processed_20_40_cpg_num >= maxCpGlimit:
                FileIndex_20_40 = FileIndex_20_40 + 1
                npz_file_name_20_40 = outputfile + "." + str(FileIndex_20_40) + ".wgbs_20_40.chunk.npz"
                np.savez_compressed(npz_file_name_20_40, X = X_batch_20_40, Y = Y_batch_20_40)
                processed_20_40_cpg_num = 0
                del(X_batch_20_40)
                del(Y_batch_20_40)
        # ----------------------------------------------
        # wgbs: 40-60
        #
        Y_batch_40_60_indices = np.where( (Y_batch > 0.4) & (Y_batch <= 0.6) )
        if len(Y_batch_40_60_indices) > 0:
            if ('X_batch_40_60' not in globals()) and ('Y_batch_40_60' not in globals()):
                X_batch_40_60 = X_batch[Y_batch_40_60_indices]
                Y_batch_40_60 = Y_batch[Y_batch_40_60_indices]
                processed_40_60_cpg_num = processed_40_60_cpg_num + len(Y_batch_40_60)
            else:
                x = X_batch[Y_batch_40_60_indices]
                y = Y_batch[Y_batch_40_60_indices]
                X_batch_40_60 = np.concatenate((X_batch_40_60, x), axis = 0)
                Y_batch_40_60 = np.concatenate((Y_batch_40_60, y), axis = 0)
                processed_40_60_cpg_num = processed_40_60_cpg_num + len(y)
            if processed_40_60_cpg_num >= maxCpGlimit:
                FileIndex_40_60 = FileIndex_40_60 + 1
                npz_file_name_40_60 = outputfile + "." + str(FileIndex_40_60) + ".wgbs_40_60.chunk.npz"
                np.savez_compressed(npz_file_name_40_60, X = X_batch_40_60, Y = Y_batch_40_60)
                processed_40_60_cpg_num = 0
                del(X_batch_40_60)
                del(Y_batch_40_60)
        # ----------------------------------------------
        # wgbs: 60-80
        #
        Y_batch_60_80_indices = np.where( (Y_batch > 0.6) & (Y_batch <= 0.8) )
        if len(Y_batch_60_80_indices) > 0:
            if ('X_batch_60_80' not in globals()) and ('Y_batch_60_80' not in globals()):
                X_batch_60_80 = X_batch[Y_batch_60_80_indices]
                Y_batch_60_80 = Y_batch[Y_batch_60_80_indices]
                processed_60_80_cpg_num = processed_60_80_cpg_num + len(Y_batch_60_80)
            else:
                x = X_batch[Y_batch_60_80_indices]
                y = Y_batch[Y_batch_60_80_indices]
                X_batch_60_80 = np.concatenate((X_batch_60_80, x), axis = 0)
                Y_batch_60_80 = np.concatenate((Y_batch_60_80, y), axis = 0)
                processed_60_80_cpg_num = processed_60_80_cpg_num + len(y)
            if processed_60_80_cpg_num >= maxCpGlimit:
                FileIndex_60_80 = FileIndex_60_80 + 1
                npz_file_name_60_80 = outputfile + "." + str(FileIndex_60_80) + ".wgbs_60_80.chunk.npz"
                np.savez_compressed(npz_file_name_60_80, X = X_batch_60_80, Y = Y_batch_60_80)
                processed_60_80_cpg_num = 0
                del(X_batch_60_80)
                del(Y_batch_60_80)
        # ----------------------------------------------
        # wgbs: 80-100
        #
        Y_batch_80_100_indices = np.where( (Y_batch > 0.8) & (Y_batch <= 1.0) )
        if len(Y_batch_80_100_indices) > 0:
            if ('X_batch_80_100' not in globals()) and ('Y_batch_80_100' not in globals()):
                X_batch_80_100 = X_batch[Y_batch_80_100_indices]
                Y_batch_80_100 = Y_batch[Y_batch_80_100_indices]
                processed_80_100_cpg_num = processed_80_100_cpg_num + len(Y_batch_80_100)
            else:
                x = X_batch[Y_batch_80_100_indices]
                y = Y_batch[Y_batch_80_100_indices]
                X_batch_80_100 = np.concatenate((X_batch_80_100, x), axis = 0)
                Y_batch_80_100 = np.concatenate((Y_batch_80_100, y), axis = 0)
                processed_80_100_cpg_num = processed_80_100_cpg_num + len(y)
            if processed_80_100_cpg_num >= maxCpGlimit:
                FileIndex_80_100 = FileIndex_80_100 + 1
                npz_file_name_80_100 = outputfile + "." + str(FileIndex_80_100) + ".wgbs_80_100.chunk.npz"
                np.savez_compressed(npz_file_name_80_100, X = X_batch_80_100, Y = Y_batch_80_100)
                processed_80_100_cpg_num = 0
                del(X_batch_80_100)
                del(Y_batch_80_100)
else:
    pass
#---------------------------------------------------------
#---------------------------------------------------------

#
def get_cglist(cglist_file):
    target_cg_list = []
    with open(cglist_file, "rt") as fp:
        for line in fp:
            target_cg_list.append(line.rstrip())
    return target_cg_list
#
#
def normZMWprediction_HCO(ZMWprefile_HiFi, ZMWprefile_ZMW, minZMWdepth=5):
    # 2024-10-15
    ZMWprobDict_Raw = defaultdict(lambda: [[], [], []])
    ZMWprobDict_norm = {}
    processed_pos = set()
    processed_HiFi_zmwid = set()
    # read in Closed_ZMW and Open_ZMW files
    with gzip.open(ZMWprefile_HiFi, "rt") as f_hifi, gzip.open(ZMWprefile_ZMW, "rt") as f_closed_open:
        # HiFi-reads ZMW prediction
        for line in f_hifi:
            if not ("nan" in line) and not line.startswith('Chr'):
                Chr, PosC, ZMWid, ZMWtype, ZMWTdep, Hap, HapBlk, strand, ZMWprob, preclass = line.rstrip("\n").split()
                linekey = f"{Chr}_{PosC}"
                zmwkey = f"{linekey}_{ZMWid}"
                if ZMWtype == "H":
                    ZMWprobDict_Raw[linekey][0].append(float(ZMWprob))
                    processed_pos.add(linekey)
                    processed_HiFi_zmwid.add(zmwkey)
        # Subreads ZMW prediction
        for line in f_closed_open:
            # ZMWtype "C "O"
            if not ("nan" in line) and not line.startswith('Chr'):
                Chr, PosC, ZMWid, ZMWtype, ZMWTdep, Hap, HapBlk, strand, ZMWprob, preclass = line.rstrip("\n").split()
                linekey = f"{Chr}_{PosC}"
                zmwkey = f"{linekey}_{ZMWid}"
                if zmwkey not in processed_HiFi_zmwid: #nonHiFi
                    if ZMWtype == "C":#nonHiFi-ClosedZMW
                        ZMWprobDict_Raw[linekey][1].append(float(ZMWprob))
                        processed_pos.add(linekey)
                    if ZMWtype == "O":#nonHiFi-OpenZMW
                        ZMWprobDict_Raw[linekey][2].append(float(ZMWprob))
                        processed_pos.add(linekey)
    #    
    for pos in ZMWprobDict_Raw:
        if len(ZMWprobDict_Raw[pos][0]) + len(ZMWprobDict_Raw[pos][1]) + len(ZMWprobDict_Raw[pos][2]) >= minZMWdepth:
            HiFi_raw_prob_list, Closed_raw_prob_list, Open_raw_prob_list = ZMWprobDict_Raw[pos]
            # -------
            # count-based prediction as a feature
            total_molecules = len(ZMWprobDict_Raw[pos][0]) + len(ZMWprobDict_Raw[pos][1]) + len(ZMWprobDict_Raw[pos][2])
            me_molecules = np.sum((np.array(ZMWprobDict_Raw[pos][0])>0.5).astype(int)) + np.sum((np.array(ZMWprobDict_Raw[pos][1])>0.5).astype(int)) + np.sum((np.array(ZMWprobDict_Raw[pos][2])>0.5).astype(int))
            count_based_DNAm_freq = me_molecules/total_molecules
            #
            HiFi_norm_prob_list = get_normalized_histo_chatGPT(HiFi_raw_prob_list)
            Closed_norm_prob_list = get_normalized_histo_chatGPT(Closed_raw_prob_list)
            Open_norm_prob_list = get_normalized_histo_chatGPT(Open_raw_prob_list)
            # -------
            # HiFi_norm + ClosedZMW_norm + OpenZMW_norm + count_based_prediction
            ZMWprobDict_norm[pos] = np.concatenate([HiFi_norm_prob_list, Closed_norm_prob_list, Open_norm_prob_list, np.array([count_based_DNAm_freq])])
    return ZMWprobDict_norm

def get_normalized_histo_chatGPT(probs):
    bin_width = 0.1
    cov = len(probs)
    if cov == 0:
        norm_hist = [0] * 10
    else:
        hist = np.histogram(probs, bins=10, range=[0, 1])[0]
        norm = np.linalg.norm(hist)
        norm_hist = hist / norm
    return norm_hist
#
def intersect_function(ZMWprobDict_norm, trueMeDict):
    to_remove = set(ZMWprobDict_norm) - set(trueMeDict)
    for key in to_remove:
        del ZMWprobDict_norm[key]
    WGBS_label = {}
    for key in ZMWprobDict_norm:
        WGBS_label[key] = trueMeDict[key]
    return ZMWprobDict_norm, WGBS_label
#
def getMethyInfo(trueMefile, cutoffdep):
    with gzip.open(trueMefile, "rt") as fp:
        return {
            f"{div[0]}_{int(div[1])}": float(div[3])
            for div in (
                line.strip().split()
                for line in fp
                if line.strip() and int(line.split()[4]) + int(line.split()[5]) >= cutoffdep
            )
        }
#
def CreateDataMatrix_distance(ZMWprobDict_norm, R):
    """
    update: 2024-10-15: sort
    update: 2023-08-22: add distance
    update: 2022-12-13
    centered on the target site
    input: ZMWprobDict_norm: key:chr_pos = [20 values]
    """
    total_normhistos = []
    total_pos = []
    ZMWprobDict_poslist = list(ZMWprobDict_norm)
    ZMWprobDict_poslist_sorted = sorted(ZMWprobDict_poslist, key=lambda x: int(x.split('_')[1]))
    for pos in ZMWprobDict_poslist_sorted:
        total_normhistos.append(ZMWprobDict_norm[pos])
        # -----------------------------
        Chr, PosC = pos.split("_")
        total_pos.append(PosC)
        # -----------------------------
    featPad = np.stack(total_normhistos)
    featuresWindow = sliding_window_view(featPad, 2*R+1, axis=0)
    N = featuresWindow.shape[0]
    # -----------------------------
    featDis=np.stack(total_pos)
    distance_win = sliding_window_view(featDis, 2*R+1, axis=0)
    # -----------------------------
    # Calculate the centered distances for each window
    centered_distances = []
    for i in range(0,N):
        window_string = distance_win[i]
        window = window_string.astype(int)
        centered_value = window[len(window) // 2]
        distances = np.abs(window - centered_value)
        #
        distances_v1 = 100.0 / (np.sqrt(distances) + 100.0)
        distances_v2 = (distances_v1 - 0.5) / 0.5
        centered_distances.append(distances_v2.tolist())  # Convert numpy array to list
    # merge
    #
    # featuresWindow_merge = np.random.random((N, 31, 2*R+1))
    # -------
    # 2023-08-28: update as follows:
    # 32 = HiFi_norm(10) + ClosedZMW_norm(10) + OpenZMW_norm(10) + Count_based_prediction(1) + distance_weights(1)
    M = featuresWindow.shape[1]# HCO=31 feature
    featuresWindow_merge = np.random.random((N, (M+1), 2*R+1))
    #
    for i in range(0,N):
        a = np.array(centered_distances[i])
        b = featuresWindow[i]
        featuresWindow_merge[i] = np.vstack((b, a))
    #
    featuresWindow_merge = np.swapaxes(featuresWindow_merge, 1, 2)
    return featuresWindow_merge, ZMWprobDict_poslist_sorted[R:len(ZMWprobDict_poslist)-R]
