import gzip
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from keras.models import model_from_json
from numpy.lib.stride_tricks import sliding_window_view
from sitepredict import *
#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hifizmw", "--HiFiZMWpre", help="HiFi reads - single molecule level prediction. File format: txt.gz")
    parser.add_argument("-allzmw", "--ALLZMWpre", default="NA", help="ZMW prediction - all molecules' predictions. File format: txt.gz")
    parser.add_argument("-minZMWdep", "--minZMWdepth",help="min ZMW_molecule/HiFi_reads/CCS_reads depth")
    parser.add_argument("-r","--windowRadius", help="how many CpG sites in a sliding window? (2 * r + 1)")
    parser.add_argument("-j", "--modeljson", help="json file of the model, json file")
    parser.add_argument("-w", "--modelweights", help="weights of the model, h5 file")
    parser.add_argument("-o", "--outputfile", help="CpG site level methylation frequency")
    return parser.parse_args()
#
if __name__ == "__main__":
    # --------
    args = get_args()
    #
    ZMWprefile_HiFi_raw = args.HiFiZMWpre
    ALLZMW_prefile = args.ALLZMWpre
    minZMWdepth = int(args.minZMWdepth)
    R = int(args.windowRadius)
    modeljsonFile = args.modeljson
    modelweightsFile = args.modelweights
    outputfile = args.outputfile
    # ------------------------------------------------------------------------------------------
    # step1: model loading
    json_file = open(modeljsonFile, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    loaded_model.load_weights(modelweightsFile)
    # ------------------------------------------------------------------------------------------
    # step2: read in read-level prediction
    ZMWprobDict_norm = {}
    ZMWprobDict_norm, predictions_count = normZMWprediction_HCO(ZMWprefile_HiFi, ALLZMW_prefile, minZMWdepth)
    X_test, X_ave_weights, ZMWprobDict_poslist_sorted = CreateDataMatrix_distance_testing(ZMWprobDict_norm, R)
    # ------------------------------------------------------------------------------------------
    # step3: do the prediction
    predictions_model = loaded_model.predict(X_test)
    # ------------------------------------------------------------------------------------------
    # step4: format to write out
    Chr_info = [x.split("_")[0] for x in ZMWprobDict_poslist_sorted]
    Pos_info = [x.split("_")[1] for x in ZMWprobDict_poslist_sorted]
    #
    predictions_count_trim = [predictions_count[x][0] for x in ZMWprobDict_poslist_sorted]
    predictions_count_trim_me_count = [predictions_count[x][1] for x in ZMWprobDict_poslist_sorted]
    predictions_count_trim_ume_count = [predictions_count[x[2]] for x in ZMWprobDict_poslist_sorted]
    #
    df_SiteLevelPre = pd.DataFrame({"Chr":Chr_info,
                                    "PosC":Pos_info,
                                    "model_based":predictions_model.flatten(),
                                    "count_based":predictions_count_trim,
                                    "Me":predictions_count_trim_me_count,
                                    "unMeth":predictions_count_trim_ume_count,
                                    "cisCpG_ave_weights":X_ave_weights.flatten(),
        })
    df_SiteLevelPre = df_SiteLevelPre[["Chr", "PosC", "model_based","count_based", "Me", "unMeth", "cisCpG_ave_weights"]]
    prefile = outputfile + ".JointCpG.pre.txt.gz"
    df_SiteLevelPre.to_csv(prefile, sep='\t', index=False, compression='gzip')
    # ------------------------------------------------------------------------------------------
else:
    pass
#
#
def normZMWprediction_HCO(ZMWprefile_HiFi,  ALLZMW_prefile, minZMWdepth):
    """
    Freeze,2024-01-19
        INPUT: HiFi read level prediction
               nonHiFi ZMW level prediction
        Return: ZMWprobDict_norm = [ [ HiFi reads ] [ nonHiFi-subreads-closedZMW ] [ nonHiFi-subreads-openZMW ] ]
        Notes: ALLZMW_prefile could be "NA", need to check before selecting one file mode or two file mode
    """
    # ZMWprobDict_Raw = {} #[ [ HiFi reads ] [ nonHiFi-subreads-closedZMW ] [ nonHiFi-subreads-openZMW ] ]
    ZMWprobDict_Raw = defaultdict(lambda: [[], [], []])
    ZMWprobDict_norm = {}
    count_based_pred_dict = {}
    #
    processed_pos = set()
    #
    processed_HiFi_zmwid = set()
    #
    with gzip.open(ZMWprefile_HiFi, "rt") as f_hifi, gzip.open(ALLZMW_prefile, "rt") as f_closed_open:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for line in f_hifi:
            if not ("nan" in line) and not line.startswith('Chr'):
                Chr, PosC, ZMWid, ZMWTdep, ZMWtype, Hap, HapBlk, strand, ZMWprob, preclass = line.rstrip("\n").split()
                linekey = f"{Chr}_{PosC}"
                zmwkey = f"{linekey}_{ZMWid}"
                ZMWprobDict_Raw[linekey][0].append(float(ZMWprob))
                processed_pos.add(linekey)
                processed_HiFi_zmwid.add(zmwkey)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for line in f_closed_open:
            if not ("nan" in line) and not line.startswith('Chr'):
                # ZMWtype "H" "C "O"
                Chr, PosC, ZMWid, ZMWTdep, ZMWtype, Hap, HapBlk, strand, ZMWprob, preclass = line.rstrip("\n").split()
                linekey = f"{Chr}_{PosC}"
                zmwkey = f"{linekey}_{ZMWid}"
                if zmwkey not in processed_HiFi_zmwid:
                    if ZMWtype == "C":
                        ZMWprobDict_Raw[linekey][1].append(float(ZMWprob))
                        processed_pos.add(linekey)
                    if ZMWtype == "O":
                        ZMWprobDict_Raw[linekey][2].append(float(ZMWprob))
                        processed_pos.add(linekey)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for pos in list(ZMWprobDict_Raw):
        # CpG site with the desired ZMW to be covered
        # ZMWprobDict_Raw = {} #[ [ HiFi reads ] [ nonHiFi-subreads-closedZMW ] [ nonHiFi-subreads-openZMW ] ]
        if len(ZMWprobDict_Raw[pos][0]) + len(ZMWprobDict_Raw[pos][1]) + len(ZMWprobDict_Raw[pos][2]) >= int(minZMWdepth):
            #
            HiFi_raw_prob_list, Closed_raw_prob_list, Open_raw_prob_list = ZMWprobDict_Raw[pos]
            # -------
            # count-based prediction as a feature
            total_molecules = len(ZMWprobDict_Raw[pos][0]) + len(ZMWprobDict_Raw[pos][1]) + len(ZMWprobDict_Raw[pos][2])
            me_molecules = np.sum((np.array(ZMWprobDict_Raw[pos][0])>0.5).astype(int)) + np.sum((np.array(ZMWprobDict_Raw[pos][1])>0.5).astype(int)) + np.sum((np.array(ZMWprobDict_Raw[pos][2])>0.5).astype(int))
            count_based_DNAm_freq = me_molecules/total_molecules
            count_based_pred_dict[pos] = [count_based_DNAm_freq, me_molecules, total_molecules-me_molecules]
            #
            HiFi_norm_prob_list = get_normalized_histo(HiFi_raw_prob_list)# len=10
            Closed_norm_prob_list = get_normalized_histo(Closed_raw_prob_list)# len=10
            Open_norm_prob_list = get_normalized_histo(Open_raw_prob_list)# len=10
            # -------
            ZMWprobDict_norm[pos] = np.concatenate([HiFi_norm_prob_list, Closed_norm_prob_list, Open_norm_prob_list, np.array([count_based_DNAm_freq])])
    return ZMWprobDict_norm, count_based_pred_dict
#
#
def get_normalized_histo(probs):
    """
    Freeze 2024-01-09
    """
    cov = len(probs)
    if cov == 0:
        norm_hist = [0] * 10
    else:
        hist = np.histogram(probs, bins=10, range=[0, 1])[0]
        norm = np.linalg.norm(hist)
        norm_hist = hist / norm
    return norm_hist
#
#
def CreateDataMatrix_distance_testing(ZMWprobDict_norm, R):
    """
    Freeze 2024-01-19:
        Input: ZMWprobDict_norm = [ [ HiFi reads ] [ nonHiFi-subreads-closedZMW ] [ nonHiFi-subreads-openZMW ] ]
        Return:
    """
    total_normhistos = []
    total_pos = []
    total_ave_weights = []
    ZMWprobDict_poslist = list(ZMWprobDict_norm)
    ZMWprobDict_poslist_sorted = sorted(ZMWprobDict_poslist, key=lambda x: int(x.split('_')[1]))
    #
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
    featDis = np.stack(total_pos)
    distance_win = sliding_window_view(featDis, 2*R+1, axis=0)
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
    M = featuresWindow.shape[1]
    featuresWindow_merge = np.random.random((N, (M+1), 2*R+1))
    #
    for i in range(0,N):
        a = np.array(centered_distances[i])
        # -----------------------------------
        ave_a = np.mean(a)
        total_ave_weights.append(ave_a)
        b = featuresWindow[i]
        featuresWindow_merge[i] = np.vstack((b, a))
    #
    featuresWindow_merge = np.swapaxes(featuresWindow_merge, 1, 2)
    return featuresWindow_merge, np.array(total_ave_weights), ZMWprobDict_poslist_sorted[R:len(ZMWprobDict_poslist)-R]
