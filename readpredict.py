import time
import pysam
import random
import array
import argparse
import numpy as np
from readpredict import *
from collections import defaultdict
from keras.models import model_from_json
#
#
def get_args():
    parser = argparse.ArgumentParser()
    # input related parameters
    parser.add_argument("-bam", "--bamfile", help="aligned BAM file with kinetic signals, fi/fp/ri/rp for hifi reads, and ip/pw for subreads")
    parser.add_argument("-cg", "--cpgbed", help="CpG position, column1 is chr, column2 is pos, such as chr1 132. Here pos is 1-based coordinate on Rreference's forward strand")
    parser.add_argument("-type", "--datatype", default='HIFI', help="data mode, SUBREAD or HIFI")
    parser.add_argument("-r", "--radius", default='10',type=int, help="cpg +/- r")
    parser.add_argument("-min", "--minsubreadsdepth", default='1', type=int, help="min subreads depth for a ZMW to be included")
    parser.add_argument("-max", "--maxsubreadsdepth", default='60', type=int, help="max subreads depth for a ZMW to be included")
    parser.add_argument("-maq", "--MappingQuality", default='30', type=int, help="min mapping quality for a read to be included")
    parser.add_argument("-bpq", "--BaseQuality", default='0',type=int, help="min base quality required")
    parser.add_argument("-chr", "--Chromosome", default="whole_genome", help="processing chromosome")
    parser.add_argument("-j", "--modeljson", help="json file of the model, json file")
    parser.add_argument("-w", "--modelweights", help="weights of the model, h5 file")
    # processing parameters
    parser.add_argument("-p", "--Processor", type=int,help="Parallelism processors")
    # output related parameters
    parser.add_argument("-sc", "--Scale", action='store_true',default=False, help="flag to indicate scale the signal of each subreads, (x-mean)/std. Without this tag is raw signal value")
    parser.add_argument("-oe", "--outputType", help = "Output Types: ModBam,Feature,MoleculeLevel. ModBam: Modification Bam file. | Feature: Feature matrix. | MoleculeLevel: Molecule level modification prediction")
    parser.add_argument("-o", "--outprefix", help = "[*outprefix*].modification.bam for ModBam outputType.  [*outprefix*].Kmat.txt.gz for Feature outputType. [*outprefix*].SingleMol.pre.txt.gz for MoleculeLevel outputType")
    parser.add_argument("-keepK","--KeepKinetics", default="remove", help="remove or keep. flag to indicate keep or remove the kinetic sigals in the output Bam file")
    #
    # TODOlater other parameters:2023-12-26
    # (1) CPU numbers
    # (2) openZMW
    # (3) models
    #
    return parser.parse_args()
#  
#
if __name__ == "__main__":
    # testing file:
    # bamfile="/storage/yangjianLab/caoyujie/project/meth/Bam_file/4_hifi_subreads/HG01109_WT_hifi/speed_testing_debug/chr22_5x.sort.HG01109_WT_hifi.aln.bam"
    # cgfile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/pbMC/cpos/hg38/chr/chr22.cpos1.list"
    # modelweightsFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/4_ancesty/7_train_by_all/2_ZMW_model/ZMW_model_5_inds.json.weights.h5"
    # modeljsonFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/4_ancesty/7_train_by_all/2_ZMW_model/ZMW_model_5_inds.json"
    args = get_args()
    #
    bamf = args.bamfile
    cgfile = args.cpgbed
    datatype = args.datatype
    R = int(args.radius)
    minsubdep = int(args.minsubreadsdepth)
    maxsubdep = int(args.maxsubreadsdepth)
    mappingQ = int(args.MappingQuality)   
    baseQ = int(args.BaseQuality)
    processingChr = args.Chromosome
    modeljsonFile = args.modeljson
    modelweightsFile = args.modelweights
    #
    proc = int(args.Processor)
    #
    scaleFlag = args.Scale
    outputType = args.outputType
    outputp = args.outprefix    
    KiFlag = args.KeepKinetics

    if outputType == "Feature":
        # Action: ZMW level feature extraction
        outfile = outputp + ".Kmat.txt.gz"
        pass
    else:
        # ZMWmeth model loading
        json_file = open(modeljsonFile, 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights(modelweightsFile)
        # cpg reference position loading
        cglist = Buildcglist_dict(cgfile, processingChr)
        # single CPU version, 2023-12-26
        HiFiReadsMeth(inbamf, cglist, outputp, outputType, minsubdep, maxsubdep, baseQ, mappingQ, R, scaleFlag, KiFlag, loaded_model)
else:
    pass
#
#
#
# Coordinates in pysam are always 0-based (following the python convention)
def HiFiReadsMeth(inbamf, cglist, outputp, outputType, minsubdep, maxsubdep, baseQ, mappingQ, R, scaleFlag, KiFlag, loaded_model):
    """
       2023-12-26
       INPUT:  aligned BAM file with fi/fp/ri/rp tags
       OUTPUT: when outputType == "ModBam": aligned BAM file with modification status coding in MM and ML tags, with/without fi/fp/ri/rp tags(KiFlag="move","keep")
               when outputType == "MoleculeLevel": Sinle molecule level prediction
    """
    samfile = pysam.AlignmentFile(inbamf, "rb")
    if outputType == "ModBam":
        # Action: modification coded in ML/MM tag of the output BAM file
        outputf_bam_file = outputp + ".modification.bam"
        outsamfile = pysam.AlignmentFile(outputf_bam_file, "wb", template=samfile)
    if outputType == "MoleculeLevel":
        # Action: Sinle molecule level prediction
        outfile = outputp + ".SingleMol.pre.txt.gz"
        pass
    winsize = 2 * R + 1
    output_reads = []# store certain number of reads before writing out
    #
    for read in samfile.fetch():
        if (not read.is_secondary) and (not read.is_supplementary) and (read.mapping_quality >= mappingQ) and (Matching_ratio_Fun(read) >= 0.85):
            fn = read.get_tag("fn")
            rn = read.get_tag("rn")
            total_subreads_dep = int(fn) + int(rn)
            if (total_subreads_dep >= minsubdep) and (total_subreads_dep <= maxsubdep) and (int(fn) >= 1) and (int(rn) >= 1):
                aln_ref_chr = read.reference_name
                aln_ref_start = read.reference_start #0-based leftmost coordinate
                aln_ref_end = read.reference_end
                overlapping_cpg = find_overlapping_cpg(cglist_sort_increasing = cglist[aln_ref_chr], Ref_chr = aln_ref_chr, Ref_start = aln_ref_start, Ref_end = aln_ref_end)
                #
                if len(overlapping_cpg) >= 1:
                    # ------------------------------------------------------------------------------------------------
                    # step1: get the alignment information for the read
                    read_fi_list = read.get_tag("fi")
                    read_fp_list = read.get_tag("fp")
                    read_ri_list = read.get_tag("ri")
                    read_rp_list = read.get_tag("rp")
                    read_SEQ_string = read.query_sequence
                    read_is_reverse = read.is_reverse
                    read_query_length = read.query_length
                    read_cigartuples = read.cigartuples
                    haplotype, haplotype_block = (str(read.get_tag("HP")), str(read.get_tag("PS"))) if read.has_tag("HP") else ("X", "X")# HP:i:2   PS:i:11457
                    # ------------------------------------------------------------------------------------------------
                    # step2: locate genome CpG to each read
                    # return: ["cpg1":pos1,"cpg2":pos2 ... ]
                    cpgPos_on_read = locate_cpgPos_on_read(aln_ref_chr = aln_ref_chr, aln_ref_start = aln_ref_start, read_cigar = read_cigartuples, overlapping_cpg = overlapping_cpg)
                    if len(cpgPos_on_read) >= 1:
                        X_reads = np.empty((len(cpgPos_on_read), 13, 2*R+1))
                        all_cpgs = list(cpgPos_on_read)
                        for index in range(len(cpgPos_on_read)):
                            cpg = all_cpgs[index]
                            Chr, pos_C = cpg.split("_")
                            pos_on_read = cpgPos_on_read[cpg]
                            #
                            # collect read level features
                            F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list = HiFiRead_cpg_K_Feature(pos_on_read, read_is_reverse, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag)
                            #
                            # construct the feature matrix
                            X_mat_cpg = Format_closedZMW(F_template_seq = F_template_seq, F_ipdlist = F_ipd_list, F_pwlist = F_pw_list, R_template_seq = R_template_seq, R_ipdlist = R_ipd_list, R_pwlist = R_pw_list, npasses = total_subreads_dep)
                            #
                            # add X_mat_cpg to the X_reads
                            X_reads[index] = X_mat_cpg
                        # -----------------------------------
                        # CpGs on a read finished, next, (1) predict methy_prob and create ML tag; (2) create MM tag; (3) iterate each tag, keeping or removing fi/fp/ri/rp tags; (4) create new alignment enquiry and write out
                        # -----------------------------------
                       
                        # (1)
                        # start_time = time.time()    
                        MLtag_value = read_Meth_prediction(X_reads, loaded_model)
                        # end_time = time.time()
                        # execution_time = end_time - start_time
                        # print(execution_time)
                        # MLtag_tuple = ("ML", MLtag_value, "B")
                        MLtag_tuple = ("ML", MLtag_value)
                        #
                        # (2)
                        leadingC_pos = cpgPos_on_read[all_cpgs[0]]
                        seq_0_lC = read_SEQ_string[0:leadingC_pos]
                        count_0_lC = seq_0_lC.count('C')
                        MMtag_value = str(count_0_lC)
                        for c1, c2 in zip(list(cpgPos_on_read.keys())[::2], list(cpgPos_on_read.keys())[1::2]):
                            c1_pos_on_read = cpgPos_on_read[c1]
                            c2_pos_on_read = cpgPos_on_read[c2]
                            seq_c1_c2 = read_SEQ_string[c1_pos_on_read:c2_pos_on_read]
                            count_c1_c2 = seq_c1_c2.count('C')
                            MMtag_value = MMtag_value + "," + str(count_c1_c2)
                        # MMtag_tuple = ("MM","C+m,"+MMtag_value+";",'Z')
                        MMtag_tuple = ("MM","C+m,"+MMtag_value+";")
                        #
                        # (3)
                        if  KiFlag == "remove":
                            read.set_tag('fi', None)
                            read.set_tag('fp', None)
                            read.set_tag('ri', None)
                            read.set_tag('rp', None)
                        read_full_tags_list = read.get_tags(with_value_type = False)
                        read_full_tags_list.append(MMtag_tuple)
                        read_full_tags_list.append(MLtag_tuple)
                        read.set_tags(read_full_tags_list)
                        #
                        # (4) write out
                        output_reads.append(read)
                        if len(output_reads) >= 10000:
                            for processed_read in output_reads:
                                outsamfile.write(processed_read)
                            output_reads.clear()
    # Write out the remaining lines if any
    if output_reads:
        for processed_read in output_reads:
            outsamfile.write(processed_read)
    return 0
#
#
def getliststat(inputlist):
    return (np.array(inputlist)).min(),(np.array(inputlist)).mean(),(np.array(inputlist)).std()
#
#
def getkineticswin(inputKlist, scaleFlag=True):
    if scaleFlag:
        Kmin, Kmean, Kstd = getliststat(inputKlist)
        Kwinlist = np.round(((np.array(inputKlist)) - Kmean) / Kstd, 2)
        return Kwinlist
    else:
        return inputKlist
#
#
# 2023-12-11: read
def Matching_ratio_Fun(read):
    TotalLen = 0
    Matching_len = 0
    cigar = read.cigartuples
    for Tuple in cigar:
        op, count = Tuple
        # 0:M, 1:I, 4:S, 7:=, 8:X ||||||||||||| sum of "M/I/S/=/X" is the SEQ length
        if ((op == 0) or (op == 1) or (op == 4) or (op == 7) or (op == 8)):
            TotalLen = TotalLen + count
            if ((op == 0) or (op == 7) or (op == 8)):
                Matching_len = Matching_len + count
    Matching_ratio = float(Matching_len / TotalLen)
    return Matching_ratio
#
#
def Buildcglist_dict(cgfile, processingChr):
    cglist_dict = defaultdict(list)
    with open(cgfile, "r") as fp:
        for line in fp:
            div = line.rstrip("\n").split()
            if processingChr == "whole_genome" or div[0] == processingChr:
                chromosome = div[0]
                value = div[0] + "_" + div[1]
                cglist_dict[chromosome].append(value)
    # Sort values numerically based on div[1]
    for chromosome in cglist_dict:
        cglist_dict[chromosome] = sorted(cglist_dict[chromosome], key=lambda x: int(x.split("_")[1]))
    return cglist_dict
#
#
def find_overlapping_cpg(cglist_sort_increasing, Ref_chr, Ref_start, Ref_end):
    chromosome_dict = {}
    overlapping_cpg = []
    # Create dictionary for quick chromosome lookup
    for cpg in cglist_sort_increasing:
        chromosome, pos_C = cpg.split("_")
        if chromosome not in chromosome_dict:
            chromosome_dict[chromosome] = []
        chromosome_dict[chromosome].append(cpg)
    # Binary search to find the first index where pos_C > Ref_end
    chromosomes = chromosome_dict.get(Ref_chr, [])
    left = 0
    right = len(chromosomes) - 1
    while left <= right:
        mid = (left + right) // 2
        pos_C = int(chromosomes[mid].split("_")[1])
        if pos_C <= Ref_end:
            left = mid + 1
        else:
            right = mid - 1
    # Iterate over the CpG sites for the specific chromosome
    for cpg in chromosomes[:left]:
        chromosome, pos_C = cpg.split("_")
        if Ref_start <= int(pos_C) <= Ref_end:
            overlapping_cpg.append(cpg)
        elif int(pos_C) > Ref_end:
            break
    return overlapping_cpg
#
#
# needs to simplify later, 2023-12-26
def locate_cpgPos_on_read(aln_ref_chr, aln_ref_start, read_cigar, overlapping_cpg):
    """
    FUN:locate ref_cpg_position on the read by decoding the cigar_string
    INPUT:overlapping_cpg = ["chr_C1","chr_C2","chr_C3" ... ]
    RETURN:cpgPos_on_read_dict = ["chr_C1":pos1,"chr_C2":pos2 ... ] on this read
    """
    cpgPos_on_read_dict = {}
    cigar_blk_num = len(read_cigar)
    # ------------------------------------------------------------------------
    for cpg in overlapping_cpg:
        read_leading_flag = 1
        read_pos_blk_start = 0
        read_pos_blk_end = 0
        cpg_ref_pos = int(cpg.split("_")[1])# 1-based position
        for index in range(cigar_blk_num):
            op, count = read_cigar[index]
            if (op == 1 or op == 4 or op == 5 or op == 6) and (read_leading_flag == 1):
                if not (op==6):#1,4 consume reads, 5_hard_clip should also be counted on reads, since the aligner doesn't change ipd/pw signal vector as the SEQ
                    read_pos_blk_start = read_pos_blk_end
                    read_pos_blk_end = read_pos_blk_start + count
            else:
                if (read_leading_flag == 1):
                    read_leading_flag = 2
                    if (op==0 or op==2 or op==3 or op==7 or op==8) and (read_leading_flag == 2):
                        Ref_walking_pos_start = aln_ref_start
                        Ref_walking_pos_end = Ref_walking_pos_start + count
                    if (op == 0 or op == 1  or op == 4 or op == 7 or op == 8) and (read_leading_flag == 2):
                        read_pos_blk_start = read_pos_blk_end
                        read_pos_blk_end =  read_pos_blk_start + count
                else:
                    if (op==0 or op==2 or op==3 or op==7 or op==8):
                        Ref_walking_pos_start = Ref_walking_pos_end
                        Ref_walking_pos_end = Ref_walking_pos_start + count
                    if (op == 0 or op == 1  or op == 4 or op == 7 or op == 8):
                        read_pos_blk_start = read_pos_blk_end
                        read_pos_blk_end =  read_pos_blk_start + count
                # capture the cpg on the moving ref
                if (Ref_walking_pos_start <= cpg_ref_pos) and (Ref_walking_pos_end >= cpg_ref_pos):
                    if (not (op == 3)) and (not (op == 2)):
                        last_op_needed = cpg_ref_pos - Ref_walking_pos_start
                        if (op==0 or op==1 or op==4 or op==7 or op==8):
                            read_pos = read_pos_blk_start + last_op_needed
                            cpgPos_on_read_dict[cpg] = read_pos - 1
                            break
    return cpgPos_on_read_dict
#
#
# 2023-12-25:update
def Format_closedZMW(F_template_seq, F_ipdlist, F_pwlist, R_template_seq, R_ipdlist, R_pwlist,npasses):
    F_seq_context = list(F_template_seq)
    F_ipdlist = np.asarray(F_ipdlist, dtype=np.float16)
    F_pwlist = np.asarray(F_pwlist, dtype=np.float16)
    #
    R_seq_context = list(R_template_seq)
    R_ipdlist = np.asarray(R_ipdlist, dtype=np.float16)
    R_pwlist = np.asarray(R_pwlist, dtype=np.float16)
    #
    total_win_size = len(F_ipdlist)  # R=10, total_win_size = 21
    npasses_list = np.full(total_win_size, npasses)
    #
    F_seq_mat_array = np.zeros((4, total_win_size))
    R_seq_mat_array = np.zeros((4, total_win_size))
    #
    seq_dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for index in range(total_win_size):
        F_seq_mat_array[seq_dict[F_seq_context[index]], index] = 1
        R_seq_mat_array[seq_dict[R_seq_context[index]], index] = 1
    mat = np.vstack([F_seq_mat_array, F_ipdlist, F_pwlist, R_seq_mat_array, R_ipdlist, R_pwlist, npasses_list])
    return mat
#
#
def str_complementary(string_a):
    translation_table = str.maketrans("ATCG", "TAGC")
    complementary_string_a = string_a.translate(translation_table)
    complementary_string_a = "".join(
        random.choice("ATCG") if c not in "ATCG" else c
        for c in complementary_string_a
    )
    return complementary_string_a
#
#
def string_checker(string_a):
    list_string_a = list(string_a)
    for index in range(len(list_string_a)):
        if list_string_a[index] not in ['A', 'T', 'C', 'G']:
            list_string_a[index] = str(
                random.sample(['A', 'T', 'C', 'G'], 1)[0])
            continue
    checkedstring = "".join(list_string_a)
    return checkedstring
#
#
def read_Meth_prediction(X_reads, loaded_model):
    """return np.array for each CpG on the reads"""
    X_testing = np.transpose(X_reads, (0, 2, 1))
    predictions = loaded_model.predict(X_testing)
    predictions_list = array.array('i',(predictions.flatten() * 256).astype(int))
    return predictions_list
#
#
def HiFiRead_cpg_K_Feature(pos_on_read, read_is_reverse, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag):
    """
    2023-12-25:
    Related knowledge: DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
    Return: (1) sequence context from the template of Ref_Forward_strand(should be in 3'->5', on reads is 5'->3') + IPD/PW on the reads while syning the template(should be in F_ipd_list and F_pw_list)
            (2) sequence context from the template of Ref_Reverse_strand + IPD/PW on the reads while syning the template(should be in R_ipd_list and R_pw_list)
    """
    if read_is_reverse:# this hifi read aligned to Ref_R_strand
        # (1)
        # For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
        # information is stored in fi, fp tags
        RefFC_on_f_list = read_query_length - pos_on_read - 1#checked, 2023-12-25
        L_RefFC_on_f_list = RefFC_on_f_list - R - 1
        R_RefFC_on_f_list = RefFC_on_f_list + R
        if (L_RefFC_on_f_list > 0) and (R_RefFC_on_f_list < read_query_length):
            F_ipd_list = getkineticswin(inputKlist=read_fi_list[L_RefFC_on_f_list:R_RefFC_on_f_list],scaleFlag=scaleFlag)
            F_pw_list = getkineticswin(inputKlist=read_fp_list[L_RefFC_on_f_list:R_RefFC_on_f_list],scaleFlag=scaleFlag)
            # 2023-12-25:F_template_seq is the template sequence context in 3'->5' direction, which is the moving direction of the DNA polyermase
            F_template_seq = read_SEQ_string[::-1][L_RefFC_on_f_list:R_RefFC_on_f_list]
            #
        # (2)
        # For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
        # information is stored in ri, rp tags
        RefRC_on_r_list = pos_on_read + 1
        L_RefRC_on_r_list = RefRC_on_r_list - R - 1
        R_RefRC_on_r_list = RefRC_on_r_list + R
        if (L_RefRC_on_r_list > 0) and (R_RefRC_on_r_list < read_query_length):
            R_ipd_list = getkineticswin(inputKlist=read_ri_list[L_RefRC_on_r_list:R_RefRC_on_r_list],scaleFlag=scaleFlag)
            R_pw_list = getkineticswin(inputKlist=read_rp_list[L_RefRC_on_r_list:R_RefRC_on_r_list],scaleFlag=scaleFlag)
            R_template_seq = (str_complementary(read_SEQ_string))[L_RefRC_on_r_list:R_RefRC_on_r_list]
            R_ipd_list = R_ipd_list[::-1]
            R_pw_list = R_pw_list[::-1]

    else:# this hifi read aligned to Ref_F_strand
        # (1)
        # For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
        # information is stored in fi, fp tags
        RefRC_on_f_list = pos_on_read + 1
        L_RefRC_on_f_list = RefRC_on_f_list - R - 1
        R_RefRC_on_f_list = RefRC_on_f_list + R
        if (L_RefRC_on_f_list>0) and (R_RefRC_on_f_list<read_query_length):
            R_ipd_list = getkineticswin(inputKlist=read_fi_list[L_RefRC_on_f_list:R_RefRC_on_f_list],scaleFlag=scaleFlag)
            R_pw_list = getkineticswin(inputKlist=read_fp_list[L_RefRC_on_f_list:R_RefRC_on_f_list],scaleFlag=scaleFlag)
            R_template_seq = (str_complementary(read_SEQ_string))[L_RefRC_on_f_list:R_RefRC_on_f_list]
            #
        # (2)
        # For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
        # information is stored in ri, rp tags
        RefFC_on_r_list = read_query_length - pos_on_read - 1
        L_RefFC_on_r_list = RefFC_on_r_list - R - 1
        R_RefFC_on_r_list = RefFC_on_r_list + R
        if (L_RefFC_on_r_list>0) and (R_RefFC_on_r_list<read_query_length):
            F_ipd_list = getkineticswin(inputKlist=read_ri_list[L_RefFC_on_r_list:R_RefFC_on_r_list],scaleFlag=scaleFlag)
            F_pw_list = getkineticswin(inputKlist=read_rp_list[L_RefFC_on_r_list:R_RefFC_on_r_list],scaleFlag=scaleFlag)
            F_template_seq = read_SEQ_string[::-1][L_RefFC_on_r_list:R_RefFC_on_r_list]
            F_ipd_list = F_ipd_list[::-1]
            F_pw_list = F_pw_list[::-1]
    return F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list
#
#
def Buildcglist_dict(cgfile, processingChr):
    cglist_dict = defaultdict(list)
    with open(cgfile, "r") as fp:
        for line in fp:
            div = line.rstrip("\n").split()
            if processingChr == "whole_genome" or div[0] == processingChr:
                chromosome = div[0]
                value = div[0] + "_" + div[1]
                cglist_dict[chromosome].append(value)
    # Sort values numerically based on div[1]
    for chromosome in cglist_dict:
        cglist_dict[chromosome] = sorted(cglist_dict[chromosome], key=lambda x: int(x.split("_")[1]))
    return cglist_dict
#
#
def Matching_ratio_Fun(read):
    TotalLen = 0
    Matching_len = 0
    cigar = read.cigartuples
    for Tuple in cigar:
        op, count = Tuple
        # 0:M, 1:I, 4:S, 7:=, 8:X ||||||||||||| sum of "M/I/S/=/X" is the SEQ length
        if ((op == 0) or (op == 1) or (op == 4) or (op == 7) or (op == 8)):
            TotalLen = TotalLen + count
            if ((op == 0) or (op == 7) or (op == 8)):
                Matching_len = Matching_len + count
    Matching_ratio = float(Matching_len / TotalLen)
    return Matching_ratio
#
#
def getliststat(inputlist):
    return (np.array(inputlist)).min(),(np.array(inputlist)).mean(),(np.array(inputlist)).std()
#
#
def getkineticswin(inputKlist, scaleFlag=True):
    if scaleFlag:
        Kmin, Kmean, Kstd = getliststat(inputKlist)
        Kwinlist = np.round(((np.array(inputKlist)) - Kmean) / Kstd, 2)
        return Kwinlist
    else:
        return inputKlist