import re
import pysam
import random
import array
import gzip
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import tensorflow as tf
from collections import Counter
import sys
from MeLoDe_Read import *
#
#
def get_args():
    parser = argparse.ArgumentParser()
    # input parameters
    parser.add_argument("-bam", "--bamfile", required=True, help="aligned/unaligned BAM file with kinetic signals, fi/fp/ri/rp for hifi reads, and ip/pw for subreads")
    parser.add_argument("-bamstatus", "--InBamStatus",required=True, help="aligned or unaligned")
    parser.add_argument("-format","--readformat",required=True, help="read type, SUBREAD or HIFI")
    parser.add_argument("-oe", "--outputType", required=True, help = "Output Types: ModBam,Feature,MoleculeLevel. ModBam: Bam file with MM/ML tags. | Feature: Feature matrix. | MoleculeLevel: Molecule level modification prediction")
    parser.add_argument("-o", "--outprefix", required=True, help = "[*outprefix*].modification.bam for ModBam outputType.  [*outprefix*].Kmat.txt.gz for Feature outputType. [*outprefix*].SingleMol_5mCpG.txt.gz for MoleculeLevel outputType")
    # ZMWmeth model files
    parser.add_argument("-H","--HDir", required=False, help="ZMWmeth-HiFi model directory")
    parser.add_argument("-C","--CDir", required=False, help="ZMWmeth-ClosedZMW model directory")
    parser.add_argument("-O","--ODir", required=False, help="ZMWmeth-OpenZMW model directory")
    # control parameters
    parser.add_argument("-keepK","--KeepKinetics", required=False, default="remove", help="remove or keep. flag to indicate keep or remove the kinetic sigals in the output Bam file")
    parser.add_argument("-r", "--radius", required=False, default='10',type=int, help="CpG +/- r base pairs will be included in the read-level prediction")
    parser.add_argument("-min", "--minsubreadsdepth", required=False, default='1', type=int, help="min subreads depth for a ZMW to be included")
    parser.add_argument("-max", "--maxsubreadsdepth", required=False, default='60', type=int, help="max subreads depth for a ZMW to be included")
    parser.add_argument("-maq", "--MappingQuality", required=False, default='30', type=int, help="min mapping quality for a read to be included")
    parser.add_argument("-p", "--Processor", required=False, default=1, type=int,help="Parallelism processors")
    parser.add_argument("-ref","--Reference", required=False, help="The reference genome in fasta format when the input is an aligned bam")
    parser.add_argument("-sc", "--Scale", required=False, action='store_true',default=False, help="flag to indicate scale the signal of each subreads, (x-mean)/std. Without this tag is raw signal value")
    parser.add_argument("-chr", "--Chromosome", required=False, default="whole_genome", help="processing chromosome")
    return parser.parse_args()
#  
#
if __name__ == "__main__":

    args = get_args()
    #
    #
    bamf = args.bamfile
    inputType = args.InBamStatus
    datatype = args.readformat
    outputType = args.outputType
    outputp = args.outprefix
    #
    #
    KiFlag = args.KeepKinetics
    R = int(args.radius)
    minsubdep = int(args.minsubreadsdepth)
    maxsubdep = int(args.maxsubreadsdepth)
    mappingQ = int(args.MappingQuality)
    proc = int(args.Processor)
    scaleFlag = args.Scale
    processingChr = args.Chromosome
    #
    if datatype == "HIFI":
        if inputType == "unaligned":
            if outputType == "Feature":
                outfile = outputp + ".Kmat.txt.gz"
                TOP_HiFiReadsFeature(Topbam=bamf, outfile=outfile, minsubdep=minsubdep, maxsubdep=maxsubdep, R=R, scaleFlag=scaleFlag)
            if outputType == "ModBam":
                if args.HDir:
                    loaded_model = tf.saved_model.load(args.HDir, None)
                    outputfile = outputp + ".5mc.mod.unaligned.bam"
                    TOP_HiFiReadsMeth(Topbam=bamf, modTopbam=outputfile, minsubdep=minsubdep, maxsubdep=maxsubdep, R=R, scaleFlag=scaleFlag, KiFlag=KiFlag, loaded_model=loaded_model)
                else:
                    print("error! --HDir is missing!")
                    sys.exit()
        elif inputType == "aligned":
            if outputType == "Feature":
                if args.Reference:
                    cglist = Buildcglist_dict_Seq(args.Reference)
                    outfile = outputp + ".Kmat.txt.gz"
                    Aligned_HiFiReadsFeature(bamf, cglist, outfile, minsubdep, maxsubdep, mappingQ, R, scaleFlag)
                else:
                    print("error! Please input --Reference file used for the aligned bam file you provided")
                    sys.exit()
            else:
                if outputType == "ModBam" or outputType == "MoleculeLevel":
                    if args.HDir and args.Reference:
                        loaded_model = tf.saved_model.load(args.HDir, None)
                        cglist = Buildcglist_dict_Seq(args.Reference)
                        Aligned_HiFiReadsMeth(bamf, cglist, outputp, outputType, minsubdep, maxsubdep, mappingQ, R, scaleFlag, KiFlag, loaded_model)
                    else:
                        print("error! --HDir is missing, or --Reference is missing, or both are missing")
                        sys.exit()
                else:
                    print("error! outputType should be ModBam or MoleculeLevel")
                    sys.exit()
        else:
            print("--InBamStatus should be aligned or unaligned")
            sys.exit()
    #
    if datatype == "SUBSEAD":
        if inputType == "aligned":
            if outputType == "Feature":
                if args.Reference:
                    cglist = Buildcglist_dict_Seq(args.Reference)
                    outfile = outputp + ".Kmat.txt.gz"
                    Aligned_subreadsBam_Feature(bamf, cglist, outputp, outputType, mappingQ, R, scaleFlag)
                else:
                    print("error! Please input --Reference file used for the aligned bam file you provided")
                    sys.exit()
            if outputType == "MoleculeLevel":
                if args.CDir and args.ODir:
                    # load closedZMW model
                    closedZMW_subreads_model = tf.saved_model.load(args.CDir, None)
                    # load openZMW model
                    openZMW_subreads_model = tf.saved_model.load(args.ODir, None)
                    #
                    Aligned_subreadsBam_Meth(bamf, cglist, outputp, outputType, mappingQ, R, scaleFlag, closedZMW_subreads_model, openZMW_subreads_model)
                else:
                    print("error! Please input ZMWmeth models for closedZMW(--CDir) and openZMW(--ODir)")
                    sys.exit()
        else:
            print("error! Subreads data should be aligned to the reference genome")
            sys.exit()
else:
    pass
#
#
def TOP_HiFiReadsFeature(Topbam, outfile, minsubdep, maxsubdep, R, scaleFlag):
    """
        Freeze,2024-01-04
        INPUT: unaligned BAM file with fi/fp/ri/rp
        OUTPUT: read-level features for each CpG(each line), but without the reference positions
    """
    output_lines = []
    samfile = pysam.AlignmentFile(Topbam, "rb", check_sq=False)
    fout = gzip.open(outfile, "wb")   
    for read in samfile.fetch(until_eof = True):
        fn = read.get_tag("fn")
        rn = read.get_tag("rn")
        total_subreads_dep = int(fn) + int(rn)
        if (total_subreads_dep >= minsubdep) and (total_subreads_dep <= maxsubdep) and (int(fn) >= 1) and (int(rn) >= 1):
            # ------------------------------------------------------------------------------------------------
            # step1: get the information for the read
            read_fi_list = read.get_tag("fi")
            read_fp_list = read.get_tag("fp")
            read_ri_list = read.get_tag("ri")
            read_rp_list = read.get_tag("rp")
            read_SEQ_string = read.query_sequence
            read_query_length = read.query_length
            zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
            read_is_reverse = False
            # ------------------------------------------------------------------------------------------------
            # step2: locate CpG on each read
            # return: [pos1, pos2, pos3, pos4, ... ] in increasing order, pos1, pos2, pos3 are 0-based index, a[pos1] = 'C', a[pos1, pos1+1]='C', a[pos1, pos1+2] = 'CG'
            cpgstring = 'CG'#5'->3' CG context. Template 5'-3'CpG is readed as 5'-3' CpG on reads by the DNA polymerase. Recall that DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
            cpgPos_on_read = [m.start() for m in re.finditer(cpgstring, read_SEQ_string)]
            # ------------------------------------------------------------------------------------------------
            # remove cpg at the head or tail of a read, which causing out of range index
            i = 0
            while i < len(cpgPos_on_read) and cpgPos_on_read[0] < R + 5:
                i += 1
            cpgPos_on_read = cpgPos_on_read[i:]
            i = len(cpgPos_on_read) - 1
            while i >= 0 and cpgPos_on_read[i] > len(read_SEQ_string) - R - 5:
                i -= 1
            cpgPos_on_read = cpgPos_on_read[:i+1]
            # ------------------------------------------------------------------------------------------------
            if len(cpgPos_on_read) >= 1:
                for index in range(len(cpgPos_on_read)):
                    pos_on_read = cpgPos_on_read[index]
                    # collect read level features
                    HiFiread_template_seq, HiFiread_template_ipd_list, HiFiread_template_pw_list, HiFiread_com_template_seq, HiFiread_com_template_ipd_list, HiFiread_com_template_pw_list = HiFiRead_cpg_K_Feature(pos_on_read, read_is_reverse, R, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag)
                    # format line parts
                    F_template_seq_string = "".join(HiFiread_template_seq)
                    F_ipd_part_line = ",".join(map(str, HiFiread_template_ipd_list))
                    F_pw_part_line = ",".join(map(str, HiFiread_template_pw_list))
                    R_template_seq_string = "".join(HiFiread_com_template_seq)
                    R_ipd_part_line = ",".join(map(str, HiFiread_com_template_ipd_list))
                    R_pw_part_line = ",".join(map(str, HiFiread_com_template_pw_list))
                    # combine line parts into a line
                    output_cpg_line = f"Top\t{zmwname}\t{pos_on_read}\t{fn}\t{rn}\t{F_template_seq_string}\t{F_ipd_part_line}\t{F_pw_part_line}\t{R_template_seq_string}\t{R_ipd_part_line}\t{R_pw_part_line}\n"
                    output_lines.append(output_cpg_line)
                    if len(output_lines) >= 10000:
                        fout.write(''.join(output_lines).encode('utf-8'))
                        output_lines.clear()
    # Write out the remaining lines if any
    if output_lines:
        fout.write(''.join(output_lines).encode('utf-8'))
    return 0
#
#
def TOP_HiFiReadsMeth(Topbam, modTopbam, minsubdep, maxsubdep, R, scaleFlag, KiFlag, loaded_model):
    """
    Bug_fix:2024-09-12, including all the reads
    Freeze with bug for ML tag,2024-01-04
        INPUT: Topbam: the original unaligned HiFi Bam file with fi/fp/ri/rp signals
        OUTPUT: modTopbam: modification called bam file tagged by MM,ML. And the fi/fp/ri/rp signals can be kept or removed according to KiFlag.
        Compared with def Aligned_HiFiReadsMeth, the input Topbam file is unmapped/unaligned HiFi reads.
        So, we don't need to decode CIGAR string to locate CpG positions on a read, but using python string function
    """
    output_reads = []
    samfile = pysam.AlignmentFile(Topbam, "rb", check_sq=False)
    outsamfile = pysam.AlignmentFile(modTopbam, "wb", template=samfile)    
    for read in samfile.fetch(until_eof = True):
        fn = read.get_tag("fn")
        rn = read.get_tag("rn")
        total_subreads_dep = int(fn) + int(rn)
        if (total_subreads_dep >= minsubdep) and (total_subreads_dep <= maxsubdep) and (int(fn) >= 1) and (int(rn) >= 1):
            # ------------------------------------------------------------------------------------------------
            # step1: get the information for the read
            read_fi_list = read.get_tag("fi")
            read_fp_list = read.get_tag("fp")
            read_ri_list = read.get_tag("ri")
            read_rp_list = read.get_tag("rp")
            read_SEQ_string = read.query_sequence
            read_query_length = read.query_length
            read_is_reverse = False# checked, 2024-01-03
            # ------------------------------------------------------------------------------------------------
            # step2: locate CpG on each read
            # return: [pos1, pos2, pos3, pos4, ... ] in decending order, pos1, pos2, pos3 are 0-based index, a[pos1] = 'C', a[pos1, pos1+1]='C', a[pos1, pos1+2] = 'CG'
            cpgstring = 'CG'#5'->3' CG context. 2023-12-28: double-checked, template 5'-3'CpG is readed as 5'-3' CpG on reads by the DNA polymerase,Recall that DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
            cpgPos_on_read = [m.start() for m in re.finditer(cpgstring, read_SEQ_string)]
            # ------------------------------------------------------------------------------------------------
            # remove cpg at the head or tail of a read, which causing out of range index
            i = 0
            while i < len(cpgPos_on_read) and cpgPos_on_read[0] < R + 5:
                i += 1
            cpgPos_on_read = cpgPos_on_read[i:]
            i = len(cpgPos_on_read) - 1
            while i >= 0 and cpgPos_on_read[i] > len(read_SEQ_string) - R - 5:
                i -= 1
            cpgPos_on_read = cpgPos_on_read[:i+1]
            # ------------------------------------------------------------------------------------------------
            if len(cpgPos_on_read) >= 1:
                X_reads = np.empty((len(cpgPos_on_read), 13, 2*R+1))
                for index in range(len(cpgPos_on_read)):
                    pos_on_read = cpgPos_on_read[index]
                    # collect read level features
                    HiFiread_template_seq, HiFiread_template_ipd_list, HiFiread_template_pw_list, HiFiread_com_template_seq, HiFiread_com_template_ipd_list, HiFiread_com_template_pw_list = HiFiRead_cpg_K_Feature(pos_on_read, read_is_reverse, R, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag)
                    # construct the feature matrix
                    X_mat_cpg = Format_closedZMW(F_template_seq = HiFiread_template_seq, F_ipdlist = HiFiread_template_ipd_list, F_pwlist = HiFiread_template_pw_list, R_template_seq = HiFiread_com_template_seq, R_ipdlist = HiFiread_com_template_ipd_list, R_pwlist = HiFiread_com_template_pw_list, npasses = total_subreads_dep)
                    # add X_mat_cpg to the X_reads
                    X_reads[index] = X_mat_cpg
                predictions_list = read_Meth_prediction(X_reads, loaded_model)
                #
                # ------------------------------------------------------------------------------------------------
                # (1) ML:B:C,scaled-probabilities, for example: 'ML:B:C,204,89,26,130', currently is 'ML:B:i,204,89,26,130'
                MLtag_value = array.array('i',(predictions_list * 256).astype(int))
                MLtag_tuple = ("ML", MLtag_value)
                # (2) MM tag
                MMtag_value = create_MMtag_value_selected_motif(top_strand = read_SEQ_string, Selected_motif_pos_on_top_strand = cpgPos_on_read, umodifiedbase='C')
                MMtag_tuple = ("MM","C+m,"+MMtag_value+";")
                # (3) kinetics tag keep or remove
                if  KiFlag == "remove":
                    read.set_tag('fi', None)
                    read.set_tag('fp', None)
                    read.set_tag('ri', None)
                    read.set_tag('rp', None)
                read_full_tags_list = read.get_tags(with_value_type = False)
                read_full_tags_list.append(MMtag_tuple)
                read_full_tags_list.append(MLtag_tuple)
                read.set_tags(read_full_tags_list)
                # (4) write out
                output_reads.append(read)
                if len(output_reads) >= 10000:
                    for processed_read in output_reads:
                        outsamfile.write(processed_read)
                    output_reads.clear()
            
            # ------------------------------------------------------------------------------------------------
            # these are reads without CpG sites captured
            else:
                # remove MM,ML tag, in case there is an MM,ML tag already
                read.set_tag('MM', None)
                read.set_tag('ML', None)
                if  KiFlag == "remove":
                    read.set_tag('fi', None)
                    read.set_tag('fp', None)
                    read.set_tag('ri', None)
                    read.set_tag('rp', None)
                read_full_tags_list = read.get_tags(with_value_type = False)
                read.set_tags(read_full_tags_list)
                # (4B) write out
                output_reads.append(read)
                if len(output_reads) >= 10000:
                    for processed_read in output_reads:
                        outsamfile.write(processed_read)
                    output_reads.clear()
        # ------------------------------------------------------------------------------------------------
        # these are reads failed the condition
        else:
            # remove MM,ML tag, in case there is an MM,ML tag already
            read.set_tag('MM', None)
            read.set_tag('ML', None)
            if  KiFlag == "remove":
                read.set_tag('fi', None)
                read.set_tag('fp', None)
                read.set_tag('ri', None)
                read.set_tag('rp', None)
            read_full_tags_list = read.get_tags(with_value_type = False)
            read.set_tags(read_full_tags_list)
            # (4B) write out
            output_reads.append(read)
            if len(output_reads) >= 10000:
                for processed_read in output_reads:
                    outsamfile.write(processed_read)
                output_reads.clear()

    # Write out the remaining lines if any
    if output_reads:
        for processed_read in output_reads:
            outsamfile.write(processed_read)
    outsamfile.close()
    return 0
#
#
#
def Aligned_HiFiReadsFeature(inbamf, cglist, outfile, minsubdep, maxsubdep, mappingQ, R, scaleFlag):
    """
        INPUT: aligned BAM file with fi/fp/ri/rp
        OUTPUT: read-level features for each CpG with the reference positions
    """
    samfile = pysam.AlignmentFile(inbamf, "rb")
    fout = gzip.open(outfile, "wb")
    output_lines = []
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
                    zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
                    # ------------------------------------------------------------------------------------------------
                    # step2: locate genome CpG to each read
                    # return: ["cpg1":pos1,"cpg2":pos2 ... ]
                    cpgPos_on_read = locate_cpgPos_on_SEQ(aln_ref_chr = aln_ref_chr, aln_ref_start = aln_ref_start, read_cigar = read_cigartuples, overlapping_cpg = overlapping_cpg)
                    # ----------------------------------------------------------------------------------------
                    # heading or tailing removing
                    cpgPos_on_read_keys = list(cpgPos_on_read)
                    for key in cpgPos_on_read_keys:
                        if cpgPos_on_read[key] < (R+5):
                            del cpgPos_on_read[key]
                        else:
                            break
                    for key in reversed(cpgPos_on_read_keys):
                        if cpgPos_on_read[key] > len(read_SEQ_string) - R - 5:
                            del cpgPos_on_read[key]
                        else:
                            break
                    # ----------------------------------------------------------------------------------------
                    if len(cpgPos_on_read) >= 1:
                        all_cpgs = list(cpgPos_on_read)
                        for index in range(len(cpgPos_on_read)):
                            cpg = all_cpgs[index]
                            Chr, pos_C = cpg.split("_")
                            pos_on_read = cpgPos_on_read[cpg]
                            #
                            # collect read level features
                            F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list = HiFiRead_cpg_K_Feature(pos_on_read, read_is_reverse, R, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag)
                            #
                            F_template_seq_string = "".join(F_template_seq)
                            F_ipd_part_line = ",".join(map(str, F_ipd_list))
                            F_pw_part_line = ",".join(map(str, F_pw_list))
                            R_template_seq_string = "".join(R_template_seq)
                            R_ipd_part_line = ",".join(map(str, R_ipd_list))
                            R_pw_part_line = ",".join(map(str, R_pw_list))
                            output_cpg_line = f"{Chr}\t{pos_C}\t{zmwname}\t{fn}\t{rn}\t{F_template_seq_string}\t{F_ipd_part_line}\t{F_pw_part_line}\t{R_template_seq_string}\t{R_ipd_part_line}\t{R_pw_part_line}\n"
                            output_lines.append(output_cpg_line)
                            if len(output_lines) >= 10000:
                                fout.write(''.join(output_lines).encode('utf-8'))
                                output_lines.clear()
    # Write out the remaining lines if any
    if output_lines:
        fout.write(''.join(output_lines).encode('utf-8'))
    return 0
#
#
#
# Coordinates in pysam are always 0-based (following the python convention)
def Aligned_HiFiReadsMeth(bamf, cglist, outputp, outputType, minsubdep, maxsubdep, mappingQ, R, scaleFlag, KiFlag, loaded_model):
    """
    Freeze with bug for ML tag,2024-01-04,2024-01-05,2024-01-09
        INPUT:  aligned BAM file with fi/fp/ri/rp tags
        OUTPUT: when outputType == "ModBam": aligned BAM file with modification status coding in MM and ML tags, with/without fi/fp/ri/rp tags(KiFlag="move","keep")
               when outputType == "MoleculeLevel": Single molecule level prediction
    """
    samfile = pysam.AlignmentFile(bamf, "rb")
    output_reads = []
    #
    # ---------
    if outputType == "ModBam":
        # Action: modification coded in ML/MM tag of the output BAM file
        outputf_bam_file = outputp + ".modification.bam"
        outsamfile = pysam.AlignmentFile(outputf_bam_file, "wb", template=samfile)
    if outputType == "MoleculeLevel":
        # Action: Sinle molecule level prediction
        outfile = outputp + ".SingleMol.pre.txt.gz"
        fout = gzip.open(outfile, "wb")
    # ---------
    #
    # prediction process, by read by cpg.
    for read in samfile.fetch():
        fn = read.get_tag("fn")
        rn = read.get_tag("rn")
        total_subreads_dep = int(fn) + int(rn)
        if (not read.is_secondary) and (not read.is_supplementary) and (read.mapping_quality >= mappingQ) and (Matching_ratio_Fun(read) >= 0.85) and (total_subreads_dep >= minsubdep) and (total_subreads_dep <= maxsubdep) and (int(fn) >= 1) and (int(rn) >= 1):
            #
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
                zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
                zmwType = "H"# Hifi
                # ------------------------------------------------------------------------------------------------
                # step2: locate genome CpG to each read
                # return: ["cpg1":pos1,"cpg2":pos2 ... ]
                cpgPos_on_SEQ = locate_cpgPos_on_SEQ(aln_ref_chr = aln_ref_chr, aln_ref_start = aln_ref_start, read_cigar = read_cigartuples, overlapping_cpg = overlapping_cpg)
                # ----------------------------------------------------------------------------------------
                # heading or tailing removing
                cpgPos_on_SEQ_keys = list(cpgPos_on_SEQ)
                for key in cpgPos_on_SEQ_keys:
                    if cpgPos_on_SEQ[key] < (R+5):
                        del cpgPos_on_SEQ[key]
                    else:
                        break
                for key in reversed(cpgPos_on_SEQ_keys):
                    if cpgPos_on_SEQ[key] > len(read_SEQ_string) - R - 5:
                        del cpgPos_on_SEQ[key]
                    else:
                        break
                # cpgPos_on_top_strand(5'->3'), still in the order of ["cpg1":top_pos1, "cpg2":top_pos2,]. such as ['chr1_100':200, 'chr1_1111':20]
                cpgPos_on_top_strand = Con_posSEQ_to_posRead(read_is_reverse = read_is_reverse, read_length = read_query_length, pos_on_SEQ_dict = cpgPos_on_SEQ)
                # ----------------------------------------------------------------------------------------                          
                if len(cpgPos_on_SEQ) >= 1:
                    X_reads = np.empty((len(cpgPos_on_SEQ), 13, 2*R+1))
                    all_cpgs = list(cpgPos_on_SEQ)
                    for index in range(len(cpgPos_on_SEQ)):
                        cpg = all_cpgs[index]
                        pos_on_SEQ = cpgPos_on_SEQ[cpg]
                        # collect read level features
                        F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list = HiFiRead_cpg_K_Feature(pos_on_SEQ, read_is_reverse, R, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag)
                        # construct the feature matrix
                        X_mat_cpg = Format_closedZMW(F_template_seq = F_template_seq, F_ipdlist = F_ipd_list, F_pwlist = F_pw_list, R_template_seq = R_template_seq, R_ipdlist = R_ipd_list, R_pwlist = R_pw_list, npasses = total_subreads_dep)
                        # add X_mat_cpg to the X_reads
                        X_reads[index] = X_mat_cpg
                    #
                    # make predictions for the current read
                    predictions_list = read_Meth_prediction(X_reads, loaded_model)
                    ###########################
                    #### write out process ####
                    ###########################
                    if outputType == "MoleculeLevel":
                        Chr_info, Pos_info = zip(*[cpg.split("_") for cpg in all_cpgs])
                        common_info = [zmwname, zmwType, total_subreads_dep, haplotype, haplotype_block, '+']
                        ZMW_info, zmwType, zmwdep_info, Hap_info, HapBlk_info, strand_sign = [[info] * len(cpgPos_on_SEQ) for info in common_info]
                        predictions_bi = (predictions_list >= 0.5).astype(int)
                        df_testing_pre = pd.DataFrame(
                            {
                                "Chr": Chr_info,
                                "Pos1": Pos_info,
                                "ZMWid": ZMW_info,
                                "ZMWtype": zmwType,
                                "subdepth": zmwdep_info,
                                "Hap": Hap_info,
                                "HapBlk": HapBlk_info,
                                "strand": strand_sign,
                                "y_pred": predictions_list,
                                "y_pred_bi": predictions_bi
                            }
                        )
                        df_testing_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                    if outputType == "ModBam":
                        # 2024-01-05: need to distinguish +- strand, since MM,ML tag is with respected to top_strand(5'->3')
                        if read_is_reverse:
                            predictions_list = predictions_list[::-1]
                            top_strand = str_complementary(read_SEQ_string[::-1])
                            cpgPos_on_top_strand_rev = dict(reversed(cpgPos_on_top_strand.items()))#reverse the key order, such as ['chr1_100':200, 'chr1_1111':20] to ['chr1_1111':20, 'chr1_100':200]
                        else:
                            top_strand = read_SEQ_string
                            cpgPos_on_top_strand_rev = cpgPos_on_top_strand
                        # (1)
                        MLtag_value = array.array('i', (predictions_list * 256).astype(int))
                        MLtag_tuple = ("ML", MLtag_value)
                        # (2)
                        MMtag_value = create_MMtag_value_selected_motif(top_strand, Selected_motif_pos_on_top_strand=cpgPos_on_top_strand_rev, umodifiedbase='C')
                        MMtag_tuple = ("MM", "C+m," + MMtag_value + ";")
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
                # ---------------------------------------------------------------------------------------------------------------
                # these are aligned_reads failed the condition
                else:
                    if outputType == "ModBam":
                        # remove MM,ML tag, in case there is an MM,ML tag already
                        read.set_tag('MM', None)
                        read.set_tag('ML', None)
                        if  KiFlag == "remove":
                            read.set_tag('fi', None)
                            read.set_tag('fp', None)
                            read.set_tag('ri', None)
                            read.set_tag('rp', None)
                        read_full_tags_list = read.get_tags(with_value_type = False)
                        read.set_tags(read_full_tags_list)
                        # (4B) write out
                        output_reads.append(read)
                        if len(output_reads) >= 10000:
                            for processed_read in output_reads:
                                outsamfile.write(processed_read)
                            output_reads.clear()
            # ---------------------------------------------------------------------------------------------------------------
            # these are aligned_reads failed the condition
            else:
                if outputType == "ModBam":
                    # remove MM,ML tag, in case there is an MM,ML tag already
                    read.set_tag('MM', None)
                    read.set_tag('ML', None)
                    if  KiFlag == "remove":
                        read.set_tag('fi', None)
                        read.set_tag('fp', None)
                        read.set_tag('ri', None)
                        read.set_tag('rp', None)
                    read_full_tags_list = read.get_tags(with_value_type = False)
                    read.set_tags(read_full_tags_list)
                    # (4B) write out
                    output_reads.append(read)
                    if len(output_reads) >= 10000:
                        for processed_read in output_reads:
                            outsamfile.write(processed_read)
                        output_reads.clear()
        # ---------------------------------------------------------------------------------------------------------------
        # these are aligned_reads failed the condition
        else:
            if outputType == "ModBam":
                # remove MM,ML tag, in case there is an MM,ML tag already
                read.set_tag('MM', None)
                read.set_tag('ML', None)
                if  KiFlag == "remove":
                    read.set_tag('fi', None)
                    read.set_tag('fp', None)
                    read.set_tag('ri', None)
                    read.set_tag('rp', None)
                read_full_tags_list = read.get_tags(with_value_type = False)
                read.set_tags(read_full_tags_list)
                # (4B) write out
                output_reads.append(read)
                if len(output_reads) >= 10000:
                    for processed_read in output_reads:
                        outsamfile.write(processed_read)
                    output_reads.clear()
    return 0

#
#
def find_overlapping_cpg(cglist_sort_increasing, Ref_chr, Ref_start, Ref_end):
    """Freeze,2024-01-04"""
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
def locate_cpgPos_on_SEQ(aln_ref_chr, aln_ref_start, read_cigar, overlapping_cpg):
    """
    Freeze,2024-01-04
    FUN:locate ref_cpg_position on the read by decoding the cigar_string
    INPUT:overlapping_cpg = ["chr_C1","chr_C2","chr_C3" ... ]
    RETURN:cpgPos_on_read_dict = ["chr_C1":pos1,"chr_C2":pos2 ... ] on this read. |-Here,pos1,pos2...are 0-based on the read-|
    # ---
    Notes: 2024-01-02: Acutally, this is position on SEQ. Function name changed accordingly
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
def Format_closedZMW(F_template_seq, F_ipdlist, F_pwlist, R_template_seq, R_ipdlist, R_pwlist, npasses):
    """Freeze,2024-01-04"""
    F_seq_context = list(F_template_seq)
    F_ipdlist = np.asarray(F_ipdlist, dtype=np.float16)
    F_pwlist = np.asarray(F_pwlist, dtype=np.float16)
    #
    R_seq_context = list(R_template_seq)
    R_ipdlist = np.asarray(R_ipdlist, dtype=np.float16)
    R_pwlist = np.asarray(R_pwlist, dtype=np.float16)
    #
    total_win_size = len(F_ipdlist)
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
def Format_openZMW(template_seq, ipdlist, pwlist, npasses):
    """Freeze, 2024-01-17"""
    seq_context = list(template_seq)
    ipdlist = np.asarray(ipdlist, dtype=np.float16)
    pwlist = np.asarray(pwlist, dtype=np.float16)
    #
    total_win_size = len(ipdlist)
    npasses_list = np.full(total_win_size, npasses)
    #
    seq_mat_array = np.zeros((4, total_win_size))
    #
    seq_dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for index in range(total_win_size):
        seq_mat_array[seq_dict[seq_context[index]], index] = 1
    mat = np.vstack([seq_mat_array, ipdlist, pwlist, npasses_list])
    return mat
#
def Buildcglist_dict(cgfile, processingChr="whole_genome"):
    """Freeze,2024-01-04"""
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
def Matching_ratio_Fun(read):
    """Freeze,2024-01-04"""
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
def getkineticswin(inputKlist, scaleFlag=True):
    """Freeze,2024-01-04"""
    if scaleFlag:
        Kmin, Kmean, Kstd = getliststat(inputKlist)
        Kwinlist = np.round(((np.array(inputKlist)) - Kmean) / Kstd, 2)
        return Kwinlist
    else:
        return inputKlist
#
def getliststat(inputlist):
    """Freeze,2024-01-04"""
    return (np.array(inputlist)).min(),(np.array(inputlist)).mean(),(np.array(inputlist)).std()
#
def str_complementary(string_a):
    """Freeze,2024-01-04"""
    complementary_string_a = string_a.translate(string_a.maketrans("ATCG", "TAGC"))
    return complementary_string_a
#
def read_Meth_prediction(X_reads, loaded_model):
    """
    Modified, 2024-10-09, create predictions_list, Done
    Modified,2024-07-23,notes: the golang version should be correct, no update is needed
    Freeze,2024-01-04
    """
    X_testing = np.transpose(X_reads, (0, 2, 1))
    X_test_tensor = tf.convert_to_tensor(X_testing)
    X_test_tensor = tf.cast(X_test_tensor, tf.float32)
    inference_func = loaded_model.signatures["serving_default"]
    predictions = inference_func(input_1=X_test_tensor)
    key = list(predictions)[0]

    predictions_array = predictions[key].numpy()
    predictions_list = predictions_array.flatten()

    return predictions_list
#
#
def HiFiRead_cpg_K_Feature(pos_on_SEQ, read_is_reverse, R, read_query_length, read_SEQ_string, read_fi_list, read_fp_list, read_ri_list, read_rp_list, scaleFlag):
    """
    Freeze, 2024-01-09
    Related knowledge: DNA polymerase moves along in 3'->5' direction on the template, generating the 5'->3' reads.
    Return: (1) TopHiFi_Template_seq in 3'->5' + IPD/PW while syning the TopHiFi_Template_seq
                HiFiread_template_seq, HiFiread_template_ipd_list, HiFiread_template_pw_list
            (2) ComHiFi_Template_seq in 3'->5' + IPD/PW while syning the ComHiFi_Template_seq
                HiFiread_com_template_seq, HiFiread_com_template_ipd_list, HiFiread_com_template_pw_list
     """
    if read_is_reverse:# this hifi read aligned to Ref_R_strand
        # (1)
        # For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
        # information is stored in fi, fp tags, when synthesis the topHiFi
        RefFC_on_f_list = read_query_length - pos_on_SEQ - 1#checked, 2023-12-25,2023-12-28
        L_RefFC_on_f_list = RefFC_on_f_list - R - 1
        R_RefFC_on_f_list = RefFC_on_f_list + R
        HiFiread_template_ipd_list = getkineticswin(inputKlist=read_fi_list[L_RefFC_on_f_list:R_RefFC_on_f_list],scaleFlag=scaleFlag)
        HiFiread_template_pw_list = getkineticswin(inputKlist=read_fp_list[L_RefFC_on_f_list:R_RefFC_on_f_list],scaleFlag=scaleFlag)
        HiFiread_template_seq = (read_SEQ_string[pos_on_SEQ-R:pos_on_SEQ+R+1])[::-1]
        # (2)
        # For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
        # information is stored in ri, rp tags, when synthesis the comHiFi
        RefRC_on_r_list = pos_on_SEQ + 1
        L_RefRC_on_r_list = RefRC_on_r_list - R
        R_RefRC_on_r_list = RefRC_on_r_list + R + 1
        HiFiread_com_template_ipd_list = getkineticswin(inputKlist=read_ri_list[L_RefRC_on_r_list:R_RefRC_on_r_list],scaleFlag=scaleFlag)
        HiFiread_com_template_pw_list = getkineticswin(inputKlist=read_rp_list[L_RefRC_on_r_list:R_RefRC_on_r_list],scaleFlag=scaleFlag)
        HiFiread_com_template_seq = (str_complementary(read_SEQ_string))[L_RefRC_on_r_list:R_RefRC_on_r_list]
        HiFiread_com_template_ipd_list = HiFiread_com_template_ipd_list[::-1]
        HiFiread_com_template_pw_list = HiFiread_com_template_pw_list[::-1]
    
    else:# this hifi read aligned to Ref_F_strand, topHiFi=SEQ
        # (1)
        # For subreads mapped to Ref_F_strand, which carrying information when syn the complementary R_strand, key 'R' #
        # information is stored in fi, fp tags, when synthesis the topHiFi
        RefRC_on_f_list = pos_on_SEQ + 1
        L_RefRC_on_f_list = RefRC_on_f_list - R
        R_RefRC_on_f_list = RefRC_on_f_list + R + 1
        HiFiread_template_ipd_list = getkineticswin(inputKlist=read_fi_list[L_RefRC_on_f_list:R_RefRC_on_f_list],scaleFlag=scaleFlag)
        HiFiread_template_pw_list = getkineticswin(inputKlist=read_fp_list[L_RefRC_on_f_list:R_RefRC_on_f_list],scaleFlag=scaleFlag)
        HiFiread_template_seq = (str_complementary(read_SEQ_string))[L_RefRC_on_f_list:R_RefRC_on_f_list]
        # (2)
        # For subreads mapped to Ref_R_strand, which carrying information when syn the complementary F_strand, key 'F' #
        # information is stored in ri, rp tags, when synthesis the comHiFi
        RefFC_on_r_list = read_query_length - pos_on_SEQ - 1
        L_RefFC_on_r_list = RefFC_on_r_list - R - 1
        R_RefFC_on_r_list = RefFC_on_r_list + R
        HiFiread_com_template_ipd_list = getkineticswin(inputKlist=read_ri_list[L_RefFC_on_r_list:R_RefFC_on_r_list],scaleFlag=scaleFlag)
        HiFiread_com_template_pw_list = getkineticswin(inputKlist=read_rp_list[L_RefFC_on_r_list:R_RefFC_on_r_list],scaleFlag=scaleFlag)
        HiFiread_com_template_seq = (read_SEQ_string[pos_on_SEQ-R:pos_on_SEQ+R+1])[::-1]
        HiFiread_com_template_ipd_list = HiFiread_com_template_ipd_list[::-1]
        HiFiread_com_template_pw_list = HiFiread_com_template_pw_list[::-1]
    return HiFiread_template_seq, HiFiread_template_ipd_list, HiFiread_template_pw_list, HiFiread_com_template_seq, HiFiread_com_template_ipd_list, HiFiread_com_template_pw_list
#
#
def Aligned_subreadsBam_Feature(bamf, cglist, outputp, outputType, mappingQ, R, scaleFlag):
    """
        2024-01-17,2024-01-16,2024-01-03
        INPUT: aligned subreads BAM file with ip/pw tags
        OUTPUT: outputType == "Feature": Single molecule level prediction
    """
    samfile = pysam.AlignmentFile(bamf, "rb")
    chromosome_order_list = list(samfile.references)
    Feature_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    if outputType == "Feature":
        outfile = outputp + ".SingleMol.KFeature.txt.gz"
        fout = gzip.open(outfile, "wb")
        processed_read_count = 0
        for read in samfile.fetch():
            if (not read.is_secondary) and (not read.is_supplementary) and (read.mapping_quality >= mappingQ) and (Matching_ratio_Fun(read) >= 0.85):
                aln_ref_chr = read.reference_name
                aln_ref_start = read.reference_start #0-based leftmost coordinate
                aln_ref_end = read.reference_end
                overlapping_cpg = find_overlapping_cpg(cglist_sort_increasing = cglist[aln_ref_chr], Ref_chr = aln_ref_chr, Ref_start = aln_ref_start, Ref_end = aln_ref_end)
                if len(overlapping_cpg) >= 1:
                    read_ipd_list = read.get_tag("ip")
                    read_pw_list = read.get_tag("pw")
                    read_SEQ_string = read.query_sequence
                    read_is_reverse = read.is_reverse
                    read_query_length = read.query_length
                    read_cigartuples = read.cigartuples
                    read_query_name = read.query_name
                    zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
                    cpgPos_on_read = locate_cpgPos_on_SEQ(aln_ref_chr = aln_ref_chr, aln_ref_start = aln_ref_start, read_cigar = read_cigartuples, overlapping_cpg = overlapping_cpg)
                    # ----------------------------------------------------------------------------------------
                    # heading or tailing removing
                    cpgPos_on_read_keys = list(cpgPos_on_read)
                    for key in cpgPos_on_read_keys:
                        if cpgPos_on_read[key] < (R+5):
                            del cpgPos_on_read[key]
                        else:
                            break
                    for key in reversed(cpgPos_on_read_keys):
                        if cpgPos_on_read[key] > len(read_SEQ_string) - R - 5:
                            del cpgPos_on_read[key]
                        else:
                            break
                    # ----------------------------------------------------------------------------------------
                    if len(cpgPos_on_read) >= 1:
                        all_cpgs = list(cpgPos_on_read)
                        for index in range(len(cpgPos_on_read)):
                            cpg = all_cpgs[index]
                            pos_on_SEQ = cpgPos_on_read[cpg]
                            # collect subread level features, should be raw value
                            strand_sign, seq_string, ipd_list, pw_list = SubRead_cpg_K_Feature(pos_on_SEQ, read_is_reverse, read_query_length, read_SEQ_string, read_ipd_list, read_pw_list)
                            Feature_dict[cpg][zmwname][strand_sign][read_query_name] = [seq_string, ipd_list, pw_list]
                        # ---
                        # after reading certain lines, check whether the safe distance reached    
                        processed_read_count = processed_read_count + 1
                        # ---
                        if processed_read_count % 2000 == 0:
                            # 2024-01-23
                            # > a = ['chr1_40', 'chr1_30', 'chr1_80', 'chr2_10', 'chr1_33', 'chr2_100']
                            # > b = ['chr1', 'chr2', 'chr3']
                            # > sorted_a = sorted(a, key=lambda x: (b.index(x.split('_')[0]), int(x.split('_')[1])))
                            # > sorted_a
                            # > ['chr1_30', 'chr1_33', 'chr1_40', 'chr1_80', 'chr2_10', 'chr2_100']
                            sorted_processed_cpgs = sorted(list(Feature_dict), key=lambda x: (chromosome_order_list.index(x.split('_')[0]), int(x.split('_')[1])))
                            first_C_chr, first_C_pos = sorted_processed_cpgs[0].split("_")
                            last_C_chr, last_C_pos = sorted_processed_cpgs[-1].split("_")
                            if first_C_chr == last_C_chr:
                                if int(last_C_pos) - int(first_C_pos) >= 500000:
                                    # output the top 50% of the CpG in Feature_dict
                                    top_index = int(len(sorted_processed_cpgs) / 2)
                                    top_cpgs = sorted_processed_cpgs[0:top_index]
                                    #
                                    for top_cpg in top_cpgs:
                                        for zmw in Feature_dict[top_cpg]:
                                            ZMWline_by_subreads = Create_ZMWline_by_Subreads(top_cpg, zmw, Feature_dict[top_cpg][zmw], scaleFlag)
                                            fout.write(ZMWline_by_subreads.encode('utf-8'))
                                            processed_read_count = 0
                                    # remove the top cpgs from the Feature_dict
                                    for top_cpg in top_cpgs:
                                        del Feature_dict[top_cpg]
                                    continue
                            else:# different chromosome
                                # should write out elements in the dict, that located at the first_C_chr
                                first_C_chr_cpgs = [ element for element in list(Feature_dict) if element.startswith(first_C_chr) ]
                                for out_cpg in first_C_chr_cpgs:
                                    for zmw in Feature_dict[out_cpg]:
                                        ZMWline_by_subreads = Create_ZMWline_by_Subreads(out_cpg, zmw, Feature_dict[top_cpg][zmw], scaleFlag)
                                        fout.write(ZMWline_by_subreads.encode('utf-8'))
                                        del Feature_dict[top_cpg]
                                processed_read_count = 0
                                continue
                        else:
                            continue
        # ------------------------
        # 2024-01-01: at last, all the things in Feature_dict should be written out
        if Feature_dict:
            for out_cpg in list(Feature_dict):
                for zmw in Feature_dict[out_cpg]:
                    ZMWline_by_subreads = Create_ZMWline_by_Subreads(out_cpg, zmw, Feature_dict[top_cpg][zmw])
                    fout.write(ZMWline_by_subreads.encode('utf-8'))
    else:
        print("error! subreads bam should beoutputType == Feature")
#
#
def Aligned_subreadsBam_Meth(bamf, cglist, outputp, outputType, mappingQ, R, scaleFlag,closedZMW_subreads_model,openZMW_subreads_model):
    """
        2024-01-16,2024-01-19
        INPUT: aligned subreads BAM file with ip/pw tags
        OUTPUT: outputType == "MoleculeLevel": Single molecule level prediction
    """
    samfile = pysam.AlignmentFile(bamf, "rb")
    chromosome_order_list = list(samfile.references)
    Feature_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    if outputType == "MoleculeLevel":
        outfile = outputp + ".SingleMol.pre.txt.gz"
        fout = gzip.open(outfile, "wb")
        processed_read_count = 0
        for read in samfile.fetch():
            if (not read.is_secondary) and (not read.is_supplementary) and (read.mapping_quality >= mappingQ) and (Matching_ratio_Fun(read) >= 0.85):
                aln_ref_chr = read.reference_name
                aln_ref_start = read.reference_start #0-based leftmost coordinate
                aln_ref_end = read.reference_end
                overlapping_cpg = find_overlapping_cpg(cglist_sort_increasing = cglist[aln_ref_chr], Ref_chr = aln_ref_chr, Ref_start = aln_ref_start, Ref_end = aln_ref_end)
                if len(overlapping_cpg) >= 1:
                    read_ipd_list = read.get_tag("ip")
                    read_pw_list = read.get_tag("pw")
                    read_SEQ_string = read.query_sequence
                    read_is_reverse = read.is_reverse
                    read_query_length = read.query_length
                    read_cigartuples = read.cigartuples
                    read_query_name = read.query_name
                    zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
                    cpgPos_on_read = locate_cpgPos_on_SEQ(aln_ref_chr = aln_ref_chr, aln_ref_start = aln_ref_start, read_cigar = read_cigartuples, overlapping_cpg = overlapping_cpg)
                    # ----------------------------------------------------------------------------------------
                    # heading or tailing removing
                    cpgPos_on_read_keys = list(cpgPos_on_read)
                    for key in cpgPos_on_read_keys:
                        if cpgPos_on_read[key] < (R+5):
                            del cpgPos_on_read[key]
                        else:
                            break
                    for key in reversed(cpgPos_on_read_keys):
                        if cpgPos_on_read[key] > len(read_SEQ_string) - R - 5:
                            del cpgPos_on_read[key]
                        else:
                            break
                    # ----------------------------------------------------------------------------------------
                    if len(cpgPos_on_read) >= 1:
                        all_cpgs = list(cpgPos_on_read)
                        for index in range(len(cpgPos_on_read)):
                            cpg = all_cpgs[index]
                            pos_on_SEQ = cpgPos_on_read[cpg]
                            # collect subread level features, should be raw value
                            strand_sign, seq_string, ipd_list, pw_list = SubRead_cpg_K_Feature(pos_on_SEQ, read_is_reverse, read_query_length, read_SEQ_string, read_ipd_list, read_pw_list)
                            Feature_dict[cpg][zmwname][strand_sign][read_query_name] = [seq_string, ipd_list, pw_list]
                        # ---
                        # after reading certain lines, check whether the safe distance reached    
                        processed_read_count = processed_read_count + 1
                        # ---
                        if processed_read_count % 2000 == 0:
                            sorted_processed_cpgs = sorted(list(Feature_dict), key=lambda x: (chromosome_order_list.index(x.split('_')[0]), int(x.split('_')[1])))
                            first_C_chr, first_C_pos = sorted_processed_cpgs[0].split("_")
                            last_C_chr, last_C_pos = sorted_processed_cpgs[-1].split("_")
                            if first_C_chr == last_C_chr:# same chromosome
                                if int(last_C_pos) - int(first_C_pos) >= 500000:
                                    # output the top 50% of the CpG in Feature_dict
                                    top_index = int(len(sorted_processed_cpgs) / 2)
                                    top_cpgs = sorted_processed_cpgs[0:top_index]
                                    # -- making predictions and write out to the file --
                                    for top_cpg in top_cpgs:
                                        Chr, Pos = top_cpg.split("_")
                                        # there could be bug when returning a none
                                        ClosedZMW_X, ClosedZMW_name, ClosedZMW_passes, OpenZMW_X, OpenZMW_name, OpenZMW_passes = Format_zmw_from_subreads(top_cpg, Feature_dict[top_cpg], R, scaleFlag)
                                        i = ClosedZMW_X.shape[0]
                                        j = OpenZMW_X.shape[0]
                                    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                        # (1) both ClosedZMW and OpenZMW are there
                                        if i > 0 and j > 0:
                                            ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                                            OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                                            df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                                            df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                                            merged_df = pd.concat([df_closedzmw_pre, df_openzmw_pre], ignore_index = True)
                                            merged_df.to_csv(fout, sep='\t', index=False, compression='gzip')
                                        # (2) only ClosedZMW are there
                                        if i > 0 and j == 0:
                                            ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                                            df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                                            df_closedzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                                        # (3) only OpenZMW are there
                                        if i == 0 and j > 0:
                                            OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                                            df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                                            df_openzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                                    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                    processed_read_count = 0
                                    # remove the top cpgs from the Feature_dict
                                    for top_cpg in top_cpgs:
                                        del Feature_dict[top_cpg]
                                    continue#continue to read lines from files
                            else:# different chromosome
                                # should write out elements in the dict, that located at the first_C_chr
                                first_C_chr_cpgs = [ element for element in list(Feature_dict) if element.startswith(first_C_chr) ]
                                for out_cpg in first_C_chr_cpgs:
                                    Chr, Pos = out_cpg.split("_")
                                    # there could be bug when returning a none
                                    ClosedZMW_X, ClosedZMW_name, ClosedZMW_passes, OpenZMW_X, OpenZMW_name, OpenZMW_passes = Format_zmw_from_subreads(out_cpg, Feature_dict[out_cpg], R, scaleFlag)
                                    i = ClosedZMW_X.shape[0]
                                    j = OpenZMW_X.shape[0]
                                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                    # (1) both ClosedZMW and OpenZMW are there
                                    if i > 0 and j > 0:
                                        ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                                        OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                                        df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                                        df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                                        merged_df = pd.concat([df_closedzmw_pre, df_openzmw_pre], ignore_index = True)
                                        merged_df.to_csv(fout, sep='\t', index=False, compression='gzip')
                                    # (2) only ClosedZMW are there
                                    if i > 0 and j == 0:
                                        ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                                        df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                                        df_closedzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                                    # (3) only OpenZMW are there
                                    if i == 0 and j > 0:
                                        OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                                        df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                                        df_openzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                for out_cpg in first_C_chr_cpgs:
                                    del Feature_dict[out_cpg]
                                processed_read_count = 0
                                continue
                        else:
                            continue
        # ------------------------
        # 2024-01-01: at last, all the remaining things in Feature_dict should be written out
        if Feature_dict:
            for out_cpg in list(Feature_dict):
                Chr, Pos =  out_cpg.split("_")
                # there could be bug when returning a none
                ClosedZMW_X, ClosedZMW_name, ClosedZMW_passes, OpenZMW_X, OpenZMW_name, OpenZMW_passes = Format_zmw_from_subreads(out_cpg, Feature_dict[out_cpg], R, scaleFlag)
                i = ClosedZMW_X.shape[0]
                j = OpenZMW_X.shape[0]
            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                # (1) both ClosedZMW and OpenZMW are there
                if i > 0 and j > 0:
                    ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                    OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                    df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                    df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                    merged_df = pd.concat([df_closedzmw_pre, df_openzmw_pre], ignore_index = True)
                    merged_df.to_csv(fout, sep='\t', index=False, compression='gzip')
                # (2) only ClosedZMW are there
                if i > 0 and j == 0:
                    ClosedZMW_X_pred = read_Meth_prediction(ClosedZMW_X, closedZMW_subreads_model)
                    df_closedzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = ClosedZMW_name, ZMW_passes = ClosedZMW_passes, ZMW_X_pred = ClosedZMW_X_pred, ZMW_Type='C')
                    df_closedzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
                # (3) only OpenZMW are there
                if i == 0 and j > 0:
                    OpenZMW_X_pred = read_Meth_prediction(OpenZMW_X, openZMW_subreads_model)
                    df_openzmw_pre = format_zmw_data(Chr = Chr, Pos = Pos, ZMW_name = OpenZMW_name, ZMW_passes = OpenZMW_passes, ZMW_X_pred = OpenZMW_X_pred, ZMW_Type='O')
                    df_openzmw_pre.to_csv(fout, sep='\t', index=False, compression='gzip')
            # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        fout.close() 
    else:
        print("error! subreads bam should beoutputType == MoleculeLeve")
#
#
def format_zmw_data(Chr, Pos, ZMW_name, ZMW_passes, ZMW_X_pred, ZMW_Type):
    """
    Freeze, 2024-01-19
    ZMW_Type = 'C' or 'O'
    """
    ZMW_X_pred_bi = (ZMW_X_pred >= 0.5).astype(int)
    zmw_count = len(ZMW_name)
    #
    df_pre = pd.DataFrame({
        "Chr": [Chr] * zmw_count,
        "Pos1": [Pos] * zmw_count,
        "ZMWid": ZMW_name,
        "ZMWtype": [ZMW_Type] * zmw_count,
        "subdepth": ZMW_passes,
        "Hap": ['X'] * zmw_count,
        "HapBlk": ['X'] * zmw_count,
        "strand": ['+/-'] * zmw_count,
        "y_pred": ZMW_X_pred,
        "y_pred_bi": ZMW_X_pred_bi
    })
    return df_pre
#
#
def Format_zmw_from_subreads(top_cpg, Feature_dict_cpg, R, scaleFlag):
    """"
    2024-01-16
    INPUT: Feature_dict_cpg[zmw][strand][read]
    OUTPUT: ClosedZMW_dict=['Info','X'], 'Info'='zmwname'
            OpenZMW_dict=['Info','X'], 'Info'='zmwname'
    """
    ClosedZMW_X = np.empty((len(Feature_dict_cpg), 13, 2*R+1))
    ClosedZMW_name = [0] * len(Feature_dict_cpg)
    ClosedZMW_passes = [0] * len(Feature_dict_cpg)
    i = 0
    #
    OpenZMW_X = np.empty((len(Feature_dict_cpg), 7, 2*R+1))
    OpenZMW_name = [0] * len(Feature_dict_cpg)
    OpenZMW_passes = [0] * len(Feature_dict_cpg)
    j = 0
    #
    for zmw in list(Feature_dict_cpg):
        zmw_features = Feature_dict_cpg[zmw]
        F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list = Sum_Subreads_ZMW(zmw_features, scaleFlag)
        #
        # closedZMW
        if ('R' in Feature_dict_cpg[zmw] and 'F' in Feature_dict_cpg[zmw]):
            npasses = len(Feature_dict_cpg[zmw]['F']) +  len(Feature_dict_cpg[zmw]['R'])
            closezmw_mat = Format_closedZMW(F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list, npasses)
            ClosedZMW_X[i] = closezmw_mat
            ClosedZMW_name[i] = zmw
            ClosedZMW_passes[i] = npasses
            i = i + 1
        #
        # openZMW
        else:
            if F_template_seq == "NA":
                npasses = len(Feature_dict_cpg[zmw]['R'])
                openZMW_mat = Format_openZMW(template_seq = R_template_seq, ipdlist = R_ipd_list, pwlist = R_pw_list, npasses = npasses)
            if R_template_seq == "NA":
                npasses = len(Feature_dict_cpg[zmw]['F'])
                openZMW_mat = Format_openZMW(template_seq = F_template_seq, ipdlist = F_ipd_list, pwlist = F_pw_list, npasses = npasses)
            OpenZMW_X[j] = openZMW_mat
            OpenZMW_name[j] = zmw
            OpenZMW_passes = npasses
            j = j + 1
    return (ClosedZMW_X[0:i,:,:], ClosedZMW_name[0:i,:,:], ClosedZMW_passes[0:i,:,:], OpenZMW_X[0:j,:,:], OpenZMW_name[0:j,:,:], OpenZMW_passes[0:j,:,:])
#
#
def SubRead_cpg_K_Feature(pos_on_SEQ, read_is_reverse, read_query_length, read_SEQ_string, read_ipd_list, read_pw_list):
    """
    Freeze, 2024-01-17
        2024-01-03: subreads keep raw ipd/pw for each of them
        2024-01-02: return Seq_context_on_Template(3'->5'), ipd, pw
    """
    if read_is_reverse:# subreads mapped to Ref(-)_strand, synthesising the Ref(+)_template, ipd/pw are therefor with respected to Ref(+). Template is Ref(+), but should be reverse into 3'-5'
        adjust_pos_on_read = read_query_length - pos_on_SEQ - 1#5'->3' position on read
        L_pos = adjust_pos_on_read - R - 1#5'->3' position on read
        R_pos = adjust_pos_on_read + R#5'->3' position on read
        strand_sign ='F'# 'Ref(+)'#Template strand sign 
        ipd_list = read_ipd_list[L_pos:R_pos]
        pw_list = read_pw_list[L_pos:R_pos]
        # template sequence context in 3'->5' direction.
        seq_string = (read_SEQ_string[pos_on_SEQ-R:pos_on_SEQ+R+1])[::-1]

        #
    else:# subreads mapped to Ref(+), synthesising the Ref(-)_template, ipd/pw are therefor with respected to Ref(-). Template is Ref(-), and should be adjusted to 3'->5'
        adjust_pos_on_read = pos_on_SEQ + 1#5'->3' position on read
        L_pos = adjust_pos_on_read - R#5'->3' position on read
        R_pos = adjust_pos_on_read + R + 1#5'->3' position on read
        strand_sign = 'R'# 'Ref(-)'
        ipd_list = read_ipd_list[L_pos:R_pos]
        pw_list = read_pw_list[L_pos:R_pos]
        seq_string = str_complementary(read_SEQ_string)[pos_on_SEQ-R+1:pos_on_SEQ+R+1+1]
        seq_string = str_complementary(read_SEQ_string)[L_pos:R_pos]
    return strand_sign, seq_string, ipd_list,pw_list
#
#
def Create_ZMWline_by_Subreads(cpg, zmw, zmw_Feature_Dict, scaleFlag=True):
    """
    2024-01-01:
    INPUT:
     - cpg,zmw,Feature_dict[top_cpg][zmw],openZMWFlag
                DO: + strand: take mean and ccs_seq
                    - strand: take mean and ccs_seq
    RETURN:
     - ZMWline
    """
    Chr,pos_C = cpg.split("_")
    zmw_features = zmw_Feature_Dict
    #
    zmwF_subdep = len(zmw_features.get('F', []))
    zmwR_subdep = len(zmw_features.get('R', []))
    #
    F_template_seq, F_ipd_list, F_pw_list, R_template_seq, R_ipd_list, R_pw_list = Sum_Subreads_ZMW(zmw_features,scaleFlag)
    #
    F_template_seq_string = "NA" if F_template_seq == "NA" else "".join(F_template_seq)
    F_ipd_part_line = "NA" if np.array_equal(F_ipd_list, np.array(["NA"])) else ",".join(map(str, F_ipd_list[0]))
    F_pw_part_line = "NA" if np.array_equal(F_pw_list, np.array(["NA"])) else ",".join(map(str, F_pw_list[0]))
    #
    R_template_seq_string = "NA" if R_template_seq == "NA" else "".join(R_template_seq)
    R_ipd_part_line = "NA" if np.array_equal(R_ipd_list, np.array(["NA"])) else ",".join(map(str, R_ipd_list[0]))
    R_pw_part_line = "NA" if np.array_equal(R_pw_list, np.array(["NA"])) else ",".join(map(str, R_pw_list[0]))

    output_cpg_ZMW_line = f"{Chr}\t{pos_C}\t{zmw}\t{zmwF_subdep}\t{zmwR_subdep}\t{F_template_seq_string}\t{F_ipd_part_line}\t{F_pw_part_line}\t{R_template_seq_string}\t{R_ipd_part_line}\t{R_pw_part_line}\n"
    return output_cpg_ZMW_line
#
#
#
def Sum_Subreads_ZMW(ZMW_dict_Feature, scaleFlag=True):
    """
    Freeze, 2024-01-19
    INPUT shape: Feature_dict_cpg[strand][read]
    2024-01-02: both OpenZMW and ClosedZMW should be processed
    INPUT: ZMW_dict_Feature is a dictory, but 'F' or 'R' could not exist, in situations of openZMW

   (1) ClosedZMW:
        'F':{
            'r1':[[0-seq_list][1-ipd_list][2-pw_list]]
            'r2':[[0-seq_list][1-ipd_list][2-pw_list]]
            ....
            }
        'R':{
            'r11':[[0-seq_list][1-ipd_list][2-pw_list]]
            'r22':[[0-seq_list][1-ipd_list][2-pw_list]]
            ....
            }
   (2) OpenZMW with 'F':
        'F':{
            'r1':[[0-seq_list][1-ipd_list][2-pw_list]]
            'r2':[[0-seq_list][1-ipd_list][2-pw_list]]
            ....
            }
    (3) OpenZMW with 'R':
        'R':{
            'r11':[[0-seq_list][1-ipd_list][2-pw_list]]
            'r22':[[0-seq_list][1-ipd_list][2-pw_list]]
            ....
            }
    , where[1-ipd_list] and [2-pw_list] are np.array, [0-seq_list] = ['AAATTTCCCGGG','AAATTTCCCGGG','AAATTTCCCGGG','TAATTTCCCGGG','TAATTTCCCGGG']               
    """
    Returned_zmw_features = ["NA", "NA", "NA", "NA", "NA", "NA"]
    for strand in ZMW_dict_Feature:
        # Initialize an empty list to store the arrays
        ipd_arrays = []
        pw_arrays = []
        seqstring_arrays = []
        #
        # pileup the reads together
        for read in ZMW_dict_Feature[strand]:
            # Create the arrays for the current read and append them to the list
            seqstring_arrays.append(ZMW_dict_Feature[strand][read][0])
            ipd_arrays.append(np.array(ZMW_dict_Feature[strand][read][1]))
            pw_arrays.append(np.array(ZMW_dict_Feature[strand][read][2]))
        #
        # now, all the reads are v-stacked
        out_ccs_seq = get_ccs(seqstring_arrays)
        # notes:  [1,2,3]
        #.........[4,5,6], then mean is [2.5,3.5,4.5] by column
        # if [1,2,3], then keep the original value [1,2,3]          
        #
        if len(ipd_arrays) > 1:
            ccs_ipd_array = np.mean(ipd_arrays, axis = 0, keepdims = True)
            ccs_pw_array = np.mean(pw_arrays, axis = 0, keepdims = True)
        else:
            ccs_ipd_array = ipd_arrays
            ccs_pw_array = pw_arrays
        #
        if scaleFlag:# need to scale the feature
            MEAN_ipd_array = np.mean(ccs_ipd_array)
            MEAN_pw_array = np.mean(ccs_pw_array)
            STD_ipd_array = np.std(MEAN_ipd_array)
            STD_pw_array = np.std(MEAN_pw_array)

            if STD_ipd_array != 0:
                out_ccs_ipd_array = (ccs_ipd_array - MEAN_ipd_array) / STD_ipd_array
            else:
                out_ccs_ipd_array = ccs_ipd_array

            if STD_pw_array != 0:
                out_ccs_pw_array = (ccs_pw_array - MEAN_pw_array) / STD_pw_array
            else:
                out_ccs_pw_array = ccs_pw_array
        #
        else:# don't need to scale the feature, use the raw value
            out_ccs_ipd_array = ccs_ipd_array
            out_ccs_pw_array = ccs_pw_array
        #
        if strand == 'F':
            Returned_zmw_features[0] = out_ccs_seq
            Returned_zmw_features[1] = out_ccs_ipd_array
            Returned_zmw_features[2] = out_ccs_pw_array
        if strand == 'R':
            Returned_zmw_features[3] = out_ccs_seq
            Returned_zmw_features[4] = out_ccs_ipd_array
            Returned_zmw_features[5] = out_ccs_pw_array
        #
    return Returned_zmw_features
#
#
def get_ccs(seqstring_arrays):
    """
    Freeze, 2024-01-17
    seqstring_arrays=['AAATTTCCCGGG','AAATTTCCCGGG','AAATTTCCCGGG','TAATTTCCCGGG','TAATTTCCCGGG']
    code is created by chatCPT
    """
    # Transpose the array to work column-wise
    seqstring_arrays_transposed = list(zip(*seqstring_arrays))
    # Calculate the consensus sequence column-wise
    consensus_sequence = ''
    for column in seqstring_arrays_transposed:
        counts = Counter(column)
        most_common = counts.most_common()
        max_count = most_common[0][1]
        most_common_chars = [char for char, count in most_common if count == max_count]
        selected_char = random.choice(most_common_chars)
        consensus_sequence += selected_char
    return consensus_sequence
#
#
def create_MMtag_value_selected_motif(top_strand, Selected_motif_pos_on_top_strand, umodifiedbase='C'):
    """
    Freeze,2024-01-05
    INPUT: top_strand is the 5'->3' raw read from the sequencing machine.
           Selected_motif_pos_on_top_strand: such as CpG position on top_strand, these are the selected motifs that should be targeted. small positions come first
           umodifiedbase: in CpG situation, it is 'C'
    Please be noted, it is NOT SEQ from alignment file, NOT TEMPLATE STRAND(3'->5',5'->3')
    Below is an example:
        Read  AGTCTAGACTCCGTAATTACTCGCCTAG...
        C        1    2 34       5 6 78
        CpG              *         *

        MM:Z:C+m,3,1,...   # CpG sites are at C #4 (1+3) and #6 (1+3+1+1)
        ML:B:C,249,4,...   # probability of methylation at the first CpG is in [249/256,250/256); second CpG is in [4/256,5/256).
    """
    leadingcpg_pos = Selected_motif_pos_on_top_strand[0]
    seq_0_lcpg = top_strand[0:leadingcpg_pos]
    count_0_lcpg = seq_0_lcpg.count(umodifiedbase)
    MMtag_value = str(count_0_lcpg)
    for c1, c2 in zip(Selected_motif_pos_on_top_strand[::2], Selected_motif_pos_on_top_strand[1::2]):
        c1_pos_on_read = c1
        c2_pos_on_read = c2
        seq_c1_c2 = top_strand[c1_pos_on_read:c2_pos_on_read]
        count_c1_c2 = seq_c1_c2.count('C') - 1 # c1 was counted in the last round
        MMtag_value = MMtag_value + "," + str(count_c1_c2)
    return MMtag_value
#
#
def Con_posSEQ_to_posRead(read_is_reverse, read_length, pos_on_SEQ_dict):
    """
    Freeze 2024-01-05
    Convert pos_on_SEQ(list) to pos_on_top_strand(5'->3')
    SEQ is from the alignment bam file
    """
    if not read_is_reverse:
        pos_on_top_strand = pos_on_SEQ_dict
        return pos_on_top_strand
    else:
        pos_on_top_strand = {}
        keys = pos_on_SEQ_dict.keys()
        for cpg in keys:
            pos_on_SEQ = pos_on_SEQ_dict[cpg]
            pos_on_top = read_length - pos_on_SEQ - 1
            pos_on_top_strand[cpg] = pos_on_top
    return pos_on_top_strand
#
#
def Buildcglist_dict_Seq(infasta):
    """
    Freeze 2024-08-22:
    Get the 'CG' 1-based position on the reference genome
    """
    cglist_dict = defaultdict(list)
    with pysam.FastxFile(infasta) as fh:
        for entry in fh:
            chr_name = entry.name
            chr_seq = entry.sequence
            position_list = pattern_search(chr_name, chr_seq, pattern='CG')
            # write to file
            cglist_dict[chr_name] = position_list
    # Sort values numerically based on div[1]
    for chromosome in cglist_dict:
        if chromosome in ['chr' + str(i) for i in range(1, 23)]:
            cglist_dict[chromosome] = sorted(cglist_dict[chromosome], key=lambda x: int(x.split("_")[1]))
    return cglist_dict
#

def pattern_search(chr_name, string, pattern):
    """
    Freeze 2024-10-09
    > ATTTCTCTCTCGAAAAAAAACGAAAAAAAA
    > 123456789012345678901234567890
    > print(pattern_search(a,'CG'))
    > [11, 21]
    """
    positions = []
    str_len = len(string)
    pat_len = len(pattern)
    i = 0
    while i <= str_len - pat_len:
        if string[i:i+pat_len] == pattern:
            position = f"{chr_name}_{i+1}"
            positions.append(position)
            i += pat_len  # Skip to the next position after the pattern

        else:
            i += 1  # Move to the next character
    return positions
