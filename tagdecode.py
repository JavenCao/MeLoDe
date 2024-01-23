import pysam
import gzip
import argparse
import numpy as np
from tagdecode import *
#
#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-bam", "--bamfile", required=True, help="Aligned HiFi Bam file with MM/ML tags")
    parser.add_argument("-o", "--outfile", required=True, help="Read levle prediction file")
    return parser.parse_args()                        
#
#
if __name__ == '__main__':
    args = get_args()
    bamf = args.bamfile
    outfile = args.outfile
    #
    Decode_tag_to_txt(bamf, outfile)
else:
    pass


def Decode_tag_to_txt(bamf, outfile):
    """
    2024-01-22
    """
    samfile = pysam.AlignmentFile(bamf, "rb")
    fout = gzip.open(outfile, "wb")
    for read in samfile.fetch():
        if (not read.is_secondary) and (not read.is_supplementary):
            if read.has_tag("MM") and read.has_tag("ML"):
                mmtag = np.array((read.get_tag('MM', with_value_type=False)).rstrip(";").split(",")[1:],dtype=int)
                mltag = np.array(read.get_tag('ML', with_value_type=False),dtype=int)
                read_SEQ_string = read.query_sequence
                read_is_reverse = read.is_reverse
                read_cigartuples = read.cigartuples
                aln_ref_chr = read.reference_name
                aln_ref_start = read.reference_start #0-based leftmost coordinate
                zmwname = str(read.query_name.split("/")[0]) + "_" + str(read.query_name.split("/")[1])
                ZMWtype = 'H'
                fn = read.get_tag("fn")
                rn = read.get_tag("rn")
                total_subreads_dep = int(fn) + int(rn)
                haplotype, haplotype_block = (str(read.get_tag("HP")), str(read.get_tag("PS"))) if read.has_tag("HP") else ("X", "X")# HP:i:2   PS:i:11457
                strand_sign = "+"
                #
                if read_is_reverse:
                    top_strand = str_complementary(read_SEQ_string[::-1])
                else:
                    top_strand = read_SEQ_string
                #
                #
                tag_on_SEQ_pos, tag_on_SEQ_prob = mm_ml_tag_on_SEQ(read_SEQ_string, read_is_reverse, mmtag, mltag)
                # tag_on_Ref_pos may have different length as tag_on_SEQ_pos
                # so tag_on_SEQ_prob should be updated to exclude these postions
                tag_on_Ref_pos, tag_on_SEQ_prob_update = Con_posSEQ_to_posRef(tag_on_SEQ_pos, aln_ref_chr, aln_ref_start, read_cigartuples, tag_on_SEQ_prob)

                # # write out # #
                for i in range(len(tag_on_Ref_pos)):
                    predictions_bi = (tag_on_SEQ_prob_update[i] >= 0.5).astype(int)
                    outline =f'{aln_ref_chr}\t{tag_on_Ref_pos[i]}\t{zmwname}\t{ZMWtype}\t{total_subreads_dep}\t{haplotype}\t{haplotype_block}\t{strand_sign}\t{tag_on_SEQ_prob_update[i]}\t{predictions_bi}\n'
                    fout.write(outline.encode('utf-8'))                    
#
#
def Con_posSEQ_to_posRef(tag_on_SEQ_pos, aln_ref_chr, aln_ref_start, read_cigartuples, tag_on_SEQ_prob):
    """
    Freeze,2024-01-22
    """
    # same length as the input tag_on_SEQ_pos
    tag_on_Ref_pos = ['NA'] * len(tag_on_SEQ_pos)
    cigar_blk_num = len(read_cigartuples)
    
    for index in range(0, len(tag_on_SEQ_pos)):
        SEQ_pos = tag_on_SEQ_pos[index]
        read_pos_blk_start = 0
        read_pos_blk_end = 0
        Ref_walking_pos_start = aln_ref_start
        Ref_walking_pos_end = 0
        flag = 1
        #
        for cigar_index in range(cigar_blk_num):
            op, count = read_cigartuples[cigar_index]
            #
            # pos on SEQ
            if (op == 0 or op == 1 or op == 4 or op == 7 or op == 8):
                read_pos_blk_end = read_pos_blk_end + count
                read_pos_blk_start = read_pos_blk_end - count
            #
            # pos on Ref
            if (op == 0 or op == 2 or op == 3 or op == 7 or op == 8):
                if (Ref_walking_pos_start == aln_ref_start) and (flag == 1):
                    Ref_walking_pos_end = Ref_walking_pos_start + count
                    flag = 2
                else:       
                    Ref_walking_pos_start = Ref_walking_pos_end + 1
                    Ref_walking_pos_end = Ref_walking_pos_end + count
            # - - - - -
            # capture the cpg on read on the moving read
            if (read_pos_blk_start <= SEQ_pos) and (read_pos_blk_end >= SEQ_pos):
                if (op == 0 or op == 1 or op == 4 or op == 7 or op == 8):
                    last_op_needed = SEQ_pos - read_pos_blk_start
                    if (op == 0 or op == 2 or op == 3 or op == 7 or op == 8):
                        ref_pos = Ref_walking_pos_start + last_op_needed
                        tag_on_Ref_pos[index] = ref_pos
                        break
    # finally, exlcude those with 'NA' and trimed from tag_on_SEQ_prob
    indices = []
    indices = [i for i in range(len(tag_on_Ref_pos)) if tag_on_Ref_pos[i] == 'NA']
    tag_on_Ref_pos_update = [tag_on_Ref_pos[i] for i in range(len(tag_on_Ref_pos)) if i not in indices]
    tag_on_SEQ_prob_update = [tag_on_SEQ_prob[i] for i in range(len(tag_on_Ref_pos)) if i not in indices]
    return tag_on_Ref_pos_update, tag_on_SEQ_prob_update
#
#
def mm_ml_tag_on_SEQ(read_SEQ_string, read_is_reverse, mmtag, mltag):
    """
    2024-01-22
            Read  AGTCTAGACTCCGTAATTACTCGCCTAG...
            C_index  1    2 34       5 6 78
   (0-based)C_index  0    1 23       4 5 67
            C_pos    3    8 01       9 1 34         Cpos_on_top_strand=[3, 8, 10, 11, 19, 21, 23, 24]#0-based
            CpG              *         *
            1-based:         4         6
            0-based:         3         5
        MM:Z:C+m,3,1,...   # CpG sites are at C #4 (1+3) and #6 (1+3+1+1)
        ML:B:C,249,4,...   # probability of methylation at the first CpG is in [249/256,250/256); second CpG is in [4/256,5/256).
    """
    total_targets = len(mmtag)
    read_query_length = len(read_SEQ_string)
    tag_on_TOP_pos = []
    tag_on_TOP_prob = []

    tag_on_SEQ_pos = []
    tag_on_SEQ_prob = []
    # ---
    if read_is_reverse:
        top_strand = str_complementary(read_SEQ_string[::-1])
    else:
        top_strand = read_SEQ_string
    # ---
    # 2024-01-22
    Cpos_on_top_strand = [i for i in range(len(top_strand)) if top_strand[i] == 'C']
    for index in range(total_targets):
        # 0-based index on the python vector
        C_cpg_index = np.sum(mmtag[0:index+1]) + index#[3,5,...]
        cpg_pos = Cpos_on_top_strand[C_cpg_index]#[11,21,...]
        tag_on_TOP_pos.append(cpg_pos)
        tag_on_TOP_prob.append(mltag[index]/256)
    # ---
    if read_is_reverse:
        tag_on_SEQ_pos = [read_query_length-x-2 for x in (tag_on_TOP_pos[::-1])]
        tag_on_SEQ_prob = tag_on_TOP_prob[::-1]
    else:
        tag_on_SEQ_pos = tag_on_TOP_pos
        tag_on_SEQ_prob = tag_on_TOP_prob
    # ---
    return tag_on_SEQ_pos,tag_on_SEQ_prob
#
#
def str_complementary(string_a):
    """Freeze,2024-01-04"""
    complementary_string_a = string_a.translate(string_a.maketrans("ATCG", "TAGC"))
    return complementary_string_a
