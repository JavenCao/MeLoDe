from motif import *
import pysam
import argparse

if __name__ == "__main__":
    # --------
    args = get_args()
    infasta = args.fastafile
    patternstring = args.pattern
    outputfile = args.outputfile
    # --------
    fout = open(outputfile, "a+")
    with pysam.FastxFile(infasta) as fh:
        for entry in fh:
            chr_name = entry.name
            chr_seq = entry.sequence
            position_list = pattern_search(chr_seq, patternstring)
            # write to file
            index = list(position_list)
            for x in index:
                newline = str(chr_name) + "\t" + str(x) + "\n"
                fout.write(newline)
else:
    pass

# ----------------------------------------------------------
def pattern_search(string, pattern):
    """
    > ATTTCTCTCTCGAAAAAAAACGAAAAAAAA
    > 123456789012345678901234567890
    > print(pattern_search(a,'CG'))
    > [11, 21]
    """
    position = []
    str_len = len(string)
    pat_len = len(pattern)
    for i in range(0,str_len-pat_len+1):
        for j in range(0,pat_len):
            if string[i+j]!=pattern[j]:
                break
        else:
            position.append(i+1)
    return position
# ----------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta","--fastafile")
    parser.add_argument("-p","--pattern",help="sub-string to find in the fasta files, say 'CG', 'GCTA' ")
    parser.add_argument("-o","--outputfile")
    return parser.parse_args()
# ----------------------------------------------------------
