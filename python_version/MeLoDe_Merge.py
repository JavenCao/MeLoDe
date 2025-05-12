import argparse
from MeLoDe_Merge import *
#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-Hfile", "--Hfilepre", required=True, help="site-level file (H_model)")
    parser.add_argument("-HCfile", "--HCfilepre", required=True, help="site-level file (HC_model)")
    parser.add_argument("-HCOfile", "--HCOfilepre", required=True, help="site-level file (HCO_model)")
    parser.add_argument("-o", "--outputfile", required=True, help="CpG site level methylation frequency")
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    H_file = args.Hfilepre
    HC_file = args.HCfilepre
    HCO_file = args.HCOfilepre

    out_melode_H = args.outputfile + ".MeLoDe_H.site.txt"
    out_melode_H_HC = args.outputfile + ".MeLoDe_H_HC.txt"
    out_melode_H_HC_HCO = args.outputfile + ".MeLoDe_H_HC_HCO.txt"

    # Initialize dictionaries to store site data

    H_dict = {}
    HC_dict = {}
    HCO_dict = {}

    # Read H_model data

    with open(H_file, "rt") as fpH:
        for line in fpH:
            div = line.rstrip("\n").split("\t")
            pos_key = div[0] + "_" + div[1]
            H_dict[pos_key] = line

    # Read HC_model data, excluding sites already in H_dict

    with open(HC_file, "rt") as fpHC:
        for line in fpHC:
            div = line.rstrip("\n").split("\t")
            pos_key = div[0] + "_" + div[1]
            if pos_key not in H_dict:
                HC_dict[pos_key] = line

    # Read HCO_model data, excluding sites in H_dict or HC_dict

    with open(HCO_file, "rt") as fpHCO:
        for line in fpHCO:
            div = line.rstrip("\n").split("\t")
            pos_key = div[0] + "_" + div[1]
            if pos_key not in H_dict and pos_key not in HC_dict:
                HCO_dict[pos_key] = line

    # Merge dictionaries

    H_HC_dict = {}
    H_HC_dict.update(H_dict)
    H_HC_dict.update(HC_dict)

    H_HC_HCO_dict = {}
    H_HC_HCO_dict.update(H_HC_dict)
    H_HC_HCO_dict.update(HCO_dict)

    # Write output files

    with open(out_melode_H, "wt") as fp_out_melode_H:
        for line in H_dict.values():
            fp_out_melode_H.write(line)

    with open(out_melode_H_HC, "wt") as fp_out_melode_H_HC:
        for line in H_HC_dict.values():
            fp_out_melode_H_HC.write(line)

    with open(out_melode_H_HC_HCO, "wt") as fp_out_melode_H_HC_HCO:
        for line in H_HC_HCO_dict.values():
            fp_out_melode_H_HC_HCO.write(line)
