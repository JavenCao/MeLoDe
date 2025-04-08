#!/bin/bash
#SBATCH --job-name=bamK_test
#SBATCH -c 5
#SBATCH --mem=50g
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -o bamK_test.out
#SBATCH -e bamK_test.err
#
#
cpunum=1
#source ~/.bashrc
#source ~/.bash_profile
#
wkd=/home/yangjianLab/yangwen/goworking/tag_decode/test
#inBam=${wkd}/chr22.golang.hifi.aligned.modification.sorted.bam
#inBam="/storage/yangjianLab/caoyujie/project/meth/5_Replication/8_CKCG_high_2Dmodel/melode/1_modbam/ds_aligned_modbam/CKCG-00309/ds_3x/ds_3x.melode.CKCG-00309.hifi.bam"
inBam="/storage/yangjianLab/caoyujie/project/CRC/2_Bam_file/hifi_by_melode/CRC0001N/phased.sort.CRC0001N.aln.bam"

cmd="/home/yangjianLab/yangwen/goworking/tag_decode/bin/MeLoDe_TagExtract"
#
start=$(date +%s.%N)

/usr/bin/time ${cmd} -bamfile ${inBam} \
             -outfile ${wkd}/CRC0001N.txt \
             -topN 1000
#
#  -bpq 0 \

dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time: %.1f seconds\n" $dur
