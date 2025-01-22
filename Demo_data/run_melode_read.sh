# # # # # # # # # # # # Please replace accordingly # # # # # # # # # # #
rundir=/path/to/MeLoDe-main
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
inbamfile=${rundir}/Demo_data/Demo.unaligned.meHG005_hifi.bam
outprefix=${rundir}/Demo_data/5mCtaged.Demo.unaligned.meHG005_hifi
MeLoDe_Read -bamfile ${inbamfile} \
            -InBamStatus unaligned \
            -readformat HIFI \
            -outputType "ModBam" \
            -outprefix ${outprefix} \
            -HDir /storage/yangjianLab/caoyujie/software/test/MeLoDe-main/Saved_model/ZMWmeth_HiFi_model/ds \
            -Scale \
            -Processor 1