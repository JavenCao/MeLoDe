cd MeLoDe-main/Demo_data
export PATH=$PATH:/storage/yangjianLab/caoyujie/software/test/MeLoDe-main/bin
export LD_LIBRARY_PATH=/storage/yangjianLab/caoyujie/software/test/MeLoDe-main/libtensorflow/lib/:$LD_LIBRARY_PATH

# # # # # # # # # # # # Please replace accordingly # # # # # # # # # # #
rundir=/path/to/MeLoDe-main
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


inbamfile=${rundir}/Demo_data/Demo.meHG005_hifi.bam
outprefix=${rundir}/Demo_data/5mCtaged.Demo.meHG005_hifi



    MeLoDe_Read -bamfile ${inbamfile} \
                -InBamStatus unaligned \
                -readformat HIFI \
                -outputType "ModBam" \
                -outprefix ${outprefix} \
                -HDir /storage/yangjianLab/caoyujie/software/test/MeLoDe-main/Saved_model/ZMWmeth_HiFi_model/ds \
                -Scale \
                -Processor 1


inbamfile=${rundir}/Demo_data/Demo.wgaHG005_hifi.bam
outprefix=${rundir}/Demo_data/5mCtaged.Demo.wgaHG005_hifi
    MeLoDe_Read -bamfile ${inbamfile} \
                -InBamStatus unaligned \
                -readformat HIFI \
                -outputType "ModBam" \
                -outprefix ${outprefix} \
                -HDir /storage/yangjianLab/caoyujie/software/test/MeLoDe-main/Saved_model/ZMWmeth_HiFi_model/ds \
                -Scale \
                -Processor 1
