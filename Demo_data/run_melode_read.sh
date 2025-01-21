cd MeLoDe-main/Demo_data
export PATH=$PATH:/storage/yangjianLab/caoyujie/software/test/MeLoDe-main/bin
export LD_LIBRARY_PATH=/storage/yangjianLab/caoyujie/software/test/MeLoDe-main/libtensorflow/lib/:$LD_LIBRARY_PATH

# # # # # # # # # # # # Please replace accordingly # # # # # # # # # # #
rundir=/path/to/MeLoDe-main
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


inbamfile=${rundir}/Demo_data/Demo.MssI.top.hifi.bam
outprefix=${rundir}/Demo_data/5mCtaged.Demo.MssI.top.hifi



    MeLoDe_Read -bamfile ${inbamfile} \
                -InBamStatus unaligned \
                -readformat HIFI \
                -outputType "ModBam" \
                -outprefix ${outprefix} \
                -HDir /storage/yangjianLab/caoyujie/software/test/MeLoDe-main/Saved_model/ZMWmeth_HiFi_model/ds \
                -Scale \
                -Processor 1


inbamfile=${rundir}/Demo_data/Demo.WGA.top.hifi.bam
outprefix=${rundir}/Demo_data/5mCtaged.Demo.WGA.top.hifi
    MeLoDe_Read -bamfile ${inbamfile} \
                -InBamStatus unaligned \
                -readformat HIFI \
                -outputType "ModBam" \
                -outprefix ${outprefix} \
                -HDir /storage/yangjianLab/caoyujie/software/test/MeLoDe-main/Saved_model/ZMWmeth_HiFi_model/ds \
                -Scale \
                -Processor 1
