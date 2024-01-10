# LDLR-Meth
Low depth long reads methylation caller

LDLR-Meth (Methylation Caller for Low Depth Long Reads Sequencing) to detect DNA methylation (DNAm) from both high and low depth long reads sequencing, including pb-HiFi reads, pb-ZMW(subreads) and ONT-reads. The pipleline is shown as follows:

<img width="451" alt="image" src="https://github.com/JavenCao/LDLR-Meth/assets/37092316/5c44fc48-e904-4416-9145-a6e08d507e91">

## step1: single-molecule 5mC identification
### HiFi reads can be treated as a unit
### non_HiFi molecules can be called by subreads
### ONT reads is called by Dorado[src=https://github.com/nanoporetech/dorado]



## step2: CpG site level 5mC frequency estimation

## step3: Use HiFi and non_HiFi molecules in the low-depth design 



## Train your own models

Users could train their own ZMWmeth model and Joint-CpG model for single-molecule 5mC calling and site level 5mC frequency estimation.

* ZMWmeth model

(step1 - feature extration) From unaligned HiFi reads, such as m64284e_230525_074225_hifi.bam you got from the sequencer with fi/fp/ri/rp tags:

        ldlr_read --bamfile m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType unaligned \
                  --outprefix Sample_Name \
                  --outputType Feature

    By running this, a file named "Sample_Name.Kmat.txt.gz" will be generated, which store read level features, from which the model can be trained.

    OR, from aligned HiFi reads with fi/fp/ri/rp tags:

        ldlr_read --bamfile aligned.m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType aligned \
                  --outprefix Sample_Name \
                  --outputType Feature \
                  -cg cpg_motif.pos

    By running this, a file named "Sample_Name.Kmat.txt.gz" will be generated, which sotre read level features with respected to each cpg motif on the reference genome.



