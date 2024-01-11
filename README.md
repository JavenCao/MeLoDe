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

Users could train their own *ZMWmeth model* and *Joint-CpG model* to determine read-level 5mC status, and estimate site-level 5mC frequency, respectively.
Actually, not only 5mC modification, but other DNA/RNA modifications fit, as long as the training data can be properly generated

**ZMWmeth model**

To train *ZMWmeth model*, ground-truth Fully-methylated sample and Fully-unmethylated sample are needed.

ZMWmeth model can be trained with respected to HiFi reads, ZMW with subreads, double-strand resolution and single-strand resolution.

**- Step1 - feature extration**

Aligned Bam and unaligned Bam are both supported, as long as fi/fp/ri/rp(HiFi reads), ip/pw(Subreads) tags are present. Please find the commands accordingly. - (1)(2)(3)

(1) From **unaligned HiFi reads**, such as m64284e_230525_074225_hifi.bam you get from the sequencing machine with fi/fp/ri/rp tags:

        ldlr_read --bamfile m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType unaligned \
                  --radius 10 \
                  --outprefix Sample_Name \
                  --outputType Feature

By running this, a file named "Sample_Name.Kmat.txt.gz" will be generated, which stores HiFi read level features.

(2) From **aligned HiFi reads** with fi/fp/ri/rp tags:

        ldlr_read --bamfile aligned.m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType aligned \
                  --radius 10 \
                  --outprefix Sample_Name \
                  --outputType Feature \
                  -cg cpg_motif.pos

By running this, a file named "Sample_Name.Kmat.txt.gz" will be generated, which stores HiFi read level features with respected to each cpg motif on the reference genome.

(3) From **aligned subreads** with ip/pw tags:

        ldlr_read --bamfile aligned.m64284e_230525_074225_subreads.bam \
                  --datatype SUBREAD \
                  --inputType aligned \
                  --radius 10 \
                  --outprefix Sample_Name \
                  --outputType Feature \
                  -cg cpg_motif.pos

By running this, a file named "Sample_Name.Kmat.txt.gz" will be generated, which sotre single molecule level features(in the unit of ZMW) with respected to each cpg motif on the reference genome.

Each line of cpg_motif.pos is a CpG motif coordinate (1-based) on the reference genome, such as:

        chr1	10469
        chr1	10471
        chr1	10484
        chr1	10489
        ...

**- Step2 - model training**

Using feature files above (Fully-methylated sample and Fully-unmethylated sample are both needed), ZMWmeth model can be trained as follows:

        python readtrain.py --MethFeature Me.Kmat.txt.gz \
                        --UnMethFeature UnMe.Kmat.txt.gz \
                        --Tempdir /temp/folder \
                        --OutModel Trained.zmwmeth.model
