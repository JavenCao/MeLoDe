# MeLoDe
 <u>**De**</u>tect DNA **Me**thylaion from **Lo**w-depth long-read sequencing

## Installation

         git clone ### (can we also include Tensorlow_C and a GO into our git, and install them automatically, so that users will not need to do this)
         cd /path/to/MeLoDe-Main
         make
         export PATH=$PATH:/path/to/MeLoDe-Main/bin

## Workflow of MeLoDe in the low-depth long reads sequencing data

MeLoDe is developped to detect DNAm for low-depth LRS data, which is an experimental design used in population-scale study.

Using data from the PacBio platform as an example, below is the workflow to obtain the 5mCpG site-level methylation frequency from such data.

Raw data required: long reads sequencing data with subreads and kinetic signals avaliable.

**Step1:** Synthesis HiFi reads from subreads: (not MeLoDe step)

        ccs --hifi-kinetics subreads.bam hifi.bam

**Step2:** Align HiFi reads and subreads to the reference genome, keeping fi/fp/ri/rp(HiFi reads) and ip/pw(Subreads) tags: (not MeLoDe step)

        pbmm2 align --preset HiFi     Ref.fa hifi.bam      Aligned.hifi.bam
       
        pbmm2 align --preset SUBREAD  Ref.fa subreads.bam  Aligned.subreads.bam

**Step3:** Molecule-level 5mCpG probabilities:

From HiFi reads data:

        MeLoDe_ReadPre --bamfile Aligned.hifi.bam \
                       --InBamStatus aligned \
                       --readformat HIFI \
                       --outputType "MoleculeLevel" \
                       --outprefix "HiFi" \
                       --HDir /path/to/ZMWmeth_HiFi_model \
                       --Reference Ref.fa
                       
"HiFi.SingleMol_5mCpG.txt.gz" will be generated for CpG sites captured by HiFi reads.

From subreads data:

        MeLoDe_ReadPre --bamfile Aligned.subreads.bam \
                       --InBamStatus aligned \
                       --readformat SUBREAD \
                       --outputType "MoleculeLevel" \
                       --outprefix "Subread" \
                       --CDir /path/to/ZMWmeth_ClosedZMW_model \
                       --ODir /path/to/ZMWmeth_OpenZMW_model \                       
                       --Reference Ref.fa

"Subread.SingleMol_5mCpG.txt.gz" will be generated for CpG sites captured by each molecule, using subreads as source data.

**Step4:** CpG site-level DNA methylation frequency prediction:

        MeLoDe_SitePre_Combo --HiFiZMWpre HiFi.SingleMol_5mCpG.txt.gz \
                             --SubreadsZMWpre Subread.SingleMol_5mCpG.txt.gz \
                             --modelDir /path/to/JointCpG_model/Combo \
                             --outputfile Site-level-5mCpG.fre.txt.gz

This step jointly use HiFi-reads and non-HiFi molecules (automatically extracted from SubreadsZMWpre) to calculate the CpG site-level DNA methylation frequency prediction.

The tab-delimited output file is in the following format:

        1. Chr: chromosome
        2. PosC: 1-based coordinate of the CpG on the reference genome
        3. molecule_cov: total read/molecule depths
        4. model_based: model-based prediction
        5. count_based: count-based prediction
        6. Me: methylated reads count
        7. unMeth: unmethylated reads count
        8. cisCpG_ave_weights: surroudning CpGs distance weights


## MeLoDe on HiFi reads
If you only have HiFi reads (say, from Revio platform), MeLoDe is compatible with HiFi reads.

**Step1: read-level 5mC identification**

From **unaligned HiFi reads**:

        MeLoDe_Read --bamfile m64284e_230525_074225_hifi.bam \
                    --InBamStatus unaligned \
                    --readformat HIFI \
                    --outputType ModBam \
                    --outprefix Sample_Name

"Sample_Name.5mc.mod.unaligned.bam" will be generated, in which the MM/ML tag is created, tagging the positions and 5mC probabilities of CpGs on reads.
Adding "--KeepKinetics keep" will keep the fi/fp/ri/rp tags, otherwise removed.

Next, align the 5mC tagged ModBam above to the reference genome using pbmm2:

        pbmm2 align --preset HiFi Ref.fa Sample_Name.5mc.mod.unaligned.bam Aligned_Sample_Name.5mc.mod.unaligned.bam

**Step2: CpG site-level 5mC frequency estimation**

Decode the mod.bam with MM/ML tag into read-level prediction:

        MeLoDe_TagExtract --modbam Aligned_Sample_Name.5mc.mod.unaligned.bam \
                          --outfile Hpre.Read.level.meth.prob.txt.gz
        
        MeLoDe_Site_HiFi --Hpre Hpre.Read.level.meth.prob.txt.gz \
                         --modelDir /path/to/JointCpG_model/HiFi \
                         --outputfile Sample_Name

The tab-delimited output file of "Sample_Name.JointCpG.MeLoDe-HiFi.txt.gz" reports DNA 5mC frequency on CpG sites in the format that is the same as above.

## (Advanced) Train your own models

By default, MeLoDe uses the pre-trained models, which are released together with the software.

        MeLoDe
        └── models
            ├── 1_ZMWmeth_model
            │   ├── double_strand
            │   │   ├── ds.HiFi
            │   │   └── ds.Subreads
            │   └── single_strand
            │       ├── ss.HiFi
            │       └── ss.Subreads
            └── 2_JointCpG_model
                ├── ONT
                │   └── ONT.jointcpg
                └── PB
                    ├── HCO.jointcpg
                    └── H.jointcpg

Users could train their own *ZMWmeth model* and *Joint-CpG model* to determine read-level 5mC status, and estimate site-level 5mC frequency, respectively.
Actually, not only 5mC modification, but other DNA/RNA modifications fit, as long as the training data can be properly generated

### Train ZMWmeth model

To train *ZMWmeth model*, ground-truth Fully-methylated sample and Fully-unmethylated sample are needed.

Say, for 5mC, the training datasets included fully-methylated DNA samples created through treatment with CpG Methyltransferase M.SssI, while fully-unmethylated DNA samples were generated via whole genome amplification.

And, ZMWmeth model can be trained using HiFi reads and subreads in double-strand and single-strand resolution.

**- Step1 - feature extration**

Both aligned Bam and unaligned Bam are supported, as long as fi/fp/ri/rp(HiFi reads), ip/pw(Subreads) tags are present. Please find the commands accordingly. - (1)(2)(3)

(1) From **unaligned HiFi reads**, such as m64284e_230525_074225_hifi.bam you get from the sequencing machine with fi/fp/ri/rp tags:

        MeLoDe ReadExtract --bamfile Aligned.hifi.bam \
                           --datatype HIFI \
                           --inputType unaligned \
                           --radius 10 \
                           --outprefix Sample_ID \
                           --outputType Feature

"Sample_ID.Kmat.txt.gz" will be generated, which stores HiFi read level features, columns are as follows:

        1. Sequence contest for top-strand on its template in 3'-5' direction
        2. IPD features when synthesizing the top-strand
        3. PW features when synthesizing the top-strand
        4. Sequence contest for bot-strand on its template in 3'-5' direction
        5. IPD features when synthesizing the bot-strand
        6. PW features when synthesizing the bot-strand

(2) From **aligned HiFi reads** with fi/fp/ri/rp tags:

        MeLoDe ReadExtract --bamfile Aligned.hifi.bam \
                           --datatype HIFI \
                           --inputType aligned \
                           --radius 10 \
                           --outprefix Sample_ID \
                           --outputType Feature \
                           --ref Ref.fa \
                           --motif "CG"

"Sample_ID.Kmat.txt.gz" will be generated, which stores HiFi read level features with respected to each cpg motif on the reference genome, header is the same as above.

(3) From **aligned subreads** with ip/pw tags:

        MeLoDe ReadExtract --bamfile aligned.Subreads.bam \
                           --datatype SUBREAD \
                           --inputType aligned \
                           --radius 10 \
                           --outprefix Sample_ID \
                           --outputType Feature \
                           --ref Ref.fa \
                           --motif "CG"

"Sample_ID.Kmat.txt.gz" will be generated, which sotre single molecule level features (in the unit of ZMW) with respected to each cpg motif on the reference genome, columns are as follows:

        1. Sequence contest for forward strand of the Reference.fa in 3'-5' direction
        2. IPD features when synthesizing the forward-strand
        3. PW features when synthesizing the forward-strand
        4. Sequence contest for reverse strand of the Reference.fa in 3'-5' direction
        5. IPD features when synthesizing the reverse-strand
        6. PW features when synthesizing the reverse-strand

**- Step2 - model training**

Using feature files above (Fully-methylated sample and Fully-unmethylated sample are both needed), ZMWmeth model can be trained as follows:

        python readtrain.py --MethFeature Me.Kmat.txt.gz \
                            --UnMethFeature UnMe.Kmat.txt.gz \
                            --OutDir /output/folder/

Specifically, training data can be pre-selected for customization. Such as:
* use HiFi reads or subreads features alone to train
* select reads/molecules located at CpG-Island/CpG-desert
* select CpGs at chromsome 10
* use 1 or multiple read(s) for each CpG
* use same number of reads from Fully-methylated sample and Fully-unmethylated sample
* merge multiple samples and/or species
* random select 10M reads
* ... ...

### Train Joint-CpG model

To train *Joint-CpG model*, wild-type samples with ground-truth site-level methylation frequency(such as determined by WGBS) are needed.

Joint-CpG model can be trained with respected to HiFi reads, combination of HiFi-reads + nonHiFi-subreads.

**- Step1 - Single-molecule level 5mC identification**

For molecules producing HiFi reads:

        ldlr_read --bamfile Wildtype.Aligned.hifi.bam \
                  --datatype HIFI \
                  --inputType aligned \
                  --outprefix Sample_Name.HiFi \
                  -cg cpg_pos.bed \
                  --outputType MoleculeLevel

**- Step2 - model training**

        python sitetrain.py --hifizmw Sample_Name.HiFi.txt.gz \
                            --TrueSiteLevel Sample_Name.wgbs.cpg.bed \
                            --OutDir /output/folder/





## Well trained models





## Full parameter table

