# MeLoDe
 <u>**De**</u>tect DNA **Me**thylaion from **Lo**w-depth long-read sequencing

## Installation

         git clone ### (can we also include Tensorlow_C and a GO into our git, and install them automatically, so that users will not need to do this)
         cd /path/to/MeLoDe-Main
         make
         export PATH=$PATH:/path/to/MeLoDe-Main/bin

## Workflow of MeLoDe in the low-depth long reads sequencing data

MeLoDe is developped to detect DNAm for low-depth long reads sequencing data, which is an experimental design commonly used in population-scale study.

Using data from the PacBio platform as an example, below is the workflow to obtain the 5mCpG site-level methylation frequency from such data.

Raw data required: long reads sequencing data with subreads and kinetic signals avaliable (MeLoDe is versatile, and can be applied to HiFi reads directly, in case you don't have access to the subreads data, say data from the Revio platform).

**Step1:** Synthesis HiFi reads from subreads: (not MeLoDe step)

        ccs --hifi-kinetics subreads.bam hifi.bam

**Step2:** Align HiFi reads and subreads to the reference genome, keeping fi/fp/ri/rp(HiFi reads) and ip/pw(Subreads) tags: (not MeLoDe step)

        pbmm2 align --preset HiFi     Ref.fa hifi.bam      Aligned.hifi.bam
       
        pbmm2 align --preset SUBREAD  Ref.fa subreads.bam  Aligned.subreads.bam

**Step3:** Molecule-level 5mCpG probabilities:

Seprately, for HiFi reads data and for subreads, run the following command lines:

For HiFi reads data:

        MeLoDe_Read -bamfile Aligned.hifi.bam \
                    -InBamStatus aligned \
                    -readformat HIFI \
                    -outputType "MoleculeLevel" \
                    -outprefix "HiFi" \
                    -HDir /path/to/ZMWmeth_HiFi_model \
                    -Reference Ref.fa \
                    -Scale \
                    -Processor 10

"HiFi.SingleMol_5mCpG.txt.gz" will be generated, in which CpG sites captured by HiFi reads are predicted.

From subreads data:

        MeLoDe_Read -bamfile Aligned.subreads.bam \
                    -InBamStatus aligned \
                    -readformat SUBREAD \
                    -outputType "MoleculeLevel" \
                    -outprefix "Subread" \
                    -CDir /path/to/ZMWmeth_ClosedZMW_model \
                    -ODir /path/to/ZMWmeth_OpenZMW_model \
                    -Reference Ref.fa \
                    -Scale \
                    -Processor 10

"Subread.SingleMol_5mCpG.txt.gz" will be generated, in which CpG sites captured by each molecule are predicted, using subreads as source data.

These molecule-level output file is tab-delimited in the following format:

        1. Chr: chromosome,
        2. Pos1: 1-based coordinate of the CpG on the reference genome
        3. ZMWid: molecule id
        4. ZMWtype: zmwType, "H": HiFi reads, "C": Closed-ZMW, "O": Open-ZMW
        5. subdepth: subreads depth for the molecule
        6. Hap: haplotype, from HP tag
        7. HapBlk: haplotype black, fro PS tag
        8. strand: strand
        9. y_pred: 5mCpG probability,
        10. y_pred_bi: methylation status, y_pred>0.5 -> 1; y_pred<0.5 -> 0

**Step4:** CpG site-level DNA methylation frequency prediction:

        MeLoDe_Combo -HiFiZMWpre HiFi.SingleMol_5mCpG.txt.gz \
                     -SubreadsZMWpre Subread.SingleMol_5mCpG.txt.gz \
                     -modelDir /path/to/JointCpG_model/Combo \
                     -outputfile Site-level-5mCpG.fre.txt.gz

This step jointly use HiFi-molecules and non-HiFi molecules (automatically extracted from SubreadsZMWpre) to predict the CpG site-level DNA methylation frequency.

The tab-delimited output file is in the following format:

        1. Chr: chromosome
        2. PosC: 1-based coordinate of the CpG on the reference genome
        3. molecule_cov: total read/molecule depths
        4. model_based: model-based prediction
        5. count_based: count-based prediction ( = methylated reads count / total reads count)
        6. Me: methylated reads count
        7. unMeth: unmethylated reads count
        8. cisCpG_ave_weights: surroudning CpGs distance weights

## MeLoDe on HiFi reads

If you only have HiFi reads (say, from Revio platform), MeLoDe is compatible with HiFi reads. Please follow the steps below:

**Step1: read-level 5mCpG prediction**

From **unaligned HiFi reads**:

        MeLoDe_Read -bamfile m64284e_230525_074225_hifi.bam \
                    -InBamStatus unaligned \
                    -readformat HIFI \
                    -outputType "ModBam" \
                    -outprefix Sample_Name \
                    -HDir /path/to/ZMWmeth_HiFi_model \
                    -Scale \
                    -Processor 10

"Sample_Name.5mc.mod.unaligned.bam" will be generated, in which the MM/ML tags are created, tagging the positions (MM) and 5mCpG probabilities (ML) of CpGs on each HiFi reads.
Adding "--KeepKinetics keep" will keep the fi/fp/ri/rp tags, otherwise removed.

**Step2: Alignment**

align the 5mC tagged ModBam above to the reference genome using pbmm2:

        pbmm2 align --preset HiFi Ref.fa Sample_Name.5mc.mod.unaligned.bam Aligned_Sample_Name.5mc.mod.unaligned.bam

**Step3: CpG site-level 5mC frequency estimation**

Decode the mod.bam with MM/ML tag into read-level prediction:

        MeLoDe_TagExtract -bamfile Aligned_Sample_Name.5mc.mod.unaligned.bam \
                          -outfile Hpre.Read.level.meth.prob.txt.gz
        
        MeLoDe_HiFi -Hpre Hpre.Read.level.meth.prob.txt.gz \
                    -modelDir /path/to/JointCpG_model/HiFi \
                    -outputfile Sample_Name

The tab-delimited output file of "Sample_Name.JointCpG.MeLoDe-HiFi.txt.gz" reports DNA 5mCpG methylation frequency on CpG sites in the format that is the same as above.

## (Advanced) Train your own models

By default, MeLoDe uses the pre-trained models, which are released together with the software.

        ├──────────────────────────────────────────────────────
        ├── Saved_model
        │   ├── JointCpG-Combo_model
        │   │   ├── saved_model.pb
        │   │   └── variables
        │   │       ├── variables.data-00000-of-00001
        │   │       └── variables.index
        │   ├── JointCpG-HiFi_model
        │   │   ├── saved_model.pb
        │   │   └── variables
        │   │       ├── variables.data-00000-of-00001
        │   │       └── variables.index
        │   ├── ZMWmeth_HiFi_model
        │   │   └── ds
        │   │       ├── saved_model.pb
        │   │       └── variables
        │   │           ├── variables.data-00000-of-00001
        │   │           └── variables.index
        │   └── ZMWmeth_Subreads_model
        │       ├── Closed-ZMW
        │       │   ├── saved_model.pb
        │       │   └── variables
        │       │       ├── variables.data-00000-of-00001
        │       │       └── variables.index
        │       └── Open-ZMW
        │           ├── saved_model.pb
        │           └── variables
        │               ├── variables.data-00000-of-00001
        │               └── variables.index
        ├──────────────────────────────────────────────────────


Users could train their own *ZMWmeth model* and *Joint-CpG model* to determine read-level 5mCpG status, and estimate site-level 5mCpG frequency, respectively.
Actually, not only 5mC modification, but other DNA/RNA modifications fit, as long as the training data can be properly generated

### Train ZMWmeth model

To train *ZMWmeth model*, ground-truth Fully-methylated sample and Fully-unmethylated sample are needed.

Say, for 5mCpG, the training datasets included fully-methylated DNA samples created through treatment with CpG Methyltransferase M.SssI, while fully-unmethylated DNA samples were generated via whole genome amplification.

And, ZMWmeth model can be trained using HiFi reads and subreads in double-strand and single-strand resolution.

**- Step1 - feature extration**

Both aligned Bam and unaligned Bam are supported, as long as fi/fp/ri/rp(HiFi reads), ip/pw(Subreads) tags are present. Please find the commands accordingly. - (1)(2)(3)

(1) From **unaligned HiFi reads**, such as m64284e_230525_074225_hifi.bam you get from the sequencing machine with fi/fp/ri/rp tags:

        MeLoDe_Read -bamfile m64284e_230525_074225_hifi.bam \
                    -InBamStatus unaligned \
                    -readformat HIFI \
                    -outputType "Feature" \
                    -outprefix Sample_Name \
                    -Scale \
                    -Processor 10

"Sample_Name.Kmat.txt.gz" will be generated, which stores HiFi read level features (tab-delimited), columns are as follows:

        1. "TOP": flag, suggeting it is from top_unaligned file
        2. zmwname: ZMWname
        3. pos_on_read: CpG position of the molecule
        4. fn: subreads depth from forward strand, from fn tag
        5. rn: subreads depth from reverse strand, from rn tag
        6. Sequence contest for top-strand on its template in 3'-5' direction
        7. IPD features when synthesizing the top-strand, "," split
        8. PW features when synthesizing the top-strand, "," split
        9. Sequence contest for bot-strand on its template in 3'-5' direction
        10. IPD features when synthesizing the bot-strand, "," split
        11. PW features when synthesizing the bot-strand, "," split

(2) From **aligned HiFi reads** with fi/fp/ri/rp tags:

        MeLoDe_Read -bamfile aln.m64284e_230525_074225_hifi.bam \
                    -InBamStatus aligned \
                    -readformat HIFI \
                    -outputType "Feature" \
                    -outprefix Sample_Name \
                    -Scale \
                    -Processor 10 \
                    -Reference GRCh38.fa \
                    -Scale \
                    -Processor 10

(3) From **aligned subreads** with ip/pw tags:

        MeLoDe_Read -bamfile aln.subreads.m64284e_230525_074225_hifi.bam \
                    -InBamStatus aligned \
                    -readformat SUBREAD \
                    -outputType "Feature" \
                    -outprefix Sample_Name \
                    -Scale \
                    -Processor 10 \
                    -Reference GRCh38.fa \
                    -Scale \
                    -Processor 10

**- Step2 - model training**

Using feature files above (Fully-methylated sample and Fully-unmethylated sample are both needed), ZMWmeth model can be trained as follows:

        python readtrain.py --MethFeature Me.Kmat.txt.gz \
                            --UnMethFeature UnMe.Kmat.txt.gz \
                            --OutDir /output/folder/

Specifically, training data can be pre-selected for customization. Such as:
* use HiFi reads or subreads features exclusively to train
* select reads/molecules located at CpG-Island/CpG-desert
* select CpGs at chromsome 10
* use 1 or multiple reads for each CpG
* use same number of reads from Fully-methylated sample and Fully-unmethylated sample
* merge multiple samples and/or species
* random select 10M reads
* ... ...

### Train Joint-CpG model

To train *Joint-CpG model*, wild-type samples with ground-truth site-level methylation frequency (such as determined by WGBS) are needed.

Joint-CpG model can be trained with respected to HiFi reads, combination of HiFi-reads + nonHiFi-subreads.

Below, we showcase how the Joint-CpG model in MeLoDe-Combo was trained, which utilized both HiFi molecules in HiFi reads, and non-HiFi molecules in subreads.

**- Step1 - Single-molecule level 5mCpG identification**

For HiFi reads data:

        MeLoDe_Read -bamfile Aligned.hifi.bam \
                    -InBamStatus aligned \
                    -readformat HIFI \
                    -outputType "MoleculeLevel" \
                    -outprefix "HiFi" \
                    -HDir /path/to/ZMWmeth_HiFi_model \
                    -Reference Ref.fa \
                    -Scale \
                    -Processor 10

"HiFi.SingleMol_5mCpG.txt.gz" will be generated, in which CpG sites captured by HiFi reads are predicted.

From subreads data:

        MeLoDe_Read -bamfile Aligned.subreads.bam \
                    -InBamStatus aligned \
                    -readformat SUBREAD \
                    -outputType "MoleculeLevel" \
                    -outprefix "Subread" \
                    -CDir /path/to/ZMWmeth_ClosedZMW_model \
                    -ODir /path/to/ZMWmeth_OpenZMW_model \
                    -Reference Ref.fa \
                    -Scale \
                    -Processor 10

"Subread.SingleMol_5mCpG.txt.gz" will be generated, in which CpG sites captured by each molecule are predicted, using subreads as source data.

**- Step2 - create CpG study-window for Joint-CpG model**

        python MeLoDe_Format_Joint_CpG.py -hifizmw HiFi.SingleMol_5mCpG.txt.gz \
                                          -subaredzmw Subread.SingleMol_5mCpG.txt.gz \
                                          -trueMe WGBS.cpg.txt \
                                          -maxCpGnum 3000 \
                                          -o Sample_CpG

**- Step3 - run training for Joint-CpG model**

Put well-formatted CpG study-windows into a .fofn (file of file input), such as:

        ls  *npz >> JointCpG.npz.fofn

        python sitetrain -npzl JointCpG.npz.fofn \
                        -r 5 \
                        -bz 5 \
                        -c ${cpunum} \
                        -o out.MeLoDe.Combo.JointCpG.model





## Well trained models





## Full parameter table

