# LDLR-Meth
Low depth long reads methylation caller

## Quick start:
If you only have HiFi reads(>5x), feel free to use, go ahead!

**Step1: single-molecule 5mC identification**

From **unaligned HiFi reads**:

        ldlr_read --bamfile m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType unaligned \
                  --outprefix Sample_Name \
                  --outputType ModBam

By running this, "Sample_Name.5mc.mod.unaligned.bam" will be generated in which the MM/ML tag is created, tagging the positions and 5mC probabilities of CpGs on reads.
Adding "--KeepKinetics keep" will keep the fi/fp/ri/rp tags, otherwise removed

Next, align the 5mC tagged ModBam to the reference genome using pbmm2:

        pbmm2 align --preset HiFi Ref.fa Sample_Name.5mc.mod.unaligned.bam Aligned_Sample_Name.5mc.mod.unaligned.bam

**Step2:CpG-site level 5mC frequency estimation**

Decode the mod.bam with MM/ML tag into read level prediction:

        ldlr_extract --modbam Aligned_Sample_Name.5mc.mod.unaligned.bam \
                     --outfile Hpre.Read.level.meth.prob.txt.gz
        
        ldlr_site --Hpre Hpre.Read.level.meth.prob.txt.gz \
                  --outfile HiFi.Site.5mC.fre

Output file of "HiFi.Site.5mC.fre" reports DNA 5mC frequency on CpG site.

## Workflow of LDLR-Meth in the low depth design

**Step1:** Synthesis HiFi reads from subreads:

        ccs --hifi-kinetics m64284e_230525_074225.subreads.bam m64284e_230525_074225_hifi.bam

**Step2:** Align HiFi reads and subreads to the reference genome, keeping fi/fp/ri/rp(HiFi reads) and ip/pw(Subreads) tags:

        pbmm2 align --preset HiFi Ref.fa m64284e_230525_074225_hifi.bam Aligned.hifi.bam
       
        pbmm2 align --preset SUBREAD Ref.fa m64284e_230525_074225.subreads.bam Aligned.subreads.bam

**Step3:** Single-molecule level 5mC identification:

For molecules producing HiFi reads:

        ldlr_read --bamfile Aligned.hifi.bam \
                  --datatype HIFI \
                  --inputType aligned \
                  --outprefix Sample_Name.HiFi \
                  -cg cpg_pos.bed \
                  --outputType MoleculeLevel

For all the molecules by subreads (ClosedZMW and OpenZMW):

        ldlr_read --bamfile Aligned.subreads.bam \
                  --datatype SUBREAD \
                  --inputType aligned \
                  --outprefix Sample_Name \
                  -cg cpg_pos.bed \
                  --outputType MoleculeLevel

Each line of "cpg_motif.pos" is a CpG motif coordinate (1-based) on the reference genome, such as:

        chr1	10469
        chr1	10471
        chr1	10484
        chr1	10489
        ...
, which could be created by:

        ldlr_motif --fastafile Ref.fa --pattern 'CG' --outputfile pg_pos.bed

**Step4:** Using molecules with and without HiFi reads for CpG site level DNAm frequency estimation:

        ldlr_site --Hpre Sample_Name.HiFi.SingleMol.pre.txt.gz \
                  --Spre Sample_Name.Subread.SingleMol.pre.txt.gz \
                  --outfile Sample_Name.low_depth.CpG.pre.txt.gz

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
                  --cpgbed cpg_motif.pos

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

**- Step2 - model training**

Using feature files above (Fully-methylated sample and Fully-unmethylated sample are both needed), ZMWmeth model can be trained as follows:

        python readtrain.py --MethFeature Me.Kmat.txt.gz \
                        --UnMethFeature UnMe.Kmat.txt.gz \
                        --OutDir /output/folder/

Specifically, reads/ZMWs/CpGs used for training can be pre-selected. Such as:
* select reads/molecules located at CpG-Island/CpG-desert
* select CpGs at chromsome 10
* use 1 or multiple read(s) for each CpG
* use same number of reads from Fully-methylated sample and Fully-unmethylated sample
* merge multiple samples and/or species
* random select 10M reads
* ... ...






## Full parameter table

