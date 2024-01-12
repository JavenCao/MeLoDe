# LDLR-Meth
Low depth long reads methylation caller

LDLR-Meth (Methylation Caller for Low Depth Long Reads Sequencing) to detect DNA methylation (DNAm) from both high and low depth long reads sequencing, including pb-HiFi reads, pb-ZMW(subreads) and ONT-reads. The pipleline is shown as follows:

<img width="451" alt="image" src="https://github.com/JavenCao/LDLR-Meth/assets/37092316/5c44fc48-e904-4416-9145-a6e08d507e91">

## General Workflow
(1) Before align HiFi reads to the reference genome, profile the CpG methylation status for each read, generating MM/ML tag and dropping fi/fp/ri/rp tags in the output ModBam file
(2) Align mod.Bam to the reference genome.
(3) Estimate CpG site level methylation frequency on the reference genome.

### Quick start: ###
**Step1: single-molecule 5mC identification**
From **unaligned HiFi reads**:

        ldlr_read --bamfile m64284e_230525_074225_hifi.bam \
                  --datatype HIFI \
                  --inputType unaligned \
                  --outprefix Sample_Name \
                  --outputType ModBam

By running this, "Sample_Name.5mc.mod.unaligned.bam" will be generated in which the MM/ML tag wil be created, and fi/fp/ri/rp tags are removed.
Adding "--KeepKinetics keep" will keep the fi/fp/ri/rp tags.

Next, align the above ModBam file to the reference genome using pbmm2:

        pbmm2 align --preset HIFI \
                    --sample "Sample_Name" \
                    "$Ref_GRCh38" Sample_Name.5mc.mod.unaligned.bam Aligned_Sample_Name.5mc.mod.unaligned.bam

By running this, "Aligned_Sample_Name.5mc.mod.unaligned.bam" will be generated.

**Step2:CpG-site level 5mC frequency estimation**

        ldlr_site --bamfile Aligned_Sample_Name.5mc.mod.unaligned.bam \
                  --cpgbed cpg_motif.pos \
                  --outfile Site.5mC.fre

By running this, "Site.5mC.fre" will be generated reporting CpG site level DNA 5mC frequency.

Each line of cpg_motif.pos is a CpG motif coordinate (1-based) on the reference genome, such as:

        chr1	10469
        chr1	10471
        chr1	10484
        chr1	10489
        ...


## Suggested Workflow of LDLR-Meth in the low depth design

Align reads to the reference genome, keeping fi/fp/ri/rp(HiFi reads), ip/pw(Subreads) in the alignment BAM file using pbmm2:
        pbmm2 align --preset CCS --sample "$sampleid" -j ${CPUnum} --log-level "INFO" "$Ref_GRCh38" ${inputfofn} ${outputbam}
        pbmm2 align --preset SUBREAD --sample "$sampleid" -j ${CPUnum} --log-level "INFO" "$Ref_GRCh38" ${inputfofn} ${outputbam}






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

By running this, model file  will be generated.
