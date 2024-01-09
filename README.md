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
