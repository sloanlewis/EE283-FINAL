# EE283-FINAL

## Background information on the project
Using a nonhuman primate model of alcohol self-administration, our lab is interested in understanding the impact of chronic alcohol consuption on circulating immune cells. For this project, we are specifically interested in discovering possible biomarkers of alcohol consumption by analysing the miRNA contents of extracellular vesicles in blood plasma. The macaques have open access to a 4%(w/v) ethanol solution for ~12 months and establish a stable drinking phenotype (grams EtOH consumed per day) over the year. 

## Experimental design
1) Extract extracellular vesicles from plasma using a membrane-based affinity kit
2) Lyse vesicles and recover total RNA
3) Build cDNA libraries of the miRNAs using the Qiaseq miRNA libaray kit
4) Multiplex and sequence the libraries on a HiSeq 2500
5) Assess the quality of the sequences (fastqc) after de-multiplexing
6) Trim off adapter sequences and set quality and length cut-offs (trim_galore)
7) Align sequences to the _Macaca mulatta_ genome and generate alignment report 
8) Examine differential miRNA expression using edgeR

## Notes
In the Messaoudi Lab, we are currently still using the UC Riverside HPC as the lab moved from there recently. My code is all written to run on that cluster as that is where I will need to run it in the future.

## Step 1: Generate fastqc reports for all sequencing files

```
fastqc your_file.fq.gz

```
look at these reports to see starting quality of the sequencing run and to see which adapters are being detected

## Step 2: Trim all of the files in a directory according to trim_galore miRNA parameters

```
sh ../../../small-RNA/trim_galore_smallRNA_dir.sh your_directory

```
_Always set a --length 18 minumum or trim_galore will automatically trim sequences less that 20bp_

_Make sure to check the file paths_

