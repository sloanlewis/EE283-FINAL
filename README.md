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
7) Align sequences to the _Macaca mulatta_ genome and generate alignment report (bowtie2)
8) Examine differential miRNA expression using edgeR

## Notes
In the Messaoudi Lab, we are currently still using the UC Riverside HPC as the lab moved from there recently. My code is all written to run on that cluster as that is where I will need to run it in the future. This pipeline is adapted from the systemPipeR pipeline availible here: https://bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html

## Step 1: Generate fastqc reports for all sequencing files

```
fastqc your_file.fq.gz
```
look at these reports to see starting quality of the sequencing run and to see which adapters are being detected

## Step 2: Trim all of the files in a directory according to trim_galore miRNA parameters

```
sh trim_galore_smallRNA_dir.sh your_directory
```
_Always set a --length 18 minumum or trim_galore will automatically trim sequences less that 20bp_

_Make sure to check the input file names and the paths to the scripts are correct_

Check the fastqc output files to make sure all trimming parameters were correct and record the number of sequences post-trimming for your records. Check for any leftover adapter sequences or over-represented nucleotides at the 5' or 3' end that may need to be clipped off.

*Note for the sequnces from 03/05/19, 1bp was clipped off the 5' end due to over-representation of T at the 5' position in almost every sequence. The script small_RNA_trimming_clip.sh will do that.*

## Step 3: Generate the targets file, data, and results directory as described in systemPipeR

## Step 4: Align trimmed files to the _Macaca mulatta_ genome
1) Download the FASTA and GTF files from ensembl
2) Build a bowtie2 index for the reference genome
3) Create a miRNA GTF file from the whole GTF:

```
grep 'miRNA' Macaca_mulatta.Mmul_8.0.1.95.gtf > Macaca_mulatta.Mmul_8.0.1.95_miRNA.gtf
```
Copy the header from the whole GTF to your miRNA GTF

Create the sqlite file you will need for counting:

```
R
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file="data/Macaca_mulatta.Mmul_8.0.1.95_miRNA.gtf", format="gtf", dataSource="ENSEMBL", organism="Macaca mulatta")
saveDb(txdb, file="./data/Macaca_mulatta.Mmul_8.0.1.95_miRNA.sqlite")
```

*03/15/19 Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa*

*03/15/19 Macaca_mulatta.Mmul_8.0.1.95.gtf*

4) Make sure small_RNA_alignemnt.R, bowtie-smallRNA.param, .BatchJobs.R, and slurm.tmpl are downloaded for this step

5) Run the script small_RNA_alignment.R in R to submit the alignments to the cluster
  
 _Run the script from your main directory_
 
 _slurm.tmpl must be in the directory you are running the script from or it will not work_
 
 6) Check alignment stats
 
 ```
R
library(systemPipeR)
library(GenomicFeatures)
#Read the targets file
targets <- read.delim("targets.txt", comment.char = "#")
#Check if targets file is correct
targets
#Create args
args <- systemArgs(sysma="bowtie-smallRNA.param", mytargets="targets.txt")
moduleload(modules(args))
#Check if alignment files exist
file.exists(outpaths(args))
read_statsDF <- alignStats(args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```
 
 ## Generate counts and RPKM files for DEG analysis
 
 Submit small_RNA_readcounting_cluster.R to the cluster using counts_batch.sh
 
 ```
 sbatch -p highmem --mem=100g --time=24:00:00 counts_batch.sh
 ```

_Check all file names and paths before running._

## Run small_RNA_edgeR.R line by line to generate differential gene expression analysis file based on comparisons specified in targets file

## Generate some preliminary graphs for the data

```
#Sample outlier analysis
library(ape)
rpkmDFmiR <- read.delim("./results/rpkmDFmiR.xls", row.names=1, check.names=FALSE)[,-19]
rpkmDFmiR <- rpkmDFmiR[rowMeans(rpkmDFmiR) > 50,]
d <- cor(rpkmDFmiR, method="spearman")
hc <- hclust(as.dist(1-d))
pdf("./results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()
```

![tree_pdf](https://github.com/sloanlewis/EE283-FINAL/blob/master/sample_tree_rlog_femaels.pdf "RPKM Tree")

```
#Group-wise and sample wise PCA following rlog transformation
library(DESeq2)
targets <- read.delim("targets.txt", comment.char = "#")
countDFmiR <- as.matrix(read.table("./results/countDFmiR.xls"))
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDFmiR, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("./results/PCA_group.pdf")
plotPCA(rld)
dev.off()
```

![pca_pdf](https://github.com/sloanlewis/EE283-FINAL/blob/master/PCA_group_14_CMH.pdf "Group PCA")

```
#Groupwise 3D PCA using rlog
library(scatterplot3d)
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
rld <- rlog(dds)
pca <- prcomp(t(assay(rld)), center=TRUE)
pca2 <- pca$x
groups <- targets$Factor
pdf("results/3D_miRNA.pdf")
scatterplot3d(pca2[,1], pca2[,2], pca2[,3], color = as.numeric(groups), pch=19, type='h')
dev.off()
```

![3d_pdf](https://github.com/sloanlewis/EE283-FINAL/blob/master/3D_14.pdf "3D PCA")


## Final wrap-up

This project allowed me to get this pipeline running smoothly all the way through as well as document all my steps so that when I get more miRNA sequencing for this project soon, I will be able to analyze them more efficiently without having to trouble shoot problems I have already solved. Before this project, I also did not have a lot of the shell scripts running on whole directories or submitting to the cluster. This will save me a lot of time in the future. 

One major issue I had with this project was figuring out which annotation file I should use for counting the reads because the macaca mulatta genome is not as well annotated as human. I tried using featureCounts with a miRBase GFF file, but the number of miRNAs counted was not quite as strong as when I used summarizeOverlaps and my ensembl GTF. Another issue has been that for some of the older plasma samples, I am struggling to get enough coverage of the genome because the alignment rates are much lower than for newer plasma. When I aligned the unmapped reads from those files to the whole GTF file, I found that about 30% were mRNAs suggesting that possibly some fragmented RNA is making it into my miRNA selection. I am still planning on trying to make this pipeline better in the future as new tools for miRNAs continue to be developed. 

Last, my major take away from the second part of this course was actually the ATACseq data analysis as that was the data set I chose to work on for many of the labs. I was hoping to have my ATACseq data to work on for this final project, but it is currently being sequenced. I chose to improve and update my miRNA pipeline for this project as it would be more useful for me to work on my own data. 

#done
