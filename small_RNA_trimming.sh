#!/bin/bash -l

fq=$1
shift
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=alignment.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="small RNA trimming - QIASeq"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

#Load the packages
module load trim_galore
module load fastqc

first="_trimmed.fq.gz"
second="_trimmed_trimmed.fq.gz"
third="_trimmed_trimmed_trimmed.fq.gz"
fourth="_trimmed_trimmed_trimmed_trimmed.fq.gz"

#Trim the small rna adapters. Reads shorter than 18 bp are removed
trim_galore --no_report_file -q 20 --small_rna $fq

#Trim the 3' adapter.
trim_galore --no_report_file -q 20 --length 0 -a AACTGTAGGCACCATCAAT $fq$first

#Trims the 5' adapter.
trim_galore -q 20 --no_report_file --length 0 -a GATCGTCGGACTGTAGAACTCTGAAC $fq$second

#Finally trim the illumina adapters and reduce the length from 18-24 bp
trim_galore -q 20 --length 18 --max_length 24 --illumina --fastqc $fq$third
