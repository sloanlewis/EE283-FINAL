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

#Finally trim the illumina adapters and reduce the length from 18-24 bp
trim_galore --clip_R1 1 --length 17 --fastqc $fq
