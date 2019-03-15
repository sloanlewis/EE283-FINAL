#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00     # 2 day and 15 minutes
#SBATCH --output=counts_test
#SBATCH --mail-user=sloanal@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="small rna counting"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


# Print current date
date


Rscript small-RNA/small_RNA_readcounting_cluster.R


date
