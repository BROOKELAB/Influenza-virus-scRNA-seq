#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=140G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=joel3@illinois.edu
#SBATCH -J CR-count-cal07
#SBATCH --output=cell-ranger-%A_%a.out
#SBATCH --array=1-8%2
#SBATCH -D /home/labs/cbrooke_lab/Joel/sc_RNA_seq_2020/src/slurm_cal07


cd ../

Sample_ID=`head Cal07_files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
Sample_name=`head Cal07_files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`
File_folder=`head Cal07_files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/7.1.0

### specifying the output folder
cd ../results/cellranger_cal07/

# use ref made with renamed gtf to avoid parsing errors

DB=/home/groups/hpcbio_shared/cbrooke_lab/single_cell/scRNA_2020_Oct/
data/genome/cellranger_4.0.0_hg38_plus_cal0709_renamed

### Run app on file
cellranger count --id=$Sample_ID \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 128 \
    --transcriptome=$DB \
    --fastqs=$File_folder \
    --sample=$Sample_name \
    --expect-cells=5000
