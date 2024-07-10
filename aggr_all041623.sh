#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=150G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=joel3@illinois.edu
#SBATCH -J aggr-all
#SBATCH --output=aggr-all-%j.out
#SBATCH -D /home/labs/cbrooke_lab/Joel/sc_RNA_seq_2020/src/slurm_aggrall


### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load cellranger/7.1.0

### specifying the output folder
cd /home/labs/cbrooke_lab/Joel/sc_RNA_seq_2020/results

cellranger aggr --id=Aggregate_all \
                  --csv=../src/AggrList_all.csv \
                  --normalize=none
