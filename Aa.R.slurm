#!/bin/bash -l
#SBATCH -A lp_dila
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --time=72:00:00
#SBATCH --cluster=wice
#SBATCH --cpus-per-task=1
#SBATCH -o out.R.s20.%J.out.txt
#SBATCH -e out.R.s20.%J.err.txt
#SBATCH --mem 124G

# we add num_threads 30 in the blast command and asked 2 nodes above instead of 1
source /user/leuven/301/vsc30196/.bashrc
eval "$(conda shell.bash hook)"
conda init bash
#conda env list
conda activate scbac 


echo "aALB015"
conda activate scbac
module load STAR
module load  SQLite



 echo "we were running /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/software/PRISM.s17.R we fixed a print noquote and USED DB_standard, added the module load  SQLite, and the wrote the blast command to a temp.file, we were running s18. with 30 nodes"


#R CMD BATCH aALB015.R
#R CMD BATCH aALB007.R
#R CMD BATCH  aV21.ALB002.R
R CMD BATCH  aV21.ALB024.R

#submission with: and the order is important
# sbatch --export=SAMPLE=sc5rALBR002_hg19  ALLv4.16.slurm

