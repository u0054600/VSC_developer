#!/bin/bash -l
#SBATCH -A lp_dila
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --time=72:00:00
#SBATCH --cluster=wice
#SBATCH --cpus-per-task=1
#SBATCH -o out.s21.%J.out.txt
#SBATCH -e out.s21.s18.%J.err.txt
#SBATCH --mem 124G

# we add num_threads 30 in the blast command and asked 2 nodes above instead of 1
source /user/leuven/301/vsc30196/.bashrc
eval "$(conda shell.bash hook)"
conda init bash
#conda env list
conda activate scbac 


echo $SAMPLE
conda activate scbac
module load STAR
module load  SQLite

sample=${SAMPLE}

Rscript /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/software/PRISM.s21.R \
    --sample ${sample}  \
    --data_path /user/leuven/301/vsc30196/vrc/Projects//Project_Microbiome/Run/samples/${sample}/ \
    --kraken_path /lustre1/project/stg_00080/vrc/software/BBCondaS/envs/scbac/bin/kraken2 \
    --kraken_db_path /lustre1/project/stg_00080/vrc/Projects/Project_Microbiome/DB//k2_standard_20240605 \
    --seqkit_path /staging/leuven/stg_00080/vrc/software/BBCondaS/envs/scbac/bin/seqkit \
    --star_path  /data/leuven/software/biomed/wice_icelake_rocky8/2021a/software/STAR/2.7.11a-GCC-10.3.0/bin/STAR \
    --star_genome_dir /lustre1/project/stg_00080/vrc/REF/Gencode.v47/STAR/ \
    --blast_path /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/software/ncbi-blast-2.16.0+/bin/blastn \
    --blast_db_path /lustre1/project/stg_00080/vrc/Projects/Project_Microbiome/DB/DB_nt/nt \
    --prism_path /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/software/ \
    --paired T \
    --fq1_end .R1.fq.gz \
    --fq2_end .R2.fq.gz \
    --barcode_only T  \
    --max_sample 1000 \
    --min_read_per 10000 \
    --min_uniq_frac 5 --model_org_taxids /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/DB/Blast/model_org_taxids.txt \


 echo "we were running /user/leuven/301/vsc30196/vrc/Projects/Project_Microbiome/software/PRISM.s17.R we fixed a print noquote and USED DB_standard, added the module load  SQLite, and the wrote the blast command to a temp.file, we were running s18. with 30 node  V21"
#submission with: and the order is important
# sbatch --export=SAMPLE=sc5rALBR002_hg19  ALLv4.16.slurm

