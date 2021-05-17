#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=all_struct

### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=999:00:00
#SBATCH --output=/gpfs/scratch/sdecroocq/Structure_2020/output/structure_collection_2019_%A.out


### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=12800MB

################################################################################

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
# Detail Information:
scontrol show job $SLURM_JOB_ID
echo '########################################'

#Preparation of the starting snp file for the whole population (containing all variants, maf 0.005). 
#This file will be later subset in smaller population and filtered to the requiered need of each analysis

#Fichier VCF de départ

ORIGIN_VCF="/gpfs/scratch/sdecroocq/mapping/EU_apricot/Against_marouch_v3.1/combine_g_vcf/apricot_collection_2019_marouch_v3.1.snps.vcf"
WORKING_DIR="/gpfs/scratch/sdecroocq/Structure_2020"


#go to working directory

cd $WORKING_DIR/data_preparation

#prepare the name of the VCF file which is used as a root name for the analysis = SAMPLE
SAMPLE=$(basename $ORIGIN_VCF | cut -d '.' -f 1 )
echo 'SAMPLE = '${SAMPLE}


echo '########################################'
echo 'load the modules'

module load bcftools/1.6
module load vcftools/0.1.16
module load plink/1.90


echo '########################################'
echo 'load the modules done'



prepare_bia() {

echo '########################################'
echo 'prepare_biallelic file'



bcftools view -m2 -M2 -v snps $ORIGIN_VCF  -o bia_${SAMPLE}.vcf







cd $WORKING_DIR/data_preparation/jan_2020

mv $WORKING_DIR/data_preparation/bia_${SAMPLE}.vcf ./bia_${SAMPLE}_1to8.vcf

echo '########################################'
echo 'prepare_biallelic file done'

}


formating_filtering() {


# We are going to only keep variants that have been successfully, ...
# genotyped in 85% of individuals (--max-missing 0.85), 
# a minimum quality score of 30 (--minQ 30), 
# and a minor allele count of 3 (--mac 3).
# MAF 0.005 (--maf 0.005) 
echo '########################################'
echo 'data filtering'


vcftools --vcf bia_${SAMPLE}_1to8.vcf --max-missing 0.85 --maf 0.005 --mac 3 --minQ 30 --recode --recode-INFO-all --out ${SAMPLE}_clean.vcf
mv ${SAMPLE}_clean.vcf.recode.vcf ${SAMPLE}_clean.vcf

echo '########################################'
echo 'data filtering done'



}






#######################################################
# MAIN
#######################################################
prepare_bia
formating_filtering

#######################################################
echo 'Job finished'
