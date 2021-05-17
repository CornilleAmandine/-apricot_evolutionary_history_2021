#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=combine

### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=999:00:00
#SBATCH --output=/gpfs/scratch/sdecroocq/mapping/EU_apricot/Against_marouch_v3.1/Output/combine_%A.out


### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=12500MB

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

echo '########################################################'
echo '# Variables'

# ENTER DIRECTORY AND FILE NAME FOR OUTPUT!!!!!!!!!!!!!!!!!!!!!!

WDIR="/gpfs/scratch/sdecroocq/mapping/EU_apricot/Against_marouch_v3.1/combine_g_vcf"
REF_GENOME="/gpfs/scratch/sdecroocq/Reference_genome/Parmeniaca_Marouch/Marouch_genome_v3.1.fasta"
NB_THREAD="1"
JAVA_MEM_OPT="-Xms12g"
COMBINED_FILE="apricot_collection_2019_marouch_v3.1"



cd $WDIR

#######################################################
#use GenotypeGVCFs to finish the variant calling


combine_files() {

echo 'combine_files'
echo '########################################'
module load GATK/3.8
java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T CombineGVCFs \
  -R $REF_GENOME \
--variant /gpfs/scratch/sdecroocq/mapping/EU_apricot/Against_marouch_v3.1/combine_g_vcf/EU_apricot_marouch_v3.1.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/CH_apricot/Against_marouch_v3.1/combine_g_vcf/CH_marouch_v3.1.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/wild_armeniaca/Against_marouch_v3.1/combine_g_vcf/wild_armeniaca_marouch_v3.1.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_1.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_2.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_3.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_4.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_5.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_6.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/SWAG/Against_marouch_v3.1/combine_g_vcf/SWAG_marouch_v3.1_7.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_7.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_6.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_5.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_4.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_3.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_2.g.vcf \
--variant /gpfs/scratch/sdecroocq/mapping/Against_marouch_v3.1/combine_g_vcf/MUME_marouch_v3.1_1.g.vcf \
 -o $WDIR/$COMBINED_FILE.g.vcf

echo 'combine_files done'
echo '########################################'
}

genotype_g_vcf_file() {

echo 'genotype_g_vcf_file'
echo '########################################'
module load GATK/3.8
java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R $REF_GENOME \
--variant $WDIR/$COMBINED_FILE.g.vcf \
 -o $WDIR/$COMBINED_FILE.vcf

echo 'genotype_g_vcf_file done'
echo '########################################'
}

SNP_calling_and_hard_filter() {

#Hard filter the snps and indels
echo 'SNP_calling_and_hard_filter'
echo '########################################'

module load GATK/3.8
java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T SelectVariants \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.vcf \
    -selectType SNP \
    -o $WDIR/$COMBINED_FILE.raw_snps.vcf

java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T VariantFiltration \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.raw_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o $WDIR/$COMBINED_FILE.filtered_snps.vcf

java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T SelectVariants \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.filtered_snps.vcf \
    -ef -o $WDIR/$COMBINED_FILE.snps.vcf

echo 'SNP_calling_and_hard_filter done'
echo '########################################'
}

INDEL_calling_and_hard_filter() {
echo 'INDEL_calling_and_hard_filter'
echo '########################################'

java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T SelectVariants \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.vcf \
    -selectType INDEL \
    -o $WDIR/$COMBINED_FILE.raw_indels.vcf

java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T VariantFiltration \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.raw_indels.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filterName "my_indel_filter" \
    -o $WDIR/$COMBINED_FILE.filtered_indels.vcf

java -Xmx8g -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
    -T SelectVariants \
      -R $REF_GENOME \
    -V $WDIR/$COMBINED_FILE.filtered_indels.vcf \
    -ef -o $WDIR/$COMBINED_FILE.indels.vcf


echo 'INDEL_calling_and_hard_filter done'
echo '########################################'
}

#######################################################
# MAIN
#######################################################
#combine_files
#genotype_g_vcf_file
SNP_calling_and_hard_filter
INDEL_calling_and_hard_filter
#######################################################
echo 'Job finished'SNP_calling_and_hard_filter
