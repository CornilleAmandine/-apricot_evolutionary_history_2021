#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=struct_NonMume555Maf0_05

### Limit run time "hours:minutes:seconds" (default: 365 days)
#SBATCH --time=999:00:00
#SBATCH --output=/gpfs/scratch/sdecroocq/Structure_2020/NonMume555_AllVariants/NonMume555_AllVariantsMaf0_05_%A.out


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



#Before starting we must define variables


ORIGIN_VCF="/scratch/sdecroocq/Structure_2020/data_preparation/jan_2020/apricot_collection_2019_marouch_v3_clean.vcf"
WORKING_DIR="/gpfs/scratch/sdecroocq/Structure_2020/NonMume555_AllVariants"

cd $WORKING_DIR
#prepare the name of the VCF file which is used as a root name for the analysis = SAMPLE
SAMPLE="NonMume555"
echo 'SAMPLE = '${SAMPLE}


echo '########################################'
echo 'load the modules'
module load bcftools/1.6
module load vcftools/0.1.16
module load plink/1.90
echo '########################################'
echo 'load the modules done'

ind_selection_NonMume555() {

echo '########################################'
echo 'ind_selection_2_rnd'
echo '########################################'

#make groups of predefined individuals (NonMume_IndToKeep.txt)

vcftools --vcf ${ORIGIN_VCF} --keep NonMume_IndToKeep.txt --max-missing 0.85 --maf 0.05 --mac 3 --minQ 30  --recode --recode-INFO-all --out $WORKING_DIR/${SAMPLE}_555_AllVariants


cp ${SAMPLE}_555_AllVariants.recode.vcf ${SAMPLE}_555_AllVariantsMaf0_05.vcf



echo '########################################' 
echo 'ind_selection_NonMume555 done'
echo '########################################'

}
vcf_to_plink() {
echo '########################################' 
echo 'vcf_to_plink'
echo '########################################'
vcftools --vcf ${SAMPLE}_555_AllVariantsMaf0_05.vcf --plink --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05
echo '########################################' 
echo 'vcf_to_plink done'
echo '########################################'
}


prepare_format1() {

echo '########################################' 
echo 'prepare_format1'
echo '########################################'

# 50 5 0.0428 (50 for (window size in variant count or kilobase (if the 'kb' modifier is present)) 5 for (a variant count to shift the window at the end of each step) 0.0428 for (pairwise r2 threshold: at each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain) )

plink --file $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05 --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp --indep-pairwise 50 5 0.0428

#--extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis
#--make-bed creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations 

plink --file $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05 --extract $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp.prune.in --make-bed --out $WORKING_DIR/b${SAMPLE}_555_AllVariantsMaf0_05_ldp 

#The --bfile flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced
#--recode creates a new text fileset, after applying sample/variant filters and other operations
#The 'tab' modifier makes the output mostly tab-delimited instead of mostly space-delimited when the format permits both delimiters


plink --bfile $WORKING_DIR/b${SAMPLE}_555_AllVariantsMaf0_05_ldp  --recode tab --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp

plink --file $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp --make-bed --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp
echo '########################################' 
echo 'prepare_format1 done'
echo '########################################'

}
prepare_sort_indiv() {

echo '########################################' 
echo 'prepare_sort_indiv'
echo '########################################'



# reorder samples in a PLINK dataset, per an order specified in a second file
#This allows you to specify how samples should be sorted when generating new datasets

#--bfile [original fileset] --indiv-sort f [file describing new order] --make-bed --out [new prefix]

plink --bfile $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp --indiv-sort f $WORKING_DIR/NonMume_IndOrder.txt --make-bed --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp_sort

plink --file $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp_sort --make-bed --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp

plink --bfile $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp_sort --recode tab --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_ldp_sorted

plink --file ${SAMPLE}_555_AllVariantsMaf0_05_ldp_sorted --make-bed --out $WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp

echo '########################################' 
echo 'prepare_sort_indiv done'
echo '########################################'
echo '########################################' 
echo 'reordered file for structure '_555_AllVariantsMaf0_05_sorted_ldp''
echo 'prepare_format done'
echo '########################################'

}


fastructure() {
echo '########################################'
echo 'fastructure'
echo '########################################'
#run structure
module load fastStructure/1.0
structure -K 2 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 3 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 4 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 5 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 6 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 7 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 8 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 9 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 10 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 11 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 12 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 13 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 14 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100
structure -K 15 --input=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp  --output=$WORKING_DIR/${SAMPLE}_555_AllVariantsMaf0_05_sorted_ldp --full --seed=100

echo '########################################'
echo 'fastructure done'
echo '########################################'
}



#######################################################
# MAIN
#######################################################
ind_selection_NonMume555
vcf_to_plink
prepare_format1
prepare_sort_indiv
fastructure

#######################################################
echo 'Job finished'
