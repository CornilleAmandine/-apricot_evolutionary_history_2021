#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=GA

# Job Array
#SBATCH --array=1-5%5
########SBATCH --array= 2-35%40
########SBATCH --array=2-3%2

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=12800MB

# Output (%j = jobid)
#SBATCH --output=/gpfs/scratch/sdecroocq/mapping/EU_apricot/Against_marouch_v3.1/Output/GATK_%j.out

# Email me
#SBATCH --mail-user=stephane.decroocq@inra.fr
#SBATCH --mail-type=ALL

################################################################################

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'TASK Id:' $SLURM_ARRAY_TASK_ID
echo 'Directory:' $(pwd)
# Detail Information:
#scontrol show job $SLURM_JOB_ID
echo '########################################################'


echo '########################################################'
echo '# Variables'

WDIR="/gpfs/scratch/sdecroocq/mapping/EU_apricot"
SEQ_BRUTES_ZIPPED="/gpfs/scratch/apricot_data_shared/Euapricot"
SEQ_TRAVAIL="/Sequence"
REF_GENOME="/gpfs/scratch/sdecroocq/Reference_genome/Parmeniaca_Marouch/Marouch_genome_v3.1.fasta"
NB_THREAD="5"

JAVA_MEM_OPT='-Xms8g'


echo '########################################################'
echo '# Get the sample'

#en fonction du nom car je ne sais pas traiter le SAMPLE avec les structures de nom differentes
#
SAMPLES=($(ls $SEQ_BRUTES_ZIPPED/INRA_Decroocq_121109_Run_Sample_*R1*.fastq.gz) )
NB_SAMPLES=$(ls $SEQ_BRUTES_ZIPPED/INRA_Decroocq_121109_Run_Sample_*R1*.fastq.gz| wc -l) #35 sequences
echo $NB_SAMPLES



# SLURM_ARRAY_TASK_ID start at 1
# Array SAMPLES start at 0
ID=$(( $SLURM_ARRAY_TASK_ID - 1 ))

SAMPLE_R1="${SAMPLES[$ID]}"
SAMPLE_R2="$(echo $SAMPLE_R1 | sed -e 's/_R1_/_R2_/g')"
#couperbasename à partir de _
#SAMPLE=$(basename $SAMPLE_R1 | cut -d '_' -f 1)
SAMPLE=$( basename $SAMPLE_R1 | cut -d '.' -f 1 | sed 's/^INRA_Decroocq_121109_Run_Sample_//g' | sed 's/_[A-Z]\{6\}_[A-Z][0-9]\{3\}_R1_001$//g' )

echo "Sample = $SAMPLE"
echo "Sample = $SLURM_ARRAY_TASK_ID / $NB_SAMPLES"
echo "Sample R1 = $SAMPLE_R1"
echo "Sample R2 = $SAMPLE_R2"

echo '########################################################'
echo '# Load the modules'
module load cutadapt/2.3
module load QC-toolkit/2.3.3
module load bwa/0.7.17
module load samtools/1.6
module load picard-tools/2.9.2
module load GATK/3.8
module load fastp/0.20.0

echo 'modules loaded'



trimming () {
echo '########################################################'
echo '# TRIMMING with Fastp'
# attention ${SAMPLE}_1 est different de $SAMPLE_1

fastp -i $SAMPLE_R1 -I $SAMPLE_R2 -o ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.trim.fastq -O ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_2.trim.fastq -h ${WDIR}${SEQ_TRAVAIL}/report/${SAMPLE}_fastp.html \
-j ${WDIR}${SEQ_TRAVAIL}/report/${SAMPLE}_fastp.json



echo "${SAMPLE} cleaning Finished"
}




sequencing-lane-identification () {
echo '########################################################'
echo '# sequencing-lane-identification'

#Then, we start to map our reads to the reference genome
# extraire la f3 et f4  pour servir d'input dans bwa (pour les SWAG)

lane=$(cat ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.trim.fastq | head -n 1 | cut -d : -f 3)
first=$(cat ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.trim.fastq | head -n 1 | cut -d : -f 4)

echo "lane = $lane"
echo "first = $first"


# extraire la f3 et f4  pour servir d'input dans bwa (pour les mume)
#lane=$(cat ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.fastq | head -n 1 | cut -d : -f 2)
#first=$(cat ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.fastq | head -n 1 | cut -d : -f 1)

#echo "$first" > ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}.txt
#run=$(cat ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}.txt | cut -d ' ' -f 2)

}






mapping () {
echo '########################################################'
echo '# MAPPING'

# ONCE (indexing)
#bwa index $REF_GENOME
#module load samtools/1.6
#samtools faidx $REF_GENOME
#module load picard-tools/2.9.2
#java -jar /module/apps/picard-tools/2.9.2/picard.jar CreateSequenceDictionary REFERENCE=/gpfs/scratch/sdecroocq/Reference_genome/Parmeniaca_Marouch/Marouch_genome_v3.1.fasta OUTPUT=/gpfs/scratch/sdecroocq/Reference_genome/Parmeniaca_Marouch/Marouch_genome_v3.1.dict
# ATTENTION: Il faut peutetre renommer Marouch_genome_v3.1.fasta.dict en Marouch_genome_v3.1.dict

bwa mem -t $NB_THREAD -R '@RG\tID:'${first}'\tPL:illumina\tLB:'${lane}'\tSM:'${SAMPLE}'_marouch_v3' $REF_GENOME ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_1.trim.fastq ${WDIR}${SEQ_TRAVAIL}/${SAMPLE}_2.trim.fastq | samtools view -S -b - > $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.bam
echo '# MAPPING DONE'
}



samtools-sort () {
echo '########################################################'
echo '# SAMTOOLS SORT'

########################################################
#then, use samtools to sort the sequence
samtools sort -@ $NB_THREAD -m 8G -O bam -o $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.bam $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.bam
echo '# SAMTOOLS SORT DONE'
}

remove-duplicates () {
echo '########################################################'
echo '#Remove the duplicates'

ulimit -c unlimited
java -Xmx8g -jar /module/apps/picard-tools/2.9.2/picard.jar MarkDuplicates \
  I=$WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.bam \
  O=$WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
  M=$WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.markdup_metrics.txt
echo '#Remove the duplicates DONE'

}


samtools-index () {
echo '########################################################'
echo '#Make a index for .sorted.markdup.bam'
samtools index $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam
}


gatk-round1-realign-without-VCF () {
echo '########################################################'
echo '#Indel Realignment, without known vcf'

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
 -o $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
 -o $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
 --targetIntervals $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals
echo '#Indel Realignment, without known vcf DONE'
}

gatk-round1-HaplotypeCaller () {
echo '#######################################################'
echo '# use HaplotypeCaller yo test the variant'

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -nct $NB_THREAD \
  -T HaplotypeCaller \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
  --emitRefConfidence GVCF \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.g.vcf
echo '# use HaplotypeCaller yo test the variant DONE'
echo '#######################################################'
echo '#use GenotypeGVCFs to finish the variant calling'

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R $REF_GENOME \
  --variant $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.g.vcf \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf
echo '#use GenotypeGVCFs to finish the variant calling DONE'
}

gatk-round1-Hard-Filter-SNP () {
echo '#######################################################'
echo '#Hard filter the snps RND1'

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf \
  -selectType SNP \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_snps.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_snps.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "my_snp_filter" \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_snps.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_snps.vcf \
  -ef -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf



echo '#######################################################'
echo '#SNP Hard filter RND1 finished'
}


gatk-round1-Hard-Filter-Indels () {
echo '#######################################################'
echo '#Hard filter the indels'

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf \
  -selectType INDEL \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_indels.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_indels.vcf \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "my_indel_filter" \
  -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_indels.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_indels.vcf \
  -ef -o $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf

echo '#######################################################'
echo "Indel Hard filter finished RND1"
}



gatk-round2 () {
echo '#######################################################'
echo '# second round of GATK'

#Indel Realignment, with known vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
 -known $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
 -known $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
 -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
 -known $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
 -known $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
 -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
 --targetIntervals $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals

#BQSR  Base Quality Score Recalibration

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T BaseRecalibrator \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
 --knownSites $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
 --knownSites $WDIR/Against_marouch_v3.1/0_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
 -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.recal_data.report.grp

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T PrintReads \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
 -rf BadCigar \
 --BQSR $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.recal_data.report.grp \
 -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.BQSR.bam


#######################################################
# use HaplotypeCaller yo test the variant

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
   -nct $NB_THREAD \
  -T HaplotypeCaller \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.BQSR.bam \
  --emitRefConfidence GVCF \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.g.vcf

#######################################################
#use GenotypeGVCFs to finish the variant calling
java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R $REF_GENOME \
  --variant $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.g.vcf \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf

#Hard filter the snps and indels
java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf \
  -selectType SNP \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_snps.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_snps.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "my_snp_filter" \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_snps.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_snps.vcf \
  -ef -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf

echo "SNP Hard filter finished RND2"

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/${SAMPLE}_marouch_v3.HC.vcf \
  -selectType INDEL \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_indels.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.raw_indels.vcf \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "my_indel_filter" \
  -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_indels.vcf

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R $REF_GENOME \
  -V $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.filtered_indels.vcf \
  -ef -o $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf

echo "Indel Hard filter RND2 finished"
}



gatk-round3 () {
echo '#######################################################'
echo 'third round of GATK'

#Indel Realignment, with known vcf
java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
  -known $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
  -known $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
  -o $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/0_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.bam \
  -known $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
  -known $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
  -o $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
  --targetIntervals $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.IndelRealigner.intervals

#BQSR  Base Quality Score Recalibration

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
 -T BaseRecalibrator \
 -R $REF_GENOME \
 -I $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
 --knownSites $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.snps.vcf \
 --knownSites $WDIR/Against_marouch_v3.1/1_vcf_calling/Hard_filter/${SAMPLE}_marouch_v3.indels.vcf \
 -o $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.recal_data.report.grp

java -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.bam \
  -rf BadCigar \
  --BQSR $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.recal_data.report.grp \
  -o $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.BQSR.bam

#######################################################
# use HaplotypeCaller yo test the variant

java $JAVA_MEM_OPT -jar /module/apps/GATK/3.8/GenomeAnalysisTK.jar \
   -nct $NB_THREAD \
  -T HaplotypeCaller \
  -R $REF_GENOME \
  -I $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.sorted.markdup.realign.BQSR.bam \
  --emitRefConfidence GVCF \
  -o $WDIR/Against_marouch_v3.1/2_vcf_calling/${SAMPLE}_marouch_v3.g.vcf
echo 'third round of GATK DONE'
}



#######################################################
# MAIN
#######################################################

trimming
sequencing-lane-identification
mapping
samtools-sort
remove-duplicates
samtools-index
gatk-round1-realign-without-VCF
gatk-round1-HaplotypeCaller
gatk-round1-Hard-Filter-SNP
gatk-round1-Hard-Filter-Indels
gatk-round2
gatk-round3


#######################################################
echo "Job finished ($(date --iso-8601=seconds))"
