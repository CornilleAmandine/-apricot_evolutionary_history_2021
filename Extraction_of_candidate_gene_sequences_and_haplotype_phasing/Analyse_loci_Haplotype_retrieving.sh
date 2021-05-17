#!/bin/bash

################################ Slurm options #################################

#SBATCH --job-name=Analyse_loci_Haplotype_retrieving

# Working Directory
#SBATCH --workdir=/scratch/apricot_data_shared/CG_European/mapped_loci

# Merge ouput and error files
#SBATCH --output=/home/ag/Analyse_loci_Haplotype_retrieving_%j.out

# Mail job status at the start and end of a job
#SBATCH --mail-user=alexis.groppi@u-bordeaux.fr
#SBATCH --mail-type=ALL

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=125000MB


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

# modules loading
module load vcftools/0.1.16
module load samtools/1.9
module load htslib/1.9
module load java/1.8.0_72
module load bcftools/1.9

export LIBRARY_PATH=/home/ag/vcflib/lib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/ag/vcflib/lib/:$LD_LIBRARY_PATH


echo "#############################" 
echo "conversion g.vcf to vcf" `date`
set -x
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
INPUT_DIR=/scratch/apricot_data_shared/CG_European/Gvcf_Input
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_g.vcf_list_temp | tail -1)
SAMPLE_NAME=$(sed 's/.\{17\}$//' <<< "$file" | sed 's/Kz/KZ/' | sed 's/US/US_/')
vcftools --vcf $INPUT_DIR/$file --min-alleles 2 --recode --recode-INFO-all --out $OUT_DIR/$SAMPLE_NAME


echo "#############################" 
echo "filter vcf" `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_recode.vcf_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{11\}$//' <<< "$file");
/home/ag/vcflib/bin/vcffilter -f "MLEAC > 0" \
$OUT_DIR/$file | \
grep -vw "\./\." | \
grep -vw "\./[0-9]" | \
grep -vw "[0-9]/\." | \
grep -vw "[0-9]/[3-9]" | \
grep -vw "[3-9]/[0-9]" | \
grep -vw "[0-9]|[3-9]" | \
grep -vw "[3-9]|[0-9]" \
> $OUT_DIR/$SAMPLE_NAME"_filtered.vcf";



echo "#############################" 
echo "extract haplotype-relevant information " `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_filtered.vcf_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{13\}$//' <<< "$file");
/home/ag/HapCUT2/build/extractHAIRS \
--bam /scratch/apricot_data_shared/CG_European/mapping/$SAMPLE_NAME.sorted.bam \
--ref /mnt/cbib/sWAGMAN/Marouch/Assemblages/Marouch_genome_v3.1.fasta \
--indels 1 \
--VCF $OUT_DIR/$file \
--out $OUT_DIR/$SAMPLE_NAME"_frag_file";



echo "#############################" 
echo "assemble fragment file into haplotype blocks" `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_frag_file_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{10\}$//' <<< "$file");
/home/ag/HapCUT2/build/HAPCUT2 \
--fragments $OUT_DIR/$file \
--VCF $OUT_DIR/$SAMPLE_NAME"_filtered.vcf" \
--output $OUT_DIR/$SAMPLE_NAME"_haplotype_out";



echo "#############################" 
echo "HapCut to vcf " `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_haplotype_out_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{14\}$//' <<< "$file");
java -jar /home/ag/fgbio/target/scala-2.13/fgbio-1.2.0-abc7fd5-SNAPSHOT.jar HapCutToVcf \
-v $OUT_DIR/$SAMPLE_NAME"_filtered.vcf" \
-i $OUT_DIR/$file \
-o $OUT_DIR/$SAMPLE_NAME"_haplotype.vcf";



echo "#############################" 
echo "compression and indexation of vcf " `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input_haplotype.vcf_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{14\}$//' <<< "$file");
bcftools view $OUT_DIR/$file -Oz -o $OUT_DIR/$SAMPLE_NAME.vcf.gz;
bcftools index $OUT_DIR/$SAMPLE_NAME.vcf.gz;



echo "#############################" 
echo "extraction of fasta haploid 1 and 2 and extraction of locus specific vcf" `date`
OUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
file=$(head -n $SLURM_ARRAY_TASK_ID $OUT_DIR/input.vcf.gz_list_temp | tail -1);
SAMPLE_NAME=$(sed 's/.\{7\}$//' <<< "$file");
while IFS=, read -r Locus_name Position;
  do
  mkdir -p $OUT_DIR/$Locus_name #make directory only if it doesn't exist
  samtools faidx /mnt/cbib/sWAGMAN/Marouch/Assemblages/Marouch_genome_v3.1.fasta $Position | bcftools consensus -H 1 $OUT_DIR/$SAMPLE_NAME.vcf.gz -o $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P1_temp.fasta;
  samtools faidx /mnt/cbib/sWAGMAN/Marouch/Assemblages/Marouch_genome_v3.1.fasta $Position | bcftools consensus -H 2 $OUT_DIR/$SAMPLE_NAME.vcf.gz -o $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P2_temp.fasta;
  sed 's/>.*/& '$SAMPLE_NAME' '$Locus_name' P1/' $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P1_temp.fasta > $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P1.fasta
  sed 's/>.*/& '$SAMPLE_NAME' '$Locus_name' P2/' $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P2_temp.fasta > $OUT_DIR/"$SAMPLE_NAME"_"$Locus_name"_P2.fasta
  bcftools filter -r $Position $OUT_DIR/$SAMPLE_NAME.vcf.gz -o $OUT_DIR/$Locus_name/"$SAMPLE_NAME"_"$Locus_name".vcf
  done < $OUT_DIR/Locus_position_500_CG_European_on_Marouch_v3_1_Dec_2020.txt


echo "#############################" 
echo "Cleaning intermediate files " `date`
rm $OUT_DIR/$SAMPLE_NAME.recode.vcf
rm $OUT_DIR/$SAMPLE_NAME"_filtered.vcf"
rm $OUT_DIR/$SAMPLE_NAME"_frag_file"
rm $OUT_DIR/$SAMPLE_NAME"_haplotype_out"
rm $OUT_DIR/$SAMPLE_NAME"_haplotype.vcf"
rm $OUT_DIR/"$SAMPLE_NAME"_*_P1_temp.fasta
rm $OUT_DIR/"$SAMPLE_NAME"_*_P2_temp.fasta
rm $OUT_DIR/$SAMPLE_NAME.vcf.gz
rm $OUT_DIR/$SAMPLE_NAME.vcf.gz
rm $OUT_DIR/$SAMPLE_NAME.vcf.gz.csi




echo "Fin traitement : "`date`


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)



