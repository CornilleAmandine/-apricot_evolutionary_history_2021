#!/bin/bash

################################ Slurm options #################################

#SBATCH --job-name=Interleave_and_Merge_Fasta_Loci

# Working Directory
#SBATCH --workdir=/scratch/apricot_data_shared/CG_European/mapped_loci

# Merge ouput and error files
#SBATCH --output=/home/ag/Interleave_and_Merge_Fasta_Loci_%j.out

# Mail job status at the start and end of a job
#SBATCH --mail-user=alexis.groppi@u-bordeaux.fr
#SBATCH --mail-type=ALL

### Specify requirements - Task (default: 1 node, 1 Core, 12.5G mem/cpu)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


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
module load samtools/1.9
module load fastxtoolkit/0.0.13


echo "#############################" 
echo "Generation of Marouch reference for all loci" `date`
INPUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
OUT_DIR=/scratch/apricot_data_shared/Interleaved_Fasta_CG_European
while IFS=, read -r Locus_name Position;
do
samtools faidx /mnt/cbib/sWAGMAN/Marouch/Assemblages/Marouch_genome_v3.1.fasta $Position -o $OUT_DIR/Marouch_"$Locus_name"_Ref_temp.fasta;
sed 's/>.*/& Marouch '$Locus_name' Reference/' $OUT_DIR/Marouch_"$Locus_name"_Ref_temp.fasta > $OUT_DIR/Marouch_"$Locus_name"_Ref.fasta;
done < $INPUT_DIR/Locus_position_500_CG_European_on_Marouch_v3_1_Dec_2020.txt
rm $OUT_DIR/*_temp.fasta


echo "#############################"
echo "Copy fasta files CG_European" `date`
INPUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
OUT_DIR=/scratch/apricot_data_shared/Interleaved_Fasta_CG_European
cp -av $INPUT_DIR/*_P*.fasta $OUT_DIR/


echo "#############################"
echo "Creation of loci directories" `date`
INPUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
OUT_DIR=/scratch/apricot_data_shared/Interleaved_Fasta_CG_European
while IFS=, read -r Locus_name Position;
do
mkdir $OUT_DIR/"$Locus_name";
done < $INPUT_DIR/Locus_position_500_CG_European_on_Marouch_v3_1_Dec_2020.txt


echo "#############################"
echo "Interleave P1/P2 and merge all fasta for the same locus with the Marouch Reference" `date`
INPUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
OUT_DIR=/scratch/apricot_data_shared/Interleaved_Fasta_CG_European
while IFS=, read -r Locus_name Position;
do
cat $OUT_DIR/*_"$Locus_name"_P1.fasta > $OUT_DIR/"$Locus_name"_P1_temp.fasta;
cat $OUT_DIR/*_"$Locus_name"_P2.fasta > $OUT_DIR/"$Locus_name"_P2_temp.fasta;
paste <(awk '/^>/ { printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' $OUT_DIR/"$Locus_name"_P1_temp.fasta ) \
<(awk '/^>/ { printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' $OUT_DIR/"$Locus_name"_P2_temp.fasta ) | \
tr "\t" "\n" | sed '/^$/d' > $OUT_DIR/"$Locus_name"_P1_P2_temp.fasta;
cat $OUT_DIR/Marouch_"$Locus_name"_Ref.fasta $OUT_DIR/"$Locus_name"_P1_P2_temp.fasta | fasta_formatter -w 60 -o $OUT_DIR/$Locus_name/"$Locus_name".fasta;
done < $INPUT_DIR/Locus_position_500_CG_European_on_Marouch_v3_1_Dec_2020.txt
rm $OUT_DIR/*_temp.fasta

echo "#############################"
echo "Copy fasta files" `date`
INPUT_DIR=/scratch/apricot_data_shared/Interleaved_Fasta_CG_European
OUT_DIR=/mnt/cbib/raw_apricot_data_shared/CG_European_mapped_loci
for dir in $(ls -l $INPUT_DIR | grep ^d | tr -s ' ' ' ' |cut -d' ' -f9)
do
cp -av $INPUT_DIR/$dir/*.fasta $OUT_DIR/$dir/
done

echo "#############################" 
echo "Cleaning intermediate files" `date`
INPUT_DIR=/scratch/apricot_data_shared/CG_European/mapped_loci
rm $INPUT_DIR/*_P1.fasta
rm $INPUT_DIR/*_P2.fasta
rm $INPUT_DIR/input_*_temp



echo "Fin traitement : "`date`

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

