Variant calling was first performed for each genotype using the  the script 'GATK_calling_per_genotype.sh'.

This script performs quality check and trimming of fastq files, the mapping on reference genome, duplicates removal.

Then was performed three rounds of aplotype calling, using GATK software, as described in GATK best practice:
-First round Indel Realignment without known VCF , then haplotype calling and hard filtering of SNP and INDELS
-Second round Indel Realignment with known VCF , followed with Base Quality Score Recalibration, then haplotype calling, genotype gVCG, and hard filtering of SNP and INDELS
-Third round Indel Realignment with known VCF , followed with Base Quality Score Recalibration, then haplotype calling.


Afterwards, all individual gVCF file where combined in a single multivariant gvcf file, on which was performed Genotyping and hard filtering of SNPs and INDELs using the script 'snp_indels_calling_on_gVCF.sh'
