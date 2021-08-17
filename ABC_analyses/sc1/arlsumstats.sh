#!/bin/bash

#####################################################################ARLSUMSTAT
#Je change le gametic phase dans le fichier temp.arp
#sed -i "s/GameticPhase=0/GameticPhase=1/g" div_carthu-temp/div_carthu-temp_1_1.arp
#Calcul des stats
#./arlsumstat cbp-temp/cbp-temp_1_1.arp summary_stats_temp.txt 0 1
./arlsumstat 3gp_apricot-temp/3gp_apricot-temp_1_1.arp 3gp_apricot.obs 0 1

#./JFS_boucle.sh

hdDNA=$(head -n1 '3gp_apricot.obs')
ssDNA=$(tail -n +2 '3gp_apricot.obs')

hdssstat_bef=$(head -n1 'summary_stats_temp6.txt')
ssstat_bef=$(tail -n +2 'summary_stats_temp6.txt')

echo -e "${hdDNA}${hdssstat_bef}" > summary_stats_temp.txt
echo -e "${ssDNA}${ssstat_bef}" >> summary_stats_temp.txt

#rm summary_stats_temp.txt
rm summary_stats_temp1.txt
rm summary_stats_temp2.txt
rm summary_stats_temp3.txt
#rm summary_stats_temp3bis.txt
#rm summary_stats_temp4.txt
#rm summary_stats_temp5.txt
rm summary_stats_temp6.txt