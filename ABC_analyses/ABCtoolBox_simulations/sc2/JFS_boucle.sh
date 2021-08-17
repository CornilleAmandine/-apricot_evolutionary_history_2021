#!/bin/bash

####BUG fastsimcoal
#head -71 cbp-temp.par > cbp-temp.par

#####################################################################JOINT FREQUENCY SPECTRUM

#####################################################################POP_0_Y

####################################MISE EN FORME POUR CALCUL
#supprimer toutes les lignes avec 1 observation

#Je prépare mon fichier .obs pour les calculs A ADAPTER EN FONCTION DU NOMBRE DE MARQUEURS DE LA POP EN COLONNE (toujours nombre de marqueurs (ici 18) +3 (ligne "1 observations", ligne des d, ligne en plus du nombre de marqueurs, donc 3 en tout)

#LOOPS

#NO=152
#N1=152
#N2=84
#N3=84
#N4=286
#N5=286

SCENARIO_ID="3gp_apricot"

for i in 1_0 2_0 2_1
            do 
            #supprime les deux premières lignes du fichier
            sed '1,2d' ${SCENARIO_ID}-temp/${SCENARIO_ID}-temp_jointDAFpop$i'.obs' > JFS-temp_2.txt
            #supprime la première colonne du fichier
            cut -f 2- JFS-temp_2.txt >JFS-temp_3.txt
            #rm JFS-temp.txt
            #je retire la colonne 1 - imprime chaque ligne du fichier après avoir effacé le premier champs, $0 correspond à l'enregistrement complet
            #awk -F "\t" '{ $1 = "" ; print $0 }' JFS-temp_2.txt > JFS-temp_3.txt

######################################### SUMMARY STATISTICS CALCULEES A PARTIR DU JOINT FREQUENCY SPECTRUM

#############POLYMORPHISME CHEZ 0 Sx0_0_1 (lignes) (= W2, Tellier 2011)

###SOMMER TOUTE LA LIGNE 1 OK
sum_first_line=$(awk 'NR==1 {for(i=1;i<=NF;i++) t+=$i; print t; t=0}'  JFS-temp_3.txt)
echo $sum_first_line

###SOMMER LA DERNIERE LIGNE OK
sum_last_line=$(awk 'END {for(i=1;i<=NF;i++) t+=$i; print t; t=0}'  JFS-temp_3.txt)
echo $sum_last_line

#somme de la première colonne OK
sum_first_column=$(awk '{s+=$1} END {print s}' JFS-temp_3.txt)
echo $sum_first_column

#somme de la dernière colonne OK
sum_last_column=$(awk '{s+=$NF} END {print s}' JFS-temp_3.txt)
echo $sum_last_column

#VALEUR PREMIERE LIGNE, PREMIERE COLONNE OK
firstLINE_firstCOL=$(awk 'NR==1 {S1+=$1} END{print S1}' JFS-temp_3.txt)
echo $firstLINE_firstCOL

#VALEUR DERNIERE LIGNE, PREMIERE COLONNE OK
lastLINE_firstCOL=$(awk 'END {S1+=$1} END{print S1}' JFS-temp_3.txt)
echo $lastLINE_firstCOL

#VALEUR PREMIERE LIGNE, DERNIERE COLONNE
firstLINE_lastCOL=$(awk 'NR==1 {S1+=$NF} END{print S1}' JFS-temp_3.txt)
echo $firstLINE_lastCOL

#VALEUR DERNIERE LIGNE, DERNIERE COLONNE
lastLINE_lastCOL=$(awk 'END {S1+=$NF} END{print S1}' JFS-temp_3.txt)
echo $lastLINE_lastCOL

#tout le tableau
ALL_table=$(awk '{for (i=1; i<=NF; i++) s=s+$i}; END{print s}' JFS-temp_3.txt)
echo $ALL_table

#########################POLYMORPHISME CHEZ 0 Sx0_0_1 (lignes) (= W2, Tellier 2011) Private polymorphism in species 0

pre_sum=$(echo $firstLINE_firstCOL+$lastLINE_firstCOL+$firstLINE_lastCOL+$lastLINE_lastCOL | bc)
pre_sum2=$(echo $sum_first_line+$sum_last_line | bc)
Sx0_0_y=$(echo $pre_sum2-$pre_sum | bc)

echo $Sx0_0_y

#############POLYMORPHISME CHEZ 1 Sx1_0_1 (colonnes) (= W1, Tellier 2011) Private polymorphism in species 1

Sx1_0_y=$(echo $sum_first_column+$sum_last_column-$firstLINE_firstCOL-$lastLINE_firstCOL-$firstLINE_lastCOL-$lastLINE_lastCOL | bc)
echo $Sx1_0_y

#############POLYMORPHISME FIXE SP Sf_0_1# (=W3, Tellier 2011) fixed differences between sp

Sf_0_y=$(echo $lastLINE_firstCOL+$firstLINE_lastCOL | bc)
echo $Sf_0_y

#############POLYMORPHISME PARTAGE ENTRE LES DEUX SP Ss_0_1 (colonne) (=W4, Tellier 2011)
sumW1_W2_W3=$(echo $firstLINE_firstCOL+$Sx0_0_y+$Sx1_0_y+$Sf_0_y+$lastLINE_lastCOL| bc)
Ss_0_y=$(echo $ALL_table-$sumW1_W2_W3 | bc)
echo $Ss_0_y

###########HEADER pour mes noms de stats et création du fichier avec les quatre summ stats calculées

hdpoly_0_y=$(echo -e "Sx1_$i\tSx0_$i\tSs_$i\tSf_$i\t")
poly_0_y=$(echo -e "${Sx1_0_y}\t${Sx0_0_y}\t${Ss_0_y}\t${Sf_0_y}\t")

##########summary_stats_temp1.txt = à la fin sera le fichier avec toutes les SumStats polymorphisme shared, fixed et Sx

hdssstat_bef=$(head -n1 'summary_stats_temp1.txt')
ssstat_bef=$(tail -n +2 summary_stats_temp1.txt)

#########concatenation des SumStats Sx Sf et Ss calcules dans cette loop, à celles de la loop précédente

echo -e "${hdssstat_bef}${hdpoly_0_y}" > summary_stats_temp1.txt
echo -e "${ssstat_bef}${poly_0_y}" >> summary_stats_temp1.txt

rm JFS-temp_3.txt
rm JFS-temp_2.txt

done

#####READAPTATION DU FICHIER AVEC TOUTES LES STATS DE POLYMORPHISME SANS LE BUG (i.e. ajout de 0 pour la première simu test faite par ABCtoolBox)

tail -n +2 summary_stats_temp1.txt > summary_stats_temp2.txt

###PREPRATION SU HEADER DES SUMM STATS

head -n1 'summary_stats_temp1.txt' > summary_stats_temp3.txt

#je ne prends que les colonnes de 41 à la fin (car bug: les noms des SumStats sont répétées deux fois)
#awk '{for (i=41; i<=NF; i++) printf "%s ", $i; printf "\n"; }' summary_stats_temp4.txt > summary_stats_temp5.txt

### HEADER ET SUMM STATS CORRIGEES ENSEMBLE dans summary_stats_temp6.txt

head -n1 'summary_stats_temp3.txt' > summary_stats_temp6.txt
head -n1 'summary_stats_temp2.txt' >> summary_stats_temp6.txt

