




# making a variable for our main pipeline path
TUBBY_PIPELINE=/data2/han_lab/richardvan/tubby_pipeline
INPUT=output/retina_3mo/COUNTS_normalized_retina_3mo_DATA.txt

set -e  #exits the script when an error occurs for safety


declare -a candidateGenes=( 



		# 							# (b/c) 5 of 13 DEGs of hypothalamus 24d
		# "ENSMUSG00000073643"		# WD repeat and FYVE domain; Wdfy1	(-2.2 logFC)			
		# "ENSMUSG00000000560"		# Gama-aminobutyric acid type A receptor; Gabra2 (-1.9 logFC)
		# "ENSMUSG00000081470"		# NA (-4.44 logFC)
		# "ENSMUSG00000096768"		# NA (8.57 logFC)
		# "ENSMUSG00000093954"		# Mitochondrial inner membrane; Gm16867 (9.6 logFC)


		# 							# (b/c) 4 DEGs of hypothalamus 3mo
		# "ENSMUSG00000030761"		# Myosin VIIA; Myo7a	(2.1 logFC)			
		# "ENSMUSG00000072844"		# NA (3.5 logFC)
		# "ENSMUSG00000098975"		# NA (6.5 logFC)
		# "ENSMUSG00000093954"		# Mitoferrin-1;Gm16867 (8.8 logFC)



		# 							# (b/c) 6 of 47 DEGs of retina 24d
		# 							# Myosin VIIA; Myo7a AGAIN	(4.1 logFC)	
		# "ENSMUSG00000096718"		# Krueppel-associated box; Zfp781 (-21 logFC)				
		# "ENSMUSG00000028635"		# Endothelin 2 toxin; Edn2 (-4.23 logFC)
		# "ENSMUSG00000055368"		# Solute carrier family 6 noradrenalin (-3.2 logFC)
		# "ENSMUSG00000093954"		# Mitoferrin-1;Gm16867 (9.5 logFC)l AGAIN
		# "ENSMUSG00000049908"		# Gap junction protein connexin; Gja8 (3.2 logFC)
		# "ENSMUSG00000027360"		# Histidine decarboxylase; Hdc (2.5 logFC)


		# 							# (b/c) 6 of 5028 DEGs of retina 3mo
		# "ENSMUSG00000074738"		# Fibronectin type 3; Fndc10 (-10.4 logFC)				
		# "ENSMUSG00000098188"		# Sosondowah ankyrin; Soahc (-9 logFC)
		# "ENSMUSG00000024912"		# Fos-like antigen, AP-1 TF ; Fosl1 (-6.77 logFC)
		# "ENSMUSG00000031966"		# Galactosidase; Glb1l3 (4.5 logFC)
		# "ENSMUSG00000030470"		# Zinc finger, Cys and Gly rich; Csrp3 (3.5 logFC)
		# "ENSMUSG00000043020"		# Dynein axonemal intermediate chain; DnaI (3.1 logFC)

		# # 							# (b/c) latest curiosity
  #       "ENSMUSG00000064177"     	# Ghrelin; Ghrl
  #       "ENSMUSG00000057722"     	# Leptin receptor; Lepr
  #       "ENSMUSG00000021255"     	# Estrogen related receptor beta; Esrrb
  #       "ENSMUSG00000095562"     	# Eryhtorid differentiation regulatorErdr1
  #       "ENSMUSG00000036437"     	# Neuropeptide Y receptor Y1; Npy1r
  #       "ENSMUSG00000020660"     	# Pro-Opiomelanocortin-alpha; Pomc



		# 							# (b/c) of early "sanity checks"
  #       "ENSMUSG00000030324"     	# rhodopsin
  #       "ENSMUSG00000058831"       	# opsin 
  #       "ENSMUSG00000005705"    	# agouti-related neuropeptide
  #       "ENSMUSG00000037727"		# arginine vasopressin
  #       "ENSMUSG00000031028"     	# tub


									# (b/c) these are housekeeping genes; recommended by AMB
        "ENSMUSG00000029580"     	# actin
        "ENSMUSG00000119584"		# 18S rRNA; 
        "ENSMUSG00000072235"     	# tubulin; tuba1a (1 of 20 genes)
        "ENSMUSG00000057666"  		# GADPH
        ) 

                  
# clear the tmp folder
rm -rf tmp/*



# list all samples we will be looking at, make sure all their paired-end read files exist
echo ""
echo "== Listing all candidateGenes we will be looking at ..."
for candidateGene in ${candidateGenes[@]}
do
	echo $candidateGene

done

echo "== Now processing each one ..."
for candidateGene in ${candidateGenes[@]}
do
	
	cat <(head -1 $INPUT) <(grep $candidateGene $INPUT) | awk '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
	    }
	}
	NF>p { p = NF }
	END {    
	    for(j=1; j<=p; j++) {
	        str=a[1,j]
	        for(i=2; i<=NR; i++){
	            str=str"\t"a[i,j];
	        }
	        print str
	    }
	}' > tmp/EXPRESSIONCHECKER_${candidateGene}.txt

done


echo "== Post-Processing checking got files expected before merge ..."
for candidateGene in ${candidateGenes[@]}
do
	head -4 tmp/EXPRESSIONCHECKER_${candidateGene}.txt
	cut -f2 tmp/EXPRESSIONCHECKER_${candidateGene}.txt > tmp/2NDCOLUMN_${candidateGene}.txt 

	cut -f1 tmp/EXPRESSIONCHECKER_${candidateGene}.txt > tmp/1STCOLUMN.txt  # inefficient, but quick code-wise
done


echo "== Merging tmp files..."
for candidateGene in ${candidateGenes[@]}
do
	paste tmp/2NDCOLUMN_*.txt > tmp/EXPRESSION_merged.txt
done
paste tmp/1STCOLUMN.txt tmp/EXPRESSION_merged.txt > tmp/EXPRESSION_merged2.txt


echo "== Transposing in bash for final output in format I like ..."

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' tmp/EXPRESSION_merged2.txt > tmp/EXPRESSION_tranposed.txt

cp -av tmp/EXPRESSION_tranposed.txt output/EXPRESSION_candidateGenes_NormalizedCounts.txt

echo "===== DONE; check the final file now ..."

