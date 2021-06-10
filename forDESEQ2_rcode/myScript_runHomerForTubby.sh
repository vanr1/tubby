# (source) http://homer.ucsd.edu/homer/microarray/index.html


### COPY AND PASTE IN THESE COMMANDS TO TERMINAL
# (default run)
findMotifs.pl tubby_inputForHomer_allDEGs.txt mouse tubbyOutput_defaultRun_allDEGs/ -start -1000 -end 100 -len 8,10 -p 4
findMotifs.pl tubby_inputForHomer_upDEGs.txt mouse tubbyOutput_defaultRun_upDEGs/ -start -1000 -end 100 -len 8,10 -p 4
findMotifs.pl tubby_inputForHomer_dnDEGs.txt mouse tubbyOutput_defaultRun_dnDEGs/ -start -1000 -end 100 -len 8,10 -p 4
pwd




# (run that looks for ERRE sequence; notice -f flag)


seq2profile.pl TCAAGGTCA 1 erre > erre.motif
findMotifs.pl tubby_inputForHomer.txt mouse tubbyOutput_findERRE_noMismatch/ -start -1000 -end 100 -len 8,10 -p 4 -find erre.motif > outputfile_genesWithMotifFor_erre.txt

wc -l outputfile_genesWithMotifFor_erre.txt && head -1 outputfile_genesWithMotifFor_erre.txt | awk -F"\t" '{print NF; exit}'