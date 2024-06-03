#!/bin/sh
# this script adds rsID column to OLINK summary statistics using hg19 genome build as reference

# go to working directory where hg19 genome build is stored (i.e., "Kaviar-160204-Public-hg19-trim.vcf.gz")
cd /Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/OLINK/Kaviar-160204-Public/vcfs || exit

# define paths
SUMSTATS=/Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/OLINK/
echo "$SUMSTATS"
BUILD=/Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/OLINK/Kaviar-160204-Public/vcfs
echo "$BUILD"

# look at the reference build that contains rsids 
gunzip -cd Kaviar-160204-Public-hg19-trim.vcf.gz | head -n 100000

# look at the GWAS sumstats
head ${SUMSTATS}VEGF.A-1.tbl

# match GWAS sumstats CHR and position with that in the genome build hg19 
gunzip -cd Kaviar-160204-Public-hg19-trim.vcf.gz |  head -n 10000000 | grep "55158855" 
gunzip -cd Kaviar-160204-Public-hg19-trim.vcf.gz |  head -n 10000000 | grep "71285495"
# turns out its the hg19 build (not hg38) as it matches the information in VEGF.A-1.tbl file 
# (same chr and position in hg38 are not present) 

# unzip the hg19 build file and then use the unzipped file 
gunzip Kaviar-160204-Public-hg19-trim.vcf.gz

## LOOP THROUGH OLINK GWAS SUMSTATS 
# (note, for 91 olink markers the loop may take several days, so might want to split it into several jobs)
# the loop will add rsID based on position and chromosome as last column to GWAS sumstats
# unzip the GWAS sumstats
find ${SUMSTATS} -name "*.tbl.gz" -exec gunzip {} \; 
# same as above for remaining traits in a loop (can also do it in sections as separate jobs)
# added OFS='\t' to make sure even lines with "." are tab delimited in output 
# replaced "." entries for missing rsIDs in hg19 build with MarkerName (e.g., "chr1:23539383_A_C")
cd "$SUMSTATS" || exit
for FILE in *.tbl.gz
do 
  TRAIT=$(basename "$FILE" .gz)
  echo "$TRAIT"
  gunzip "$TRAIT".gz
  cd "$BUILD" || exit
  awk 'NR==FNR{a[$1,$2]=$3;next} NR!=FNR{if (($1,$2) in a) printf("%s\t%s\n", $0, a[$1,$2]); else if (FNR==1) printf("%s\trsID\n", $0); }' Kaviar-160204-Public-hg19-trim.vcf ${SUMSTATS}"$TRAIT" | gzip > "$TRAIT".rsid.txt.gz
  gunzip -cd "$TRAIT".rsid.txt.gz | awk 'NR==1 {print $0} NR!=1 {if ($NF == ".") $NF = $3; print}' OFS='\t' | gzip > tmp && mv tmp "$TRAIT".rsid.txt.gz
  cd "$SUMSTATS" || exit
done
# I doubled checked correct variants have been found using http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl 

# create directory for completed files
mkdir completed
mv *.tbl completed/ 
# can compress these to take up less space
tar -czvf completed.tar.gz completed # super slow
zip -r completed.zip completed # super slow

