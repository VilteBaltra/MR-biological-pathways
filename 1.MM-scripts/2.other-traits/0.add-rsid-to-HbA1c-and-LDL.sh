#!/bin/sh
# this script adds rsID column to OLINK summary statistics using hg19 genome build as reference

# go to working directory where hg19 genome build is stored (i.e., "Kaviar-160204-Public-hg19-trim.vcf.gz")
cd /Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/OLINK/Kaviar-160204-Public/vcfs || exit

# define paths
SUMSTATS=/Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/sumstats-2024/
echo "$SUMSTATS"
BUILD=/Users/vb506/Documents/Projects/MR-mediation-CM-MM/summary-stats/OLINK/Kaviar-160204-Public/vcfs
echo "$BUILD"

### LDL ###
# look at the reference build that contains rsids 
head -n 100 Kaviar-160204-Public-hg19-trim.vcf

# look at the GWAS sumstats
# for LDL
gunzip -cd ${SUMSTATS}GCST90018961_buildGRCh37.tsv.gz | grep "55158855"
gunzip -cd ${SUMSTATS}GCST90018961_buildGRCh37.tsv.gz | grep "71285495"

# match GWAS sumstats CHR and position with that in the genome build hg19 
head -n 10000000 Kaviar-160204-Public-hg19-trim.vcf | grep "55158855" 
head -n 10000000 Kaviar-160204-Public-hg19-trim.vcf | grep "71285495"
# doubled check and it's hg19 build (as also indicated by file name 'GRCh37') 
# (same chr and position in hg38 are not present) 

# unzip the hg19 build file and then use the unzipped file 
gunzip Kaviar-160204-Public-hg19-trim.vcf.gz

# Add rsID based on position and chromosome as last column to GWAS sumstats
# unzip the GWAS sumstats
cd "$SUMSTATS" || exit
TRAIT=GCST90018961_buildGRCh37.tsv
echo "$TRAIT"
gunzip "$TRAIT".gz
cd "$BUILD" || exit
awk 'NR==FNR{a[$1,$2]=$3;next} NR!=FNR{if (($1,$2) in a) printf("%s\t%s\n", $0, a[$1,$2]); else if (FNR==1) printf("%s\trsID\n", $0); }' Kaviar-160204-Public-hg19-trim.vcf ${SUMSTATS}"$TRAIT" | gzip > "$TRAIT".rsid.txt.gz
gunzip -cd "$TRAIT".rsid.txt.gz | awk 'NR==1 {print $0} NR!=1 {if ($NF == ".") $NF = $3; print}' OFS='\t' | gzip > tmp && mv tmp "$TRAIT".rsid.txt.gz
cd "$SUMSTATS" || exit
# added OFS='\t' to make sure even lines with "." are tab delimited in output 
# replaced "." entries for missing rsIDs in hg19 build with MarkerName (e.g., "chr1:23539383_A_C")

# inspect 
wc -l "$SUMSTATS"GCST90018961_buildGRCh37.tsv
# rows original: 20,521,448
gunzip -cd "$TRAIT".rsid.txt.gz | wc -l 
# rows new:  18,726,199 (some less likely due to absence in ref panel)

# but looks good 
# view top rows
head "$SUMSTATS"GCST90018961_buildGRCh37.tsv
gunzip -cd "$TRAIT".rsid.txt.gz | head
# view bottom rows
tail "$SUMSTATS"GCST90018961_buildGRCh37.tsv
gunzip -cd "$TRAIT".rsid.txt.gz | tail

mv GCST90018961_buildGRCh37.tsv.rsid.txt.gz $SUMSTATS


### HbA1c ###
# HbA1c is on GRCh38, so needs different reference panel)
# will use Kaviar-160204-Public-hg38-trim.vcf.gz

# double check positions match (looks good)
gunzip -cd ${SUMSTATS}GCST90014006_buildGRCh38.tsv.gz | grep "2628150"
head -n 10000000 Kaviar-160204-Public-hg38-trim.vcf | grep "2628150" 

# Add rsID based on position and chromosome as last column to GWAS sumstats
# unzip the GWAS sumstats
cd "$SUMSTATS" || exit
TRAIT=GCST90014006_buildGRCh38.tsv
echo "$TRAIT"
gunzip "$TRAIT".gz
cd "$BUILD" || exit
awk 'NR==FNR{a[$1,$2]=$3;next} NR!=FNR{if (($1,$2) in a) printf("%s\t%s\n", $0, a[$1,$2]); else if (FNR==1) printf("%s\trsID\n", $0); }' Kaviar-160204-Public-hg38-trim.vcf ${SUMSTATS}"$TRAIT" | gzip > "$TRAIT".rsid.txt.gz
gunzip -cd "$TRAIT".rsid.txt.gz | awk 'NR==1 {print $0} NR!=1 {if ($NF == ".") $NF = $3; print}' OFS='\t' | gzip > tmp && mv tmp "$TRAIT".rsid.txt.gz
cd "$SUMSTATS" || exit
# added OFS='\t' to make sure even lines with "." are tab delimited in output 
# replaced "." entries for missing rsIDs in hg19 build with MarkerName (e.g., "chr1:23539383_A_C")

# inspect 
wc -l "$SUMSTATS"GCST90014006_buildGRCh38.tsv
# rows original: 11,367,992
gunzip -cd "$TRAIT".rsid.txt.gz | wc -l 
# rows new:  11,320,166 (almost all present)

# looks good 
# view top rows
head "$SUMSTATS"GCST90014006_buildGRCh38.tsv
gunzip -cd "$TRAIT".rsid.txt.gz | head
# view bottom rows
tail "$SUMSTATS"GCST90014006_buildGRCh38.tsv
gunzip -cd "$TRAIT".rsid.txt.gz | tail

mv GCST90014006_buildGRCh38.tsv.rsid.txt.gz $SUMSTATS

# completed