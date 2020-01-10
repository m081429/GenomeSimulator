#!/bin/sh
EXPECTED=3
if [ $# -ne $EXPECTED ]
then
	echo "Usage: SV_Generator.sh <SCRIPTS_DIR> <Genome_DIR> <OUTPUT_VCF> <NUM_BASES_FASTAFILE_EACHLINE>"
	exit
fi
OUTPUT_VCF=All_out
SCRIPTS_DIR=$1
GENOME_DIR=$2
OUTPUT_VCF=$3
FINAL_VCF=$OUTPUT_VCF

ls "$GENOME_DIR/chr"[0-9]*.fa | perl -MList::Util=shuffle -e 'print shuffle(<STDIN>);'> reference_files

for i in {1..4}
do
CHROMOSOME=$(awk -v i="$i" '(NR==i)' reference_files)
x=${CHROMOSOME##*/}
x=${x%%.fa}
#echo $CHROMOSOME
perl $SCRIPTS_DIR/Novel_insertions.pl $CHROMOSOME $OUTPUT_VCF 10
done
for i in {5..8}
do
CHROMOSOME=$(awk -v i="$i" '(NR==i)' reference_files)
perl $SCRIPTS_DIR/Large_Deletions.pl $CHROMOSOME $OUTPUT_VCF 20
done

for i in {9..12}
do
CHROMOSOME=$(awk -v i="$i" '(NR==i)' reference_files)
perl $SCRIPTS_DIR/Copy_number_changes.pl $CHROMOSOME $OUTPUT_VCF 10
done

for i in {16..18}
do
CHROMOSOME=$(awk -v i="$i" '(NR==i)' reference_files)
perl $SCRIPTS_DIR/Inversions.pl $CHROMOSOME $OUTPUT_VCF 5
done	

for i in 19 
do
CHROMOSOME1=$(awk -v i="$i" '(NR==i)' reference_files)
j=$(($i+2))
CHROMOSOME2=$(awk -v i=$j '(NR==i)' reference_files)
perl $SCRIPTS_DIR/Translocations.pl $CHROMOSOME1 $CHROMOSOME2 $OUTPUT_VCF 2
done

rm reference_files
