#how to run the snp mutator
for x in `seq 1 22` X Y
do
perl ./SNP_mutator_v2.pl -f /data2/bsi/reference/sequence/human/ncbi/37.1/chr${x}.fa -o chr${x}.random.snp.indels -m 5000 -i 1000 -d 1000

echo "Done with $x"
done
#how to generate the sv's
sh ./SV_generator.sh /data4/bsi/epibreast/m087494.couch/Development/Genome_Smasher/Scripts/ /data2/bsi/reference/sequence/human/ncbi/37.1/ All_out
