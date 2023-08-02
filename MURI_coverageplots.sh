#Consolidate the upregulated cases, duplicates and redundancies removed, from both datasets

cat 63_coverage_candidates.txt timeseries/12_increases_withstrand.txt | egrep -v "INS_2103918_CCAGCCAG|CON_3549957_4|3015771|4415710" > 70_bothdataset_upregs.txt

#Make a .bed file for the cases

cat ancestor_upregs.bed 12_increases.bed | awk -F '\t' '{OFS=FS}{$1="REL606"}1' | grep -f 70_bothdataset_upregs.txt > 70_bothdataset_upregs.bed

#make a gtf file with both +100 and +200 bp regions, so 140 total

cat timeseries/100bp_timeseries/time_series_mutation_ids_ancestor.gtf 100bp/*ancestor.gtf | grep -f 70_bothdataset_upregs.txt > 70_100bp.gtf
cat timeseries/200bp_timeseries/time_series_mutation_ids_ancestor.gtf 200bp/*ancestor.gtf | grep -f 70_bothdataset_upregs.txt >> 70_100bp.gtf

#Run htseq using prior MURI code with 70_100bp.gtf

X=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 3); Y=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' *_coverage_.tsv >> something_right
X=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 2); Y=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' *_coverage_.tsv >> something_left

#replace something with name of the genes
