#Consolidate the upregulated cases, duplicates and redundancies removed, from both datasets

cat 63_coverage_candidates.txt timeseries/12_increases_withstrand.txt | egrep -v "INS_2103918_CCAGCCAG|CON_3549957_4|3015771|4415710" > 70_bothdataset_upregs.txt

#Make a .bed file for the cases

cat ancestor_upregs.bed 12_increases.bed | awk -F '\t' '{OFS=FS}{$1="REL606"}1' | grep -f 70_bothdataset_upregs.txt > 70_bothdataset_upregs.bed

#make a gtf file with both +100 and +200 bp regions, so 140 total

cat timeseries/100bp_timeseries/time_series_mutation_ids_ancestor.gtf 100bp/*ancestor.gtf | grep -f 70_bothdataset_upregs.txt > 70_100bp.gtf
cat timeseries/200bp_timeseries/time_series_mutation_ids_ancestor.gtf 200bp/*ancestor.gtf | grep -f 70_bothdataset_upregs.txt >> 70_100bp.gtf

#Run htseq using prior MURI code with 70_100bp.gtf

#All unambiguous proto-genes (n=28) that have more than 3 reads in either their 100 or 200bp downstream regions

cat test* | grep -f /stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/unambiguous_protogenes.txt | awk -F '\t' '($2>3)' | cut -f1 | cut -f 2 -d ',' | sed "s/_100$//g" | sed "s/_200$//g" | sort -u > MURI_coverage_candidates_step1.txt
cat test* | grep -f MURI_coverage_candidates_step1.txt | awk -F '\t' '($2>50)' | cut -f1 | cut -f 2 -d ',' | sed "s/_100//g" | sed "s/_200//g" | sort -u > over50_reads.txt
grep -v -f over50_reads.txt MURI_coverage_candidates_step1.txt > MURI_coverage_candidates_step2.txt
cat test* | grep -f MURI_coverage_candidates_step2.txt | awk -F '\t' '($2>3)' | cut -f 1 -d ',' | sort -u | cut -f 1 -d '.' > 41_MURI_coveragecases.txt

#Make genomecov files out of the 41 MURI files

#Since they're paired-end, make bam files just with R1

for i in $(cat byname/41*); do bowtie2 -x /stor/work/Ochman/hassan/LTEE_analysis/trimmomatic_29bp_trimmed/indices/REL606.fasta -U fastqs/"$i"_1_paired.fastq --very-sensitive-local -p 104 -S $i_R1.sam; done
for i in $(ls *sam | cut -f 1 -d '.'); do samtools view -b -@ 72 -S "$i".sam > "$i".bam; done
for i in $(ls *sorted.bam | rev | cut -f 3- -d "_" | rev);

#Genomecov for these 41 files

do
bedtools genomecov -d -strand + -ibam "$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t+/g" > "$i"_plus
bedtools genomecov -d -strand - -ibam "$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t-/g" > "$i"_minus
cat "$i"_plus "$i"_minus > "$i"_coverage.tsv
done



X=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 3); Y=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' *_coverage_.tsv >> something_right
X=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 2); Y=$(grep "something" ../70_bothdataset_upregs.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' *_coverage_.tsv >> something_left

#replace something with name of the genes

