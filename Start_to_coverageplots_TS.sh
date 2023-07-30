#Get all mutation IDs for the TS dataset
#Put all genomes of these clones in one file too

for i in ZDB409 ZDB425 ZDB445 ZDB467 ZDB483 ZDB486 ZDB488 ZDB309 ZDB17 ZDB199 ZDB200 ZDB25 ZDB564 CZB154; do grep -v "#" ../*"$i".gd | awk -F '\t' '{OFS=""}{print $1,"_",$5,"_",$6}'; done | sort -u > time_series_mutation_ids
for i in $(ls ../../../all_genomes/Ara-3_*fasta | rev | cut -f 2- -d '.' | rev | cut -f 5 -d '/'); do sed "s/>REL606/>"$i"/g" ../../../all_genomes/"$i"*fasta >> Ara-3_genomes.faa; done

#for i in ZDB409 ZDB425 ZDB445 ZDB467 ZDB483 ZDB486 ZDB488 ZDB309 ZDB17 ZDB199 ZDB200 ZDB25 ZDB564 CZB154; do grep -v "#" ../*"$i".gd | sed "s/^/"$i"\t/g"; done | awk -F '\t' '{OFS=""}{print $1,"_",$2,"_",$6,"_",$7,"_",$8,"_",$9,"_",$10,"_",$11,"_",$12,"_",$13}' | sed "s/_*$//g" > time_series_mutations_marked

#The only thing these mutations share is their position relative to the ancestor
#Get their coordinates in the ancestor first

egrep "SNP|INS|MOB" ../time_series_mutation_ids | awk -F '_' '{OFS=""}{print "REL606\t.\tCDS\t",$2,"\t",$2+200,"\t.\t+\t0\ttranscript_id \"",$0,"_plus_down\";gene_id \"",$0,"_plus_down\";"}' > time_series_mutation_ids_ancestor.gtf
egrep "SNP|INS|MOB" ../time_series_mutation_ids | awk -F '_' '{OFS=""}{print "REL606\t.\tCDS\t",$2-200,"\t",$2,"\t.\t-\t0\ttranscript_id \"",$0,"_minus_up\";gene_id \"",$0,"_minus_up\";"}' >> time_series_mutation_ids_ancestor.gtf
egrep "DEL|SUB|CON|AMP" ../time_series_mutation_ids | awk -F '_' '{OFS=""}{print "REL606\t.\tCDS\t",$2+$3,"\t",$2+$3+200,"\t.\t+\t0\ttranscript_id \"",$0,"_plus_down\";gene_id \"",$0,"_plus_down\";"}' >> time_series_mutation_ids_ancestor.gtf
egrep "DEL|SUB|CON|AMP" ../time_series_mutation_ids | awk -F '_' '{OFS=""}{print "REL606\t.\tCDS\t",$2-200,"\t",$2,"\t.\t-\t0\ttranscript_id \"",$0,"_minus_up\";gene_id \"",$0,"_minus_up\";"}' >> time_series_mutation_ids_ancestor.gtf
sort -nk 4 -o time_series_mutation_ids_ancestor.gtf time_series_mutation_ids_ancestor.gtf

#Only retain non-genic ones

bedtools intersect -s -a time_series_mutation_ids_ancestor.gtf -b ../../../../converted_gtf_files/REL606.gtf -wo | awk -F '\t' '(($12=="CDS"||$12=="gene"||$12=="repeat_region"||$12~"RNA"||$12~"fCDS")&&($19>10))' | cut -f 9 | sort -u > genic_elements
grep -v -f genic_elements time_series_mutation_ids_ancestor.gtf | sort -nk4 > time_series_mutation_ids_ancestor_nongenic.gtf

#blast against ancestral genome (make blastdb out of REL606), parse to get rid of duplicates

grep -v -f genic_elements time_series_mutation_ids_ancestor.gtf | sort -nk4 | gtf2bed | bedtools getfasta -name -fi ../../../../all_genomes/REL606.fasta -bed - > time_series_mutation_ids_ancestor_nongenic.faa
blastn -query time_series_mutation_ids_ancestor_nongenic.faa -db ../REL606 -num_threads 72 -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -max_target_seqs 500 -max_hsps 500 -out ts_mutation_vs_REL606.tsv
awk -F '\t' '($3>80&&$5>80)' ts_mutation_vs_REL606.tsv | cut -f 1 | sort | uniq -c | awk '($1==1)' | rev | cut -f 1 -d ' ' | rev > nondup_mutation_ids.txt
grep -f nondup_mutation_ids.txt time_series_mutation_ids_ancestor.gtf > nondup_time_series_mutation_ids_ancestor.gtf
grep --no-group-separator -A1 -f nondup_mutation_ids.txt time_series_mutation_ids_ancestor_nongenic.faa > nondup_time_series_mutation_ids_ancestor.faa

#Map them in evolved clones to get coordinates. For some reason gmap doesn't work unless the genome indices are found in the same directory, so cd.. before running

for i in $(ls ../../../all_genomes/Ara-3_*fasta | rev | cut -f 2- -d '.' | rev | cut -f 5 -d '/'); do /stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 200bp_timeseries/nondup_time_series_mutation_ids_ancestor.faa > 200bp_timeseries/"$i"_timeseries_nongenic.gff3; done

#back to working directory

for i in $(ls *nongenic*gff3 | cut -f 1-3 -d "_"); do awk -F '\t' '($3=="mRNA")' "$i"*gff3 | grep "mrna1" | cut -f 1 -d ';' | sed "s/ID=/transcript_id \"/g" | sed "s/\.mrna1/\";/g" | sed "s/mRNA/CDS/g" | awk -F '\t' '{OFS=FS}{ if ($0 ~ /minus/) $7 = "-"; print }' | awk -F '\t' '{OFS=FS}{print $0,$9}' | sed "s/;\ttranscript/;gene/g" > "$i"_nongenic.gtf; done

Prepping for htseq:
cat ../../../../insremoved_gtf/Ara-3*insremoved.gtf ../../../../insremoved_gtf/REL606_insremoved.gtf | cut -f 2 -d "\"" | sort | uniq -c | grep " 16 " | rev | cut -f 1 -d ' ' | rev | sort -u > common_annotated.txt
for i in $(ls *Ara*nongenic*gtf | cut -f 1-3 -d "_"); do grep -f common_annotated.txt ../../../../insremoved_gtf/"$i"*gtf > "$i"_htseqready.gtf && cat "$i"*nongenic*gtf >> "$i"_htseqready.gtf; done
grep -f common_annotated.txt ../../../../insremoved_gtf/REL606_insremoved.gtf > REL606_htseqready.gtf && cat nondup_time_series_mutation_ids_ancestor.gtf >> REL606_htseqready.gtf

#htseq:
for i in $(ls Ara-3*htseqready.gtf | cut -f 1-3 -d "_"); do ls ../../../../../trimmomatic_29bp_trimmed/bamfiles/namesorted/*sorted.bam | grep "$i" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_htseqready.gtf | head -n -5 > T",$11}' | sed "s/T..\/..\/..\/..\/..\/trimmomatic_29bp_trimmed\/bamfiles\/namesorted\///g" | sed "s/$/t/g" | sed "s/_namesorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_htseq_raw.tsv"}'; done > timeseries_htseq.sh
ls ../../../../../trimmomatic_29bp_trimmed/bamfiles/namesorted/*sorted.bam | egrep "REL60" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk '{OFS=""}{print $0," REL606_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T\.\.\/\.\.\/\.\.\/\.\.\/..\/trimmomatic_29bp_trimmed\/bamfiles\/namesorted\///g" | sed "s/$/t/g" | sed "s/_namesorted.bamt//g" | awk '{OFS=""}{print $0,"_"}' | sed "s/$/_htseq_raw.tsv/g" >> timeseries_htseq.sh
ls ../../../../../50k_bowtie_indices_bamfiles/*sorted.bam | grep "Ara-3" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk '{OFS=""}{print $0," Ara-3_50000gen_11364_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk '{OFS=""}{print $0,"_htseq_raw.tsv"}' > 50k_ts_htseq.sh
ls ../../../../../50k_bowtie_indices_bamfiles/*sorted.bam | grep "REL60" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk '{OFS=""}{print $0," REL606_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk '{OFS=""}{print $0,"_htseq_raw.tsv"}' >> 50k_ts_htseq.sh

#Now add in deleted genes so the number of genes are equal across all, ascribe 0 counts to them

for i in $(ls Ara-3*htseqready.gtf | cut -f 1-3 -d "_")
do
cut -f 2 -d "\"" "$i"_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> "$i"_BiolRep1*raw.tsv
cut -f 2 -d "\"" "$i"_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> "$i"_BiolRep2*raw.tsv
cut -f 2 -d "\"" "$i"_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> "$i"_BiolRep3*raw.tsv
done

cut -f 2 -d "\"" Ara-3_50000*_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> rep1-ribo-Ara-3_htseq_raw.tsv
cut -f 2 -d "\"" Ara-3_50000*_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> rep2-ribo-Ara-3_htseq_raw.tsv
cut -f 2 -d "\"" Ara-3_50000*_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> rep1-rna-Ara-3_htseq_raw.tsv
cut -f 2 -d "\"" Ara-3_50000*_htseqready.gtf | grep -v "ECB" > temp && grep -v -f temp nondup_time_series_mutation_ids_ancestor.gtf | cut -f 2 -d "\"" | sed "s/$/\t0/g" >> rep2-rna-Ara-3_htseq_raw.tsv

#Add in line, replicate and generation information

ls Ara-3*BiolRep*tsv | cut -f 1-4 -d "_" | awk -F '_' '{OFS=""}{print "sed -i \"s\/^\/",$1,"\\t",$2,"\\t",$3,"\\t",$4,"\\t\/g\""}' | paste -d '\t' - <(ls Ara-3*BiolRep*tsv) | awk '{print $1,$2,$3,$4}' | bash
ls REL606*BiolRep*tsv | cut -f 1-4 -d "_" | awk -F '_' '{OFS=""}{print "sed -i \"s\/^\/\REL606\\t0gen\\t",$1,"\\t",$2,"\\t\/g\""}' | paste -d '\t' - <(ls REL606*BiolRep*tsv) | awk '{print $1,$2,$3,$4}' | bash

#Run DeSeq2 on the TS dataset!

###

#Get all increases (candidate proto-genes)

cat 100bp_timeseries/TIMESERIES_deseq2_results.csv 200bp_timeseries/TIMESERIES_deseq2_results.csv | awk -F ',' '($8<0.05&&$5>0&&$3!~"ECB")' | cut -f 3 -d ',' | sort -u > 12_increases.txt

