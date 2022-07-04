#overlap >15bp on same strand = annotated
#everything else = noncoding
#no overlap at all = intergenic
#noncoding - intergenic = antisense
#any noncoding within 100bp of an annotated element = readthrough

#get all four text files first (annotated,antisense,intergenic,readthrough)

cut -f 2 -d "\"" REL606_400_commonwindows.gtf | sort -u | sed "s/^/\"/g" | sed "s/$/\"/g" > all_commonwindows_400.txt
../../../../../RNAseq_denovo/tools/bedtools intersect -a REL606_400_commonwindows.gtf -b ../../../converted_gtf_files/REL606.gtf -s -wo | awk -F '\t' '($19>15)' | cut -f 2 -d "\"" | sed "s/^/\"/g" | sed "s/$/\"/g" | sort -u > annot_commonwindows_400.txt
grep -v -f annot_commonwindows_400.txt all_commonwindows_400.txt > noncoding_commonwindows_400.txt
../../../../../RNAseq_denovo/tools/bedtools intersect -a REL606_400_commonwindows.gtf -b ../../../converted_gtf_files/REL606.gtf -v | cut -f 2 -d "\"" | sed "s/^/\"/g" | sed "s/$/\"/g" | sort -u > intergenic_commonwindows_400_step1.txt
grep -v -f intergenic_commonwindows_400_step1.txt noncoding_commonwindows_400.txt > antisense_commonwindows_400_step1.txt
grep -f antisense_commonwindows_400_step1.txt REL606_400_commonwindows.gtf | sed "s/\t/%/g" | sort -t '%' -nk 4 | sed "s/%/\t/g" > antisense_commonwindows_400_step1.gtf
../../../../../RNAseq_denovo/tools/bedtools intersect -a antisense_commonwindows_400.gtf -b ../../../converted_gtf_files/REL606.gtf -v -S | cut -f 2 -d "\"" | sed "s/^/\"/g" | sed "s/$/\"/g" > interim
grep -v -f interim antisense_commonwindows_400_step1.txt > antisense_commonwindows_400.txt
cat interim intergenic_commonwindows_400_step1.txt > intergenic_commonwindows_400.txt
rm *step1.txt
grep -f noncoding_commonwindows_400.txt REL606_400_commonwindows.gtf | sed "s/\t/%/g" | sort -t '%' -nk 4 | sed "s/%/\t/g" > noncoding_commonwindows_400.gtf
../../../../../RNAseq_denovo/tools/bedtools closest -a noncoding_commonwindows_400.gtf -b ../../../converted_gtf_files/REL606.gtf -s -d | awk -F '\t' '($19<100)' | cut -f 2 -d "\"" | sort -u > readthrough_commonwindows_400.txt
sed -i "s/^/,/g" readthrough_commonwindows_400.txt 
sed -i "s/$/,/g" readthrough_commonwindows_400.txt 
for file in *commonwindows_400*txt; do sed -i "s/\"/,/g" $file; done
sed -i "s/,,/,/g" readthrough_commonwindows_400.txt

#Mark the files

echo 'grep -f annot_commonwindows_400.txt htseq_time_series_400bp_namesorted_foldchanges.csv | sed "s/$/,annotated,annotated/g" > htseq_time_series_400bp_namesorted_foldchanges_mark1.csv
grep -f antisense_commonwindows_400.txt htseq_time_series_400bp_namesorted_foldchanges.csv | grep -f readthrough_commonwindows_400.txt | sed "s/$/,antisense,readthrough/g" >> htseq_time_series_400bp_namesorted_foldchanges_mark1.csv
grep -f antisense_commonwindows_400.txt htseq_time_series_400bp_namesorted_foldchanges.csv | grep -v -f readthrough_commonwindows_400.txt | sed "s/$/,antisense,antisense/g" >> htseq_time_series_400bp_namesorted_foldchanges_mark1.csv
grep -f intergenic_commonwindows_400.txt htseq_time_series_400bp_namesorted_foldchanges.csv | grep -f readthrough_commonwindows_400.txt | sed "s/$/,intergenic,readthrough/g" >> htseq_time_series_400bp_namesorted_foldchanges_mark1.csv
grep -f intergenic_commonwindows_400.txt htseq_time_series_400bp_namesorted_foldchanges.csv | grep -v -f readthrough_commonwindows_400.txt | sed "s/$/,intergenic,intergenic/g" >> htseq_time_series_400bp_namesorted_foldchanges_mark1.csv
sed -i "1 i\generation,strain,target_id,baseMean,log2FoldChange,lfcSE,pvalue,padj,type,rt_status" htseq_time_series_400bp_namesorted_foldchanges_mark1.csv' > ts_fc_code.sh
sed "s/time_series/50000gen/g" ts_fc_code.sh | sed "s/_namesorted//g" | sed "s/generation,strain,/seqtype,line,/g" > 50000gen_fc_code.sh
sed "s/time_series/50000gen/g" ts_fc_code.sh | sed "s/_namesorted//g" | sed "s/foldchanges/count/g" | sed "s/generation,strain,target_id,baseMean,log2FoldChange,lfcSE,pvalue,padj/replicate,seqtype,line,target_id,est_counts/g" > 50000gen_counts_code.sh
sed "s/foldchanges/count/g" ts_fc_code.sh | sed "s/generation,strain,target_id,baseMean,log2FoldChange,lfcSE,pvalue,padj/generation,strain,replicate,target_id,est_counts/g" > ts_count_code.sh
cat *.sh > marking_code.sh
bash marking_code.sh

for file in *mark1.csv; do sed -i "s/AraM/Ara-/g" $file; done
for file in *mark1.csv; do sed -i "s/AraP/Ara+/g" $file; done
