#Construct ancestral transcriptomes using .gd (not applied.gd) files, as the latter refer to evolved coordinates and are therefore useless.

#For 500:

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
grep -v "#" "$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || i == 5 || i == 6 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$6,"\t",$6+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_plus_down_500\";"}' | egrep -v "IS.*[0-9]{1,3}_-1_plus_down" > "$i"_evolved.gtf
grep -v "#" "$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || i == 5 || i == 6 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_minus_up_500\";"}' | egrep -v "IS.*[0-9]{1,3}_1_minus_up" >> "$i"_evolved.gtf
grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5,"\t",$5+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";"}' | egrep -v "IS.*[0-9]{1,3}_-1_plus_down" > "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";"}' | egrep -v "IS.*[0-9]{1,3}_1_minus_up" >> "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5+$6,"\t",$5+$6+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";"}' >> "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";"}' >> "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_ancestor.gtf -o "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_evolved.gtf -o "$i"_evolved.gtf
done

#manually check to make sure no bad records

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
bedtools intersect -s -v -a "$i"_evolved.gtf -b ../../converted_gtf_files/"$i"_50000gen*gtf -wo | cut -f 2 -d "\"" | awk -F '_' '{OFS=FS}{print $1,$2,$3,".*",$6,$7}' > "$i"_nongenic_windows.txt
grep -f "$i"_nongenic_windows.txt "$i"_evolved.gtf > "$i"_evolved_nongenic.gtf
grep -f "$i"_nongenic_windows.txt "$i"_ancestor.gtf > "$i"_ancestor_nongenic.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
cat ../../insremoved_gtf/"$i"_50000gen*insremoved.gtf ../../insremoved_gtf/REL606_insremoved.gtf | cut -f 2 -d "\"" | sort | uniq -c | grep " 2 " | rev | cut -f 1 -d ' ' | rev > "$i"_common_annotated.txt
grep -f "$i"_common_annotated.txt ../../insremoved_gtf/REL606*gtf > "$i"_ancestor_htseqready.gtf
grep -f "$i"_common_annotated.txt ../../insremoved_gtf/"$i"_50000gen*gtf > "$i"_evolved_temp.gtf
cat "$i"_ancestor_nongenic.gtf >> "$i"_ancestor_htseqready.gtf
cat "$i"_evolved_nongenic.gtf >> "$i"_evolved_temp.gtf #gene names need to match
paste "$i"_ancestor_htseqready.gtf "$i"_evolved_temp.gtf | awk -F '\t' '{OFS=FS}{print $10,$11,$12,$13,$14,$15,$16,$17,$9}' > "$i"_evolved_htseqready.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
ls ../../../50k_bowtie_indices_bamfiles/*sorted.bam | egrep "REL60" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_ancestor_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' > "$i"_htseq_500.sh
ls ../../../50k_bowtie_indices_bamfiles/*sorted.bam | grep "$i" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_evolved_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' >> "$i"_htseq_500.sh
done

ls rep* | awk '{OFS=""}{print $1,"_",$1}' | awk -F '_' '{OFS=""}{print "sed -i \"s\/^\/",$1,"\\t/g\" ",$3,"_",$4}' | bash
for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
cat *"$i" > "$i"_500_htseq.csv
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6
do
sed "s/-ribo-/\tribo\t/g" "$i"_500_htseq.csv | sed "s/-rna-/\trna\t/g" | sed "1s/^/repl\tseqtype\tline\ttarget_id\test_counts\n/g" > "$i"_500_htseq.tsv
done

#DESEQ2

library(tidyverse)
library(DESeq2)
library(matrixStats)
library(apeglm)

kdf2 <- read_tsv("Ara-6_500_htseq.tsv") %>% 
  select(repl, seqtype, line, target_id, est_counts) %>% 
  unite("sample", c(seqtype, line, repl), sep = "_")
line.names <- unique(str_extract(unique(kdf2$sample),"[a-z]{3,4}_Ara[-+][1-6]"))
line.names <- line.names[!is.na(line.names)]
names(line.names) <- line.names
wide.kdf2 <- kdf2 %>% 
  pivot_wider(names_from = sample, values_from = est_counts)
wide.kdf2 <- wide.kdf2 %>% rowwise %>% mutate(row.sum = sum(c_across(starts_with("rna")))) %>% filter(row.sum>5) %>% select(1:13)
df.list <- lapply(line.names, function(x){
  seqtype <- unlist(str_split(x, "_"))[1]
  anc.to.pick <- paste(seqtype, "REL", sep = "_")
  wide.kdf2 %>%
    select(target_id, starts_with(x), starts_with(anc.to.pick)) %>% # pick the correct columns, one evo and both anc
    mutate_if(is.double, as.integer) %>%                            # convert cols to integers
    column_to_rownames("target_id") %>%                             # gene ids to rownames
    as.matrix()                                                     # change to a matrix
})
(conds.df <- data.frame(condition = factor(c(rep("evo", 2), rep("anc", 4)), levels = c("anc", "evo"))))

deseq.list2 <- lapply(df.list, function(x){
  d1 <- DESeqDataSetFromMatrix(countData = x,
                               colData = conds.df,
                               design = ~condition)  
  d2 <- DESeq(d1)
  d3 <- lfcShrink(d2,
                  coef = 2,
                  type = "apeglm")
  #d3 <- results(d2)
  # return the results as a data frame with no rownames
  d4 <- as_tibble(d3, rownames = "target_id")
  return(d4)
})

deseq.df2 <- bind_rows(deseq.list2, .id ="sample") %>% 
  separate(sample, into = c("seqtype", "line"), sep = "_")
write_csv(deseq.df2, "Ara-6_deseq2_results.csv")

cat Ara*deseq2*csv | awk -F ',' '($1=="rna"&&$8<0.05&&$5>0)' | grep -v "ECB_" | cut -f 3 -d ',' | sort -u > 366_upregs

#BUILD GMAP

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6
do
cat *evolved*nongenic*gtf | grep -f 366_upregs | grep "$i" > "$i"_upreg.gtf
gffread -E -w "$i"_upreg.faa -g ../../all_genomes/"$i"*50000*fasta "$i"_upreg.gtf
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 "$i"_upreg.faa > "$i"_upreg.gff3 #map to evolved genomes
done

#Even after getting rid of duplicates and cases that ahd >10 reads in ancestor, there are 262 cases of increase.
#Of these, 249 had a log2fold change exceeding 2
#What other filter could be used?
#Replicate difference could be a huge issue
#Readthrough
#Promoter formation could be the major one
#In my original paper I only found a very few cases of promoter formation or transfer
#This could be a case where I can demonstrate this!
#That's why this is promising
#Do this with other length windows too and see what I get
