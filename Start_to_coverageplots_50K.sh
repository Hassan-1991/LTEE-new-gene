#This code describes the proto-gene pipeline until coverage plot generation
#Dataset: 50K

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
grep -v "#" ../"$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$3,"\t",$3+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_plus_down_500\";"}' > "$i"_evolved.gtf
grep -v "#" ../"$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$4-500,"\t",$4,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_minus_up_500\";"}' >> "$i"_evolved.gtf
grep -v "#" ../"$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5,"\t",$5+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_plus_down_500\";"}' > "$i"_ancestor.gtf
grep -v "#" ../"$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_minus_up_500\";"}' >> "$i"_ancestor.gtf
grep -v "#" ../"$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5+$6,"\t",$5+$6+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_plus_down_500\";"}' >> "$i"_ancestor.gtf
grep -v "#" ../"$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_minus_up_500\";"}' >> "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_ancestor.gtf -o "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_evolved.gtf -o "$i"_evolved.gtf
done

#manually check to make sure no bad records - cat *gtf | sort -nk4

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
bedtools intersect -s -a "$i"_evolved.gtf -b ../../../converted_gtf_files/"$i"_50000gen*gtf -wo | awk -F '\t' '(($12=="CDS"||$12=="gene"||$12=="repeat_region"||$12~"RNA"||$12~"fCDS")&&($19>10))' | cut -f 2 -d "\"" > "$i"_genic_windows.txt
grep -v -f "$i"_genic_windows.txt "$i"_evolved.gtf > "$i"_evolved_nongenic.gtf
grep -v -f "$i"_genic_windows.txt "$i"_ancestor.gtf > "$i"_ancestor_nongenic.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
cat ../../../insremoved_gtf/"$i"_50000gen*insremoved.gtf ../../../insremoved_gtf/REL606_insremoved.gtf | cut -f 2 -d "\"" | sort | uniq -c | grep " 2 " | rev | cut -f 1 -d ' ' | rev > "$i"_common_annotated.txt
grep -f "$i"_common_annotated.txt ../../../insremoved_gtf/REL606*gtf > "$i"_ancestor_htseqready.gtf
grep -f "$i"_common_annotated.txt ../../../insremoved_gtf/"$i"_50000gen*gtf > "$i"_evolved_temp.gtf
cat "$i"_ancestor_nongenic.gtf >> "$i"_ancestor_htseqready.gtf
cat "$i"_evolved_nongenic.gtf >> "$i"_evolved_temp.gtf #gene names need to match
paste "$i"_ancestor_htseqready.gtf "$i"_evolved_temp.gtf | awk -F '\t' '{OFS=FS}{print $10,$11,$12,$13,$14,$15,$16,$17,$9}' > "$i"_evolved_htseqready.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
ls ../../../../50k_bowtie_indices_bamfiles/*sorted.bam | egrep "REL60" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_ancestor_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' > "$i"_htseq_500.sh
ls ../../../../50k_bowtie_indices_bamfiles/*sorted.bam | grep "$i" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_evolved_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' >> "$i"_htseq_500.sh
done

#[put em in master.sh, edit as such:
#sed "s/+500/+100/g" master.sh | sed "s/-500/-100/g" | sed "s/down_500/down_100/g" | sed "s/up_500/up_100/g" | sed "s/htseq_500/htseq_100/g"]

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

variable_names <- c("Ara+1", "Ara+2", "Ara+3", "Ara+4", "Ara+5", "Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara-5", "Ara-6")

for (variable in variable_names) {
  kdf2 <- read_tsv(paste0(variable, "_500_htseq.tsv")) %>% 
    select(repl, seqtype, line, target_id, est_counts) %>% 
    unite("sample", c(seqtype, line, repl), sep = "_")
  
  line.names <- unique(str_extract(unique(kdf2$sample), "[a-z]{3,4}_[A-Za-z]+[-+][1-6]"))
  line.names <- line.names[!is.na(line.names)]
  names(line.names) <- line.names
  
  wide.kdf2 <- kdf2 %>% 
    pivot_wider(names_from = sample, values_from = est_counts)
  
  wide.kdf2 <- wide.kdf2 %>% 
    rowwise() %>% 
    mutate(row.sum = sum(c_across(starts_with("rna")))) %>% 
    filter(row.sum > 5) %>% 
    select(1:13)
  
  df.list <- lapply(line.names, function(x) {
    seqtype <- unlist(str_split(x, "_"))[1]
    anc.to.pick <- paste(seqtype, "REL", sep = "_")
    
    wide.kdf2 %>%
      select(target_id, starts_with(x), starts_with(anc.to.pick)) %>% 
      mutate_if(is.double, as.integer) %>% 
      column_to_rownames("target_id") %>% 
      as.matrix()
  })
  
  conds.df <- data.frame(condition = factor(c(rep("evo", 2), rep("anc", 4)), levels = c("anc", "evo")))
  
  deseq.list2 <- lapply(df.list, function(x) {
    d1 <- DESeqDataSetFromMatrix(countData = x,
                                 colData = conds.df,
                                 design = ~condition)  
    d2 <- DESeq(d1)
    d3 <- lfcShrink(d2,
                    coef = 2,
                    type = "apeglm")
    d4 <- as_tibble(d3, rownames = "target_id")
    return(d4)
  })
  
  deseq.df2 <- bind_rows(deseq.list2, .id ="sample") %>% 
    separate(sample, into = c("seqtype", "line"), sep = "_")
  
  write_csv(deseq.df2, paste0(variable, "_deseq2_results.csv"))
}

cat 100bp/*deseq*csv 200bp/*deseq*csv | awk -F ',' '($1=="rna"&&$8<0.05&&$5>0)' | grep -v "ECB_" | cut -f 3 -d ',' | sort -u > upregs

#BUILD GMAP

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6
do sed "s/REL606/"$i"/g" /stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/all_genomes/"$i"*50000gen*fasta >> all_genomes_modded.faa
done
cat /stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/all_genomes/REL606.fasta >> all_genomes_modded.faa
cat 100bp/*evolved*htseqread*gtf 200bp/*evolved*htseqread*gtf | grep -f upregs | gtf2bed | awk '{OFS=""}{print $4,"_",$0}' | cut -f 1,7- -d "_" | sed "s/_REL606\t/\t/g" | bedtools getfasta -fi all_genomes_modded.faa -bed - -s -name > upregs.faa

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6
do
grep --no-group-separator -A1 "$i" upregs.faa > "$i"_upregs.faa
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 "$i"_upregs.faa > "$i"_upregs.gff3 #map to evolved genomes
done

cat *upregs.gff3 | grep "mRNA.*path2" | cut -f 1 -d '(' | cut -f 2 -d "=" | sort -u > duplicates

#Manually check whether closely neighboring mutations (within <20 bp of each other) have been double-counted, and remove them:

cat 100bp/*evolved*htseqread*gtf 200bp/*evolved*htseqread*gtf | grep -f upregs | gtf2bed | awk '{OFS=""}{print $4,"_",$0}' | cut -f 1,7- -d "_" | sed "s/_REL606\t/\t/g" | sort -nk 4

Ara+3_986_DEL
Ara+3_987_INS
Ara+3_205_SNP
Ara-2_164_INS
Ara-4_770_DEL
Ara-4_771_INS

#Manually add them to the "duplicates" file.

#Finally, remove cases of overlap with repeat regions on the opposite strand:

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do sed "s/^/"$i"\t/g" /stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/converted_gtf_files/"$i"*50000*gtf | grep "repeat_region" | cut -f 1,3- >> all_repeat_regions.gtf; done
cat 100bp/*evolved*htseqread*gtf 200bp/*evolved*htseqread*gtf | grep -f upregs | gtf2bed | awk '{OFS=""}{print $4,"_",$0}' | cut -f 1,7- -d "_" | sed "s/_REL606\t/\t/g" | sort -nk 4 | bedtools intersect -wo -a - -b all_repeat_regions.gtf | awk -F '\t' '($20>90)' | cut -f 2 -d "\"" | rev | cut -f 2- -d "_" | rev | sort -u >> duplicates

grep -v -f duplicates upregs | rev | cut -f 3- -d "_" | sort -u | rev > 63_coverage_candidates.txt
grep -v -f duplicates upregs | rev | cut -f 2- -d "_" | sort -u | rev > 63_coverage_candidates)_withoutstrand.txt

#Make coverage plots from this point on.
