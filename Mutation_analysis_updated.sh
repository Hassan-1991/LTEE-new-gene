#For 500:

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
grep -v "#" "$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || i == 5 || i == 6 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$6,"\t",$6+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_plus_down_500\";"}' > "$i"_evolved.gtf
grep -v "#" "$i"_50000gen-applied.gd | grep -v "INV" | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | awk '{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || i == 5 || i == 6 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk -v var="$i" -F ' ' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$3,"_",$4,"_minus_up_500\";"}' >> "$i"_evolved.gtf
grep -v "#" "$i"_50000gen-applied.gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5,"\t",$5+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";"}' > "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen-applied.gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "SNP|INS|MOB" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";"}' >> "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen-applied.gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5+$6,"\t",$5+$6+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_plus_down_500\";"}' >> "$i"_ancestor.gtf
grep -v "#" "$i"_50000gen-applied.gd | sed "s/IS1\t-1/IS1_-1/g" | sed "s/IS1\t1/IS1_1/g" | sed "s/IS3\t-1/IS3_-1/g" | sed "s/IS3\t1/IS3_1/g" | sed "s/IS4\t-1/IS4_-1/g" | sed "s/IS4\t1/IS4_1/g" | sed "s/IS150\t-1/IS150_-1/g" | sed "s/IS150\t1/IS150_1/g" | sed "s/IS186\t-1/IS186_-1/g" | sed "s/IS186\t1/IS186_1/g" | egrep "DEL|CON|SUB|AMP" | awk -v var="$i" -F '\t' '{OFS=""}{print "REL606\t.\tCDS\t",$5-500,"\t",$5,"\t.\t-\t0\ttranscript_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";gene_id \"",var,"_",$2,"_",$1,"_",$5,"_",$6,"_minus_up_500\";"}' >> "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_ancestor.gtf -o "$i"_ancestor.gtf
sort -t "_" -nk 3 "$i"_evolved.gtf -o "$i"_evolved.gtf
done

#manually check to make sure no bad records

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
bedtools intersect -s -v -a "$i"_evolved.gtf -b ../../converted_gtf_files/"$i"_50000gen*gtf -wo | cut -f 2 -d "\"" > "$i"_nongenic_windows.txt
grep -f "$i"_nongenic_windows.txt "$i"_evolved.gtf > "$i"_evolved_nongenic.gtf
grep -f "$i"_nongenic_windows.txt "$i"_ancestor.gtf > "$i"_ancestor_nongenic.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
cat ../../insremoved_gtf/"$i"_50000gen*insremoved.gtf ../../insremoved_gtf/REL606_insremoved.gtf | cut -f 2 -d "\"" | sort | uniq -c | grep " 2 " | rev | cut -f 1 -d ' ' | rev > "$i"_common_annotated.txt
grep -f "$i"_common_annotated.txt ../../insremoved_gtf/REL606*gtf > "$i"_ancestor_htseqready.gtf
grep -f "$i"_common_annotated.txt ../../insremoved_gtf/"$i"_50000gen*gtf > "$i"_evolved_htseqready.gtf
cat "$i"_ancestor_nongenic.gtf >> "$i"_ancestor_htseqready.gtf
cat "$i"_evolved_nongenic.gtf >> "$i"_evolved_htseqready.gtf
done

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
ls ../../../50k_bowtie_indices_bamfiles/*sorted.bam | egrep "REL60" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_ancestor_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' > "$i"_htseq_500.sh
ls ../../../50k_bowtie_indices_bamfiles/*sorted.bam | grep "$i" | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | awk -v var="$i" '{OFS=""}{print $0," ",var,"_evolved_htseqready.gtf \| head -n -5 \> T",$11}' | sed "s/T..\/..\/..\/50k_bowtie_indices_bamfiles\///g" | sed "s/$/t/g" | sed "s/_sorted.bamt//g" | awk -v var=$i '{OFS=""}{print $0,"_",var}' >> "$i"_htseq_500.sh
done
