#This code describes how transcriptomes were constructed containing mutation-adjacent regions for each evolved genome.
#Mutation information for each population are available in the form of .gd files from here: https://zenodo.org/record/7447457

#Convert the mutation information into gtf files, which picks out adjacent regions of the mutations in the ancestral genome sequence:

for i in Ara+1_50000gen_11392 Ara+2_50000gen_11342 Ara+3_50000gen_11345 Ara+4_50000gen_11348 Ara+5_50000gen_11367 Ara-1_50000gen_11330 Ara-2_50000gen_11333 Ara-3_50000gen_11364 Ara-4_50000gen_11336 Ara-5_50000gen_11339 Ara-6_50000gen_11389;
do
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "MOB" | awk -F ',' '($4=="1")' | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2,"\t",$2+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$0,"\";gene_id \"",var,"_",$0,"\";"}' | sed "s/,/_/g" > ancestor_"$i"_mutation_adj_regions.gtf
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "MOB" | awk -F ',' '($4=="-1")' | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2-500,"\t",$2,"\t.\t-\t0\ttranscript_id \"",var,"_",$0,"\";gene_id \"",var,"_",$0,"\";"}' | sed "s/,/_/g" >> ancestor_"$i"_mutation_adj_regions.gtf
for j in SNP INS
do
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "$j" | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2,"\t",$2+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$0,"\_plus_down\";gene_id \"",var,"_",$0,"\_plus_down\";"}' | sed "s/,/_/g" >> ancestor_"$i"_mutation_adj_regions.gtf
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "$j" | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2-500,"\t",$2,"\t.\t-\t0\ttranscript_id \"",var,"_",$0,"\_minus_up\";gene_id \"",var,"_",$0,"\_minus_up\";"}' | sed "s/,/_/g" >> ancestor_"$i"_mutation_adj_regions.gtf
for k in DEL CON SUB AMP
do
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "$k" | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2+$3,"\t",$2+$3+500,"\t.\t+\t0\ttranscript_id \"",var,"_",$0,"_plus_down\";gene_id \"",var,"_",$0,"_plus_down\";"}' | sed "s/,/_/g" >> ancestor_"$i"_mutation_adj_regions.gtf
grep -v "^#" "$i".gd | cut -f 1,5-7 | sed "s/\t/,/g" | grep "$k" | awk -v var="$i" -F ',' '{OFS=""}{print "REL606\t.\tCDS\t",$2-500,"\t",$2,"\t.\t-\t0\ttranscript_id \"",var,"_",$0,"_minus_up\";gene_id \"",var,"_",$0,"_minus_up\";"}' | sed "s/,/_/g" >> ancestor_"$i"_mutation_adj_regions.gtf
done
done
done

for i in *gtf; do sort -u -o $i $i; done
for i in *gtf; do sort -nk 4 -o $i $i; done

#manually remove records that start with values <0, or replace them with "1". Check to see which files contain such records:

cat *.gtf | sort -k4,4n | less

#Convert the .gtf files to .faa files using gffread

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do gmap_build -D . -d "$i" all_genomes/"$i".fasta; done

for i in Ara+1_50000gen_11392 Ara+2_50000gen_11342 Ara+3_50000gen_11345 Ara+4_50000gen_11348 Ara+5_50000gen_11367 Ara-1_50000gen_11330 Ara-2_50000gen_11333 Ara-3_50000gen_11364 Ara-4_50000gen_11336 Ara-5_50000gen_11339 Ara-6_50000gen_11389;
do
gffread -E -w ancestor_"$i"_mutation_adj_regions.faa -g ../../../../all_genomes/"$i".fasta ancestor_"$i"_mutation_adj_regions.gtf
gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 ancestor_"$i"_mutation_adj_regions.faa > "$i"_mutation_adj_regions.gff3
grep "mRNA" "$i"_mutation_adj_regions.gff3 | grep "path2" | cut -f 1 -d ';' | cut -f 9 | sed "s/ID=//g" | sed "s/.mrna2//g" > "$i"_duplicates.txt
grep -v -f "$i"_duplicates.txt "$i"_mutation_adj_regions.gff3 | grep "path1" | grep "mRNA" | sed "s/mRNA/CDS/g" | sed "s/ID=/transcript_id \"/g" | sed "s/.mrna1;Name=/\";gene_id \"/g" | sed "s/;Parent/\";Parent/g"| cut -f 1,2 -d ";" | sed "s/$/;/g" > interim/"$i"_mutation_adj_regions.gtf
cp interim/"$i"_mutation_adj_regions.gtf .
sed "s/+/%/g" interim/"$i"_mutation_adj_regions.gtf | sed "s/\t-\t/\t+\t/g" | sed "s/%/-/g" | sed "s/id \"/id \"control_/g" >> "$i"_mutation_adj_regions.gtf
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10=$5-$4}' "$i"_mutation_adj_regions.gtf | awk -F '\t' '($10<600)' > "$i"_mutation_adj_regions_final.gtf
grep -v "control" "$i"_mutation_adj_regions_final.gtf | cut -f 2 -d "\"" > "$i"_IDs.txt
grep -f "$i"_IDs.txt ancestor_"$i"_mutation_adj_regions_final.gtf | sort -u > interim/ancestor_"$i"_mutation_adj_regions.gtf
cp interim/ancestor_"$i"_mutation_adj_regions.gtf .
sed "s/\t+\t/\t%\t/g" interim/ancestor_"$i"_mutation_adj_regions.gtf |  sed "s/\t-\t/\t+\t/g" | sed "s/%/-/g" | sed "s/id \"/id \"control_/g" | sort -u > ancestor_"$i"_mutation_adj_regions_final.gtf
done

for file in *.gtf; do sort -u -o $file $file; done

for i in Ara+1_50000gen_11392 Ara+2_50000gen_11342 Ara+3_50000gen_11345 Ara+4_50000gen_11348 Ara+5_50000gen_11367 Ara-1_50000gen_11330 Ara-2_50000gen_11333 Ara-3_50000gen_11364 Ara-4_50000gen_11336 Ara-5_50000gen_11339 Ara-6_50000gen_11389;
do
bedtools intersect -a ../../../insremoved_gtf/"$i"*gtf -b "$i"*gtf -v -s | cut -f 2 -d "\"" | sort -u > "$i"_coding_ids.txt
bedtools intersect -a ../../../insremoved_gtf/REL606_insremoved.gtf -b ancestor_"$i"*gtf -v -s | cut -f 2 -d "\"" | sort -u  > ancestor_"$i"_coding_ids.txt
cat ancestor_"$i"_coding_ids.txt "$i"_coding_ids.txt | sort | uniq -c | grep " 2 " | rev | cut -f 1 -d ' ' | rev | sed "s/^/\"/g" | sed "s/$/\"/g" > "$i"_common_IDs.txt
grep -f "$i"_common_IDs.txt ../../../insremoved_gtf/"$i"*gtf | sort -u >> "$i"_mutation_adj_regions_4b.gtf
grep -f "$i"_common_IDs.txt ../../../insremoved_gtf/REL606_insremoved.gtf | sort -u >> ancestor_"$i"*gtf
done
