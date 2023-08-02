#First, make bedfiles out of the mutations for easier grepping

###########FINISH THE REST TOMORROW##########

#To make coverage plots out of paired-end data, bam files with only R1 reads need to be generated

for j in 1 2 3; do for i in $(ls ../indices/ | cut -f 1 -d '.' | sort -u | sed "s/Ara-3_//g" | grep -v "REL606"); do bowtie2 -x ../indices/Ara-3_"$i".fasta -U "$i"*BiolRep"$j"*_paired --very-sensitive-local -p 72 -S "$i"_BiolRep"$j".sam; done; done
for j in 1 2; do for i in 1 2 3; do bowtie2 -x ../indices/REL606.fasta -U REL606"$j"_BiolRep"$i"_Cit-_R1_paired --very-sensitive-local -p 72 -S REL606"$j"_BiolRep"$i".sam; done; done
for i in $(ls *sam | cut -f 1 -d '.'); do samtools view -b -S -@ 72 "$i".sam > "$i".bam; done
for i in $(ls *bam | cut -f 1 -d '.'); do samtools sort -@ 72 -o "$i"_sorted.bam "$i".bam; done
ls *gen* | awk '{OFS=""}{print "mv ",$1," Ara-3_",$1}' | bash

#First, make a repo of all Ara-3 mutations

for i in Ara-3_5000gen_ZDB409 Ara-3_10000gen_ZDB425 Ara-3_15000gen_ZDB445  Ara-3_20000gen_ZDB467  Ara-3_25000gen_ZDB483  Ara-3_25000gen_ZDB488 Ara-3_25000gen_ZDB486 Ara-3_27000gen_ZDB309 Ara-3_30000gen_ZDB17 Ara-3_31500gen_ZDB200 Ara-3_31500gen_ZDB199 Ara-3_31500gen_ZDB564 Ara-3_31500gen_ZDB25 Ara-3_33000gen_CZB154 Ara-3_50000gen_11364; do grep -v "#" "$i".gd | sed "s/^/"$i"\t/g"; done | awk -F '\t' '{OFS=""}{print $1,"\t",$2,"_",$6,"_",$7,"\t",$0}' | cut -f 1,2,4- >> Ara-3_all_mutations.gd

#Make a bed file with the surrounding 200bp of each mutation leading to increase:

grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "SNP|MOB|INS" | awk -F '\t' '{OFS=""}{print "REL606\t",$2,"\t",$2+200,"\t",$1,"_plus_down\t0\t+"}' > 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "SNP|MOB|INS" | awk -F '\t' '{OFS=""}{print "REL606\t",$2-200,"\t",$2,"\t",$1,"_plus_up\t0\t+"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "SNP|MOB|INS" | awk -F '\t' '{OFS=""}{print "REL606\t",$2-200,"\t",$2,"\t",$1,"_minus_down\t0\t-"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "SNP|MOB|INS" | awk -F '\t' '{OFS=""}{print "REL606\t",$2,"\t",$2+200,"\t",$1,"_minus_up\t0\t-"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "CON|DEL" | awk -F '\t' '{OFS=""}{print "REL606\t",$2+$3,"\t",$2+$3+200,"\t",$1,"_plus_down\t0\t+"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "CON|DEL" | awk -F '\t' '{OFS=""}{print "REL606\t",$2-200,"\t",$2,"\t",$1,"_plus_up\t0\t+"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "CON|DEL" | awk -F '\t' '{OFS=""}{print "REL606\t",$2-200,"\t",$2,"\t",$1,"_minus_down\t0\t-"}' >> 12_anc.bed
grep -f 12_total_increases.txt *gd | cut -f 2,7,8 | sort -u | egrep "CON|DEL" | awk -F '\t' '{OFS=""}{print "REL606\t",$2+$3,"\t",$2+$3+200,"\t",$1,"_minus_up\t0\t-"}' >> 12_anc.bed

#Only retain ones with correct strand, get fasta file, run gmap

grep -f 12_increases_withstrand.txt 12_anc.bed > test1 && mv test1 12_anc.bed
bedtools getfasta -fi ./REL606.fasta -s -name -bed 12_anc.bed > 12_anc.faa

###########Gmap######### (also to get rid of duplicates, manscan and gmap)

#Following gmap, the goal is to extract their coordinates in each subsequent time point

grep "something" 12_anc.bed > something.bed
grep "something" *gff3 | grep "mRNA.*mrna1" | cut -f 2,4,5,7,9 | cut -f 1 -d '(' | sed "s/ID=//g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,"0",$4}' >> something.bed

X=$(grep "REL606.*down" something.bed | cut -f 3); Y=$(grep "REL606.*down" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$4>X && $4<X+200 && $6==Y { print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $4-X}' ../coverage/timeseries/*REL60* | sed "s/_Biol/\tBiol/g" > something_down_coverage_REL606
X=$(grep "REL606.*up" something.bed | cut -f 2); Y=$(grep "REL606.*up" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$4>X && $4<X+200 && $6==Y { print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $4-X}' ../coverage/timeseries/*REL60*  | sed "s/_Biol/\tBiol/g" > something_up_coverage_REL606

for i in Ara-3_5000gen_ZDB409 Ara-3_10000gen_ZDB425 Ara-3_15000gen_ZDB445 Ara-3_20000gen_ZDB467 Ara-3_25000gen_ZDB483 Ara-3_25000gen_ZDB486 Ara-3_25000gen_ZDB488 Ara-3_27000gen_ZDB309 Ara-3_30000gen_ZDB17 Ara-3_31500gen_ZDB199 Ara-3_31500gen_ZDB200 Ara-3_31500gen_ZDB25 Ara-3_31500gen_ZDB564 Ara-3_33000gen_CZB154
do 
X=$(grep ""$i".*down" something.bed | cut -f 3); Y=$(grep ""$i".*down" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$4>X && $4<X+200 && $6==Y { print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $4-X}' ../coverage/timeseries/*"$i"* | sed "s/_Biol/\tBiol/g" > something_down_coverage_"$i"
X=$(grep ""$i".*up" something.bed | cut -f 2); Y=$(grep ""$i".*up" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$4>X && $4<X+200 && $6==Y { print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $4-X}' ../coverage/timeseries/*"$i"*  | sed "s/_Biol/\tBiol/g" > something_minus_up_coverage_"$i"
done

X=$(grep "Ara-3_50000gen_11364.*down" something.bed | cut -f 3); Y=$(grep "Ara-3_50000gen_11364.*down" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' ../coverage/*rna*Ara-3* | sed "s/-rna-/\t/g" |  awk -F '\t' '{OFS=FS}{print $2,$1,$3,$4,$5}' > something_minus_down_coverage_Ara-3_50000gen_11364
X=$(grep "Ara-3_50000gen_11364.*up" something.bed | cut -f 2); Y=$(grep "Ara-3_50000gen_11364.*up" something.bed | cut -f 6); awk -F '\t' -v X="$X" -v Y="$Y" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' ../coverage/*rna*Ara-3* | sed "s/-rna-/\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1,$3,$4,$5}' > something_minus_up_coverage_Ara-3_50000gen_11364

