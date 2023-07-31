#50K dataset

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/^/"$i"\t/g" | awk -F '\t' '{OFS=""}{print $1,"_",$3,"_",$2,"\t",$0}'; done | cut -f 1,6- > 50k_ancestorcoords.gd
for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do grep -v "#" "$i"_50000gen*applied.gd | sed "s/^/"$i"\t/g" | awk -F '\t' '{OFS=""}{print $1,"_",$3,"_",$2,"\t",$0}'; done > 50k_applied.gd

#Making the coordinates easier to extract by putting them in bed files
#The easiest way to do this is specify mutation start and end coordinates on both strands, and then just retain the proto-gene-relevant strand

grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4,"\t",$2,"_plus\t",0,"\t+"}' | grep -v "DEL" > ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4,"\t",$2,"_minus\t",0,"\t-"}' | grep -v "DEL" >> ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | grep "DEL" | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4+$5,"\t",$2,"_plus\t",0,"\t+"}' >> ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | grep "DEL" | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4+$5,"\t",$2,"_minus\t",0,"\t-"}' >> ancestor_upregs
grep -f 63_coverage_candidates.txt ancestor_upregs > ancestor_upregs.bed

grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep "SNP|DEL" | awk -F '\t' '{OFS=""}{print $7,"\t",$8,"\t",$8,"\t",$2,"_plus\t0\t+"}' > evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep "SNP|DEL" | awk -F '\t' '{OFS=""}{print $7,"\t",$8,"\t",$8,"\t",$2,"_minus\t0\t-"}' >> evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep -v "SNP|DEL" | awk '{OFS="\t"}{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk '{OFS=""}{print $1,"\t",$4,"\t",$3,"\t",$2,"_plus\t0\t+"}' >> evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep -v "SNP|DEL" | awk '{OFS="\t"}{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk '{OFS=""}{print $1,"\t",$4,"\t",$3,"\t",$2,"_minus\t0\t-"}' >> evolved_upregs
grep -f 63_coverage_candidates.txt evolved_upregs > evolved_upregs.bed

#Make genome coverage files for each RNA and Riboseq run

for i in $(ls /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/*sorted.bam | rev | cut -f 1 -d '/' | cut -f 2- -d "_" | rev | sort -u);
do
bedtools genomecov -d -strand + -ibam /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/"$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t+/g" | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" > "$i"_plus
bedtools genomecov -d -strand - -ibam /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/"$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t-/g" | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" > "$i"_minus
cat "$i"_plus "$i"_minus > "$i"_coverage.tsv
done

#Make coverage plots for each mutation

mapfile -t values < 63_coverage_candidates.txt

#ancestor_right_200

for value in "${values[@]}"
do
    X_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 3)
    Y_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_ancestor" -v Y="$Y_ancestor" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' coverage/*REL60* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_ancestor"_right
done

#ancestor_left_200

for value in "${values[@]}"
do
    X_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 2)
    Y_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_ancestor" -v Y="$Y_ancestor" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*REL60* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_ancestor_left"
done

#evolved_right_200

for value in "${values[@]}"
do
    X_evolved=$(grep "$value" evolved_upregs.bed | cut -f 3)
    Y_evolved=$(grep "$value" evolved_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_evolved" -v Y="$Y_evolved" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' coverage/*"$(cut -f 1 -d "_" <<< "$value")"* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_right"
done

#evolved_Left_200

for value in "${values[@]}"
do
    X_evolved=$(grep "$value" evolved_upregs.bed | cut -f 2)
    Y_evolved=$(grep "$value" evolved_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_evolved" -v Y="$Y_evolved" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*"$(cut -f 1 -d "_" <<< "$value")"* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_left"
done

#For coverage plots in the minus strand, one more fix is required

for i in $(grep "minus" ../63_coverage_candidates.txt | cut -f 1-3 -d "_"); do awk -i inplace '{OFS='\t'}{$6 = 200 - $6}1' "$i"*coverage*; done
sed -i "1s/^/replicate\tseqtype\tline\tcount\tposition_abs\tposition\n/g" *coverage

#Generating plots in R

