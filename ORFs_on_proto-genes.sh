#Goal: Extract all open reading frames encoded in the -200bp to +500bp region surrounding proto-gene-causing mutations
#Start codons: ATG, GTG, TTG
#Minimum length: 30bp
#Same strand as the direction of novel transcription
#Exclude ORFs contained entirely within the upstream region of the mutation

#Copy-paste the names of 9 proto-genes to the file 9_proto_genes.txt
#63_coverage_candidates.txt is a list of proto-gene candidates extracted from the DeSEQ2 output
#evolved_seqs.gtf contains the +200bp regions of all mutations in evolved genomes

#Get the -200, +500 regions from this gtf
#Done separately for plus- and minus-stranded proto-genes
grep -f 9_proto_genes.txt 63_coverage_candidates.txt |
grep -f - evolved_seqs.gtf |
awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-200,$5+300,$6,$7,$8,$9}' > 9_proto_genes_plus500bp_minus200bp.gtf

grep -f 9_proto_genes.txt 63_coverage_candidates.txt |
grep -f - evolved_seqs.gtf |
awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-300,$5+200,$6,$7,$8,$9}' >> 9_proto_genes_plus500bp_minus200bp.gtf

#Extract their DNA sequences from the genome file
#all_genomes.fasta contains all LTEE genomes under consideration
#Modify the first column of the gtf file so it contains the correct genome name

cat 9_proto_genes_plus500bp_minus200bp.gtf | gtf2bed |
bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/all_genomes/all_genomes.fasta -bed - > 9_proto_genes_plus500bp_minus200bp.faa

#Extract all open reading frames from these regions with getorf
getorf -sequence 9_proto_genes_plus500bp_minus200bp.faa -outseq 9_proto_genes_plus500bp_minus200bp.getorf -table 1 -minsize 30 -find 3

#Filter by criteria mentined above
cat 9_proto_genes_plus500bp_minus200bp.getorf | paste - - |
grep -v "REVERSE" | #Get rid of opposite strand ORFs
grep -P -v "\tCTG" | #Get rid of ORFs with CTG start codon
sed "s/] / ] /g" | awk '($4>199)'| #Get rid of ORFs entirely contained within the upstream 200bp
cut -f 1 -d ']' | tr -d ">" | sed "s/ - / /g" | tr -d "[" | sed "s/(+)_/ /g"  | sed "s/(-)_/ /g" | cut -f 1,3,4 -d ' ' |
sed "s/ /\t/g" | sort -k1 > 9_proto_genes_plus500bp_minus200bp_ORF_lengths.tsv

#Modify the gtf file to join them together and then pull genome coordinates
#For the gtf file, we add another column that contains the proto-gene name
awk -F "\"" '{OFS=FS}{print $1,$2,$3,$4,$5,"@",$2}' 9_proto_genes_plus500bp_minus200bp.gtf |
sed "s/\"@\"/\t/g" | sort -k10 > 9_proto_genes_plus500bp_minus200bp_modified.gtf

#For plus-stranded proto-genes:
join -t "	" -1 1 -2 10 9_proto_genes_plus500bp_minus200bp_ORF_lengths.tsv 9_proto_genes_plus500bp_minus200bp_modified.gtf | #the join delimiter is a tab, not space
cut -f 1,2,3,7,8,10 | grep "plus" |
awk  -F '\t' '{OFS=FS}{print $0, $4+$2,$4+$3, $3-$2+1}' #printed on display

#For minus-stranded proto-genes:
join -t "	" -1 1 -2 10 9_proto_genes_plus500bp_minus200bp_ORF_lengths.tsv 9_proto_genes_plus500bp_minus200bp_modified.gtf | #same
cut -f 1,2,3,7,8,10 | grep "minus" |
awk -F '\t' '{OFS=FS}{print $0, $5-$2,$5-$3,$3-$2+1}' #same
