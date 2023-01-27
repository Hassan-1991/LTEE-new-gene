awk -F '\t' '($3=="mRNA")' Ara+1_50000gen_11392_transterm_motifs_100bpoffset.gff3 | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \""$9,"\";transcript_id \""$9,"\";"}' | sed "s/ID=//g" | sed "s/\t\";/\";/g" > Ara+1_50000gen_11392_transterm_motifs_100bpoffset.gtf

#For all:

#put this in basic_code.sh using nano:

awk -F '\t' '($3=="mRNA")' something.gff3 | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \""$9,"\";transcript_id \""$9,"\";"}' | sed "s/ID=//g" | sed "s/\t\";/\";/g" > something.gtf
ls Ara*gff3 | rev | cut -f 2- -d '.' | rev > somethings
sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh
bash interim_code.sh
bash final_code.sh
