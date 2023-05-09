#Extract evidence for lack of upstream promoters

awk -F '\t' '($7=="+")' linegraph_29_IDs_REL606.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-300,$4,$6,$7,$8,$9}' > linegraph_29_IDs_REL606_upstream300.gtf
awk -F '\t' '($7=="-")' linegraph_29_IDs_REL606.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+300,$6,$7,$8,$9}' >> linegraph_29_IDs_REL606_upstream300.gtf
gffread -E -w linegraph_29_IDs_REL606_upstream300.faa -g ../../../all_genomes/REL606.fasta linegraph_29_IDs_REL606_upstream300.gtf
