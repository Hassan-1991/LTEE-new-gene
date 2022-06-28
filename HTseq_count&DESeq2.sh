#Write a script to run HTseq on all lines/strains
#We're making two files, one for the time series dataset, and one for the 50000gen one

ls *commonwindows.gtf >> test1
ls *commonwindows.gtf >> test1
ls *commonwindows.gtf >> test1
grep -v "50000" test1 > time_series_test
sort -t '_' -nk 2 -o time_series_test time_series_test
ls *commonwindows.gtf >> test1
egrep "50000|REL60" test1 > test2
sort -u test2 > test3
cat test2 test3 > 50000gen_test
sort -t '_' -nk 2 -o 50000gen_test 50000gen_test
rm test*

#For time series
ls ../../../../trimmomatic_29bp_trimmed/bamfiles/*sorted.bam | sed "s/gen/gen_/g" | sort -t '_' -nk 4 | sed "s/gen_/gen/g" > time_series_bams
paste time_series_bams time_series_test | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | sed "s/\t/ /g" | sed "s/$/ \| head -n -5 /g" > time_series_htseqstep3
rev time_series_bams | cut -f 1 -d '/' | rev | cut -f 1-3 -d '_' | sed "s/_/,/g" > time_series_htseqstep4
paste time_series_htseqstep3 time_series_htseqstep4 | sed "s/\t/ \| sed \"s\/\^\//g" | sed "s/$/,\/g\"/g" | sed "s/$/ \| sed \"s\/TAB\/\,\/g\" >> htseq_time_series_400bp_count.csv/g" | sed "s/TAB/\t/g" > htseq_time_series_400bp_code.sh
bash htseq_time_series_400bp_code.sh

#For 50000gen
ls ../../../../50k_bowtie_indices_bamfiles/*_sorted.bam | sort -t '_' -nk 4 | sed "s/R/1_R/g" | sed "s/P/2_P/g" | sed "s/M/3_M/g" | sed "s/Ara/Ara_/g" | rev | sort -t "_" -nk 3 | rev | sed "s/_1_//g" | sed "s/_2_//g" | sed "s/_3_//g" > 50000gen_bams
paste 50000gen_bams 50000gen_test | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | sed "s/\t/ /g" | sed "s/$/ \| head -n -5 /g" > 50000gen_htseqstep3
rev 50000gen_bams | cut -f 1 -d '/' | rev | cut -f 1 -d '_' | sed "s/-/,/g" | sed "s/am//g" | sed "s/ap//g" | sed "s/AraR6/REL606/g" | sed "s/AraR7/REL607/g" > 50000gen_htseqstep4
paste 50000gen_htseqstep3 50000gen_htseqstep4 | sed "s/\t/ \| sed \"s\/\^\//g" | sed "s/$/,\/g\"/g" | sed "s/$/ \| sed \"s\/TAB\/\,\/g\" >> htseq_50000gen_400bp_count.csv/g" | sed "s/TAB/\t/g" > htseq_50000gen_400bp_code.sh
bash htseq_50000gen_400bp_code.sh
