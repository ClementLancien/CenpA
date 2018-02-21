#!/bin/bash
 
sample=$1
rawdir=$2
outdir=$3

R1="SRR633614.R1.fastq.gz"
R2="SRR633614.R1.fastq.gz"

tmp="SRR633614.tmp.fastq"

if [ ! -f $tmp ]; then

    paste <(gzcat $R1) <(gzcat $R2) | awk '{printf("%s", $0); n++; if (n%4 == 0) printf("\n"); else printf("\t\t");}' > $tmp

fi

exit
nreads=10000000

for i in `seq 1 5`
    
    do
    
        out1="$outdir/$sample.$i.R1.fastq"
        out2="$outdir/$sample.$i.R2.fastq"
        
        shuf $tmp | head -n $nreads | sed 's/\t\t/\n/g' | awk -v out1="$out1" -v out2="$out2" -F '\t' '{print $1 > out1; print $2 > out2}'
                
        gzip $out1
        gzip $out2
        
    done

if [ $(ls $outdir/$sample.*.gz | wc -l) -eq $(( 5 * 2 )) ] && [ -f $tmp ]; then
       
    rm $tmp

fi

exit