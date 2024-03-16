#! /bin/bash

inDir=/scratch/melanie/Benchmarking/Quantif

#for file in $inDir/*.gtf ; do
#	/home/apps/gffread/gffread -w ${file%*.gtf}.fa -g /scratch/melanie/chrIS.fa $file
#done

#for file in $inDir/*.fa ; do
#	if [[ $file =~ LSK114 ]] ; then
#		inFile=/scratch/melanie/inputs/LSK114_chrIS_mixA_cDNA_sub150k.fastq
#	else
#		inFile=/scratch/melanie/inputs/RNA004_chrIS_mixA_dRNA_sub150k.fastq
#	fi
#	/home/apps/minimap2-2.24/minimap2 -a -t 20 -N 100 $file $inFile | /home/apps/samtools-1.17/samtools view -@ 6 -m 16G -b - | /home/apps/samtools-1.17/samtools sort -@ 6 -m 64G - > ${file%*.fa}.bam
#	/home/apps/salmon/salmon-1.10.0/bin/salmon quant -l SF -t $file -p 20 -a ${file%*.fa}.bam -o ${file%*.bam}_salmon
#done

/home/apps/salmon/salmon-1.10.0/bin/salmon quantmerge --quants *_salmon --column numreads -o Transcript_count_matrix.tsv 
