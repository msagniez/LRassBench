#!/bin/bash
#Script to run all cDNA and dRNA assembly pipelines | Guided ; De novo ; Ab initio | for one single dataset
#Fill in the 4 lines of information in the main
#Comment out sections you don't want to run

preprocess(){
#Minimap2 - Li H. 2021. New strategies to improve minimap2 alignment accuracy. Bioinformatics 37: 4572–4574
#https://github.com/lh3/minimap2 - v2.24
#Samtools - https://www.htslib.org/
#https://github.com/samtools/samtools - v1.19.2
	~/apps/minimap2/minimap2-2.24/minimap2 -ax splice -t $threads --secondary=no $reference_fa $InDir/${Sample}.fastq > $InDir/${Sample}.sam
	samtools view -b $InDir/${Sample}.sam | samtools sort - > $InDir/${Sample}.bam
	~/apps/minimap2/minimap2-2.24/minimap2 -ax splice -t $threads --secondary=no --MD $reference_fa $InDir/${Sample}.fastq > $InDir/${Sample}_MD.sam
}

baseline(){
#spliced_bam2gff - (c) 2020 Oxford Nanopore Technologies Ltd.
#https://github.com/nanoporetech/spliced_bam2gff - vJan. 2021
	mkdir -p $InDir/$Sample/control/
        cd $InDir/$Sample/control/
        echo "Starting Control"
	spliced_bam2gff -M $InDir/${Sample}.bam > ${Sample}_converted.gff
	sed 's/;0/-0/g' ${Sample}_converted.gff | sed 's/;16/-16/g' > ${Sample}_converted-corrected.gff
	rm ${Sample}.bam
	rm ${Sample}_converted.gff
}

stringtie2(){
#Stringtie2 - Pertea M, Pertea GM, Antonescu CM, Chang T-C, Mendell JT, Salzberg SL. 2015. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33: 290–295.
#https://github.com/gpertea/stringtie - v2.2.1-2.2.2
        mkdir -p $InDir/$Sample/stringtie2/
        cd $InDir/$Sample/stringtie2/
        echo "Starting Stringtie2"
	/usr/bin/time -v stringtie -p $threads -L -G $reference_gtf -o ./${Sample}_strg2def.gtf $InDir/${Sample}.bam
	/usr/bin/time -v stringtie -p $threads -L -o ./${Sample}_strg2def-noRef.gtf $InDir/${Sample}.bam
	echo "Stringtie2 terminates"
}

Flair(){
#FLAIR - Tang AD, Soulette CM, van Baren MJ, Hart K, Hrabeta-Robinson E, Wu CJ, Brooks AN. 2020. Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nat Commun 11: 1438.
#https://github.com/BrooksLabUCSC/flair - v2.0
	mkdir -p $InDir/$Sample/flair-noRef/
	cd $InDir/$Sample/flair-noRef/
	echo "Starting FLAIR"
	samtools index $InDir/${Sample}.bam $InDir/${Sample}.bam.bai
	flair 13 -r $InDir/${Sample}.fastq -g $reference_fa -q ./${Sample}_flair-noRef.bed -o ${Sample}_flair-noRef -t $threads

	mkdir -p $InDir/$Sample/flair/
	cd $InDir/$Sample/flair/

	flair 123 -r $InDir/${Sample}.fastq -g $reference_fa -f $reference_gtf -o ${Sample}_flair -t $threads
	echo "FLAIR terminates"
}

isoQuant(){
#IsoQuant - Prjibelski AD, Mikheenko A, Joglekar A, Smetanin A, Jarroux J, Lapidus AL, Tilgner HU. 2023. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol 41: 915–918.
#https://github.com/ablab/IsoQuant - v3.3.1
	mkdir -p $InDir/$Sample/isoQuant/
	cd $InDir/$Sample/isoQuant/
	echo "Starting isoQuant"

	isoquant.py -t $threads --clean_start --reference $reference_fa --genedb $reference_gtf --complete_genedb --fastq $InDir/${Sample}.fastq --data_type nanopore -o ./

	echo "isoQuant terminates"

}

isoQuant-noRef(){
#IsoQuant - Prjibelski AD, Mikheenko A, Joglekar A, Smetanin A, Jarroux J, Lapidus AL, Tilgner HU. 2023. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol 41: 915–918.
#https://github.com/ablab/IsoQuant - v3.3.1
	mkdir -p $InDir/$Sample/isoQuant-noRef/
	cd $InDir/$Sample/isoQuant-noRef/
	echo "Starting isoQuant"

	isoquant.py -t $threads --clean_start --reference $reference_fa --fastq $InDir/${Sample}.fastq --data_type nanopore -o ./

	echo "isoQuant-noRef terminates"

}

TALON(){
#TALON - Wyman D, Balderrama-Gutierrez G, Reese F, Jiang S, Rahmanian S, Forner S, Matheos D, Zeng W, Williams B, Trout D, et al. 2020. A technology-agnostic long-read analysis pipeline for transcriptome discovery and quantification. bioRxiv 672931.
#https://github.com/mortazavilab/TALON - v5.0
	mkdir -p $InDir/$Sample/talon/
	cd $InDir/$Sample/talon/
	echo "Starting TALON"

	#add a SAM label
	talon_label_reads --f=$InDir/${Sample}_MD.sam --g=$reference_fa --t=$threads --o=./${Sample}.labeled.sam
	#initialize transcriptome database
	talon_initialize_database --f $reference_gtf_talon --g $chr --a ${chr}_annot --o ./${Sample}_talon
	#Configure samples
	echo "$chr,talon_def,nanopore,$InDir/$Sample/talon/${Sample}.labeled.sam_labeled.sam" > $InDir/$Sample/talon/dataset.config
	#Run TALON
	talon --f dataset.config --db ${Sample}_talon.db --build $chr --threads $threads --o ${Sample}_TALON
	#filter
	talon_filter_transcripts --db ${Sample}_talon.db --datasets $chr -a ${chr}_annot --maxFracA 0.5 --minCount 5 --minDatasets 1 --o filtered_transcripts.csv
	#Extract gtf
	talon_create_GTF --db ${Sample}_talon.db --whitelist filtered_transcripts.csv -a ${chr}_annot --build $chr --o ${Sample}_gtf
	echo "TALON terminates"
}

TALON_reco(){
#TALON - Wyman D, Balderrama-Gutierrez G, Reese F, Jiang S, Rahmanian S, Forner S, Matheos D, Zeng W, Williams B, Trout D, et al. 2020. A technology-agnostic long-read analysis pipeline for transcriptome discovery and quantification. bioRxiv 672931.
#https://github.com/mortazavilab/TALON - v5.0
	mkdir -p $InDir/$Sample/talon_reco/
	cd $InDir/$Sample/talon_reco/
	echo "Starting TALON"

	#SAM file generation
	~/apps/minimap2/minimap2-2.24/minimap2 -t $threads -ax splice -uf -k14 --MD $reference_fa $InDir/${Sample}.fastq > $InDir/${Sample}_mm2talonreco.sam
	#TranscriptClean
	/usr/bin/time -v  transcriptclean -t $threads --sam $InDir/${Sample}_mm2talonreco.sam --genome $reference_fa --outprefix ./
	#Subset sam to simulate 3 replicates
	~/apps/samtools-1.19.2/samtools view -s 0.33 ./TC_clean.sam > ./Preprep1.sam
	~/apps/samtools-1.19.2/samtools view -H ./TC_clean.sam | cat - ./Preprep1.sam > rep1.sam
	join -v1 <( sort ./TC_clean.sam ) <( sort ./rep1.sam ) | sed 's/ /\t/g' > Prepreps.sam
	~/apps/samtools-1.19.2/samtools view -H ./TC_clean.sam | cat - ./Prepreps.sam > reps.sam
	~/apps/samtools-1.19.2/samtools view -s 0.67 ./reps.sam > ./Preprep2.sam
	~/apps/samtools-1.19.2/samtools view -H ./TC_clean.sam | cat - ./Preprep2.sam > rep2.sam
	join -v1 <( sort ./reps.sam ) <( sort ./rep2.sam ) | sed 's/ /\t/g' > Preprep3.sam
	~/apps/samtools-1.19.2/samtools view -H ./TC_clean.sam | cat - ./Preprep3.sam > rep3.sam
	rm ./Prep*.sam
	#initialize transcriptome database
	/usr/bin/time -v talon_initialize_database --f $reference_gtf_talon --g $chr --a ${chr}_annot --o ./${Sample}_talon_reco
	#add a SAM labels
	/usr/bin/time -v talon_label_reads --f=./rep1.sam --g=$reference_fa --t=$threads --o=$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep1.labeled.sam
	/usr/bin/time -v talon_label_reads --f=./rep2.sam --g=$reference_fa --t=$threads --o=$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep2.labeled.sam
	/usr/bin/time -v talon_label_reads --f=./rep3.sam --g=$reference_fa --t=$threads --o=$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep3.labeled.sam
	#Configure samples
	echo "Rep1,$chr,nanopore,$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep1.labeled.sam_labeled.sam" > ./dataset.config
	echo "Rep2,$chr,nanopore,$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep2.labeled.sam_labeled.sam" >> ./dataset.config
	echo "Rep3,$chr,nanopore,$InDir/$Sample/talon_reco/${Sample}_mm2talonreco.rep3.labeled.sam_labeled.sam" >> ./dataset.config
	#Run TALON
	/usr/bin/time -v talon --f dataset.config --db ${Sample}_talon_reco.db --build $chr --threads $threads --o ${Sample}_TALONpip
	#Filter
	/usr/bin/time -v talon_filter_transcripts --db ${Sample}_talon_reco.db --datasets Rep1,Rep2,Rep3 -a ${chr}_annot --maxFracA 0.5 --minCount 5 --minDatasets 2 --o filtered_transcripts-A05C5.csv
	#Extract gtf
	/usr/bin/time -v talon_create_GTF --db ${Sample}_talon_reco.db --whitelist filtered_transcripts-A05C5.csv -a ${chr}_annot --build $chr --o ${Sample}_gtf
	echo "TALON terminates"

}

FLAMES(){
#FLAMES - Tian L, Jabbari JS, Thijssen R, Gouil Q, Amarasinghe SL, Voogd O, Kariyawasam H, Du MRM, Schuster J, Wang C, et al. 2021. Comprehensive characterization of single-cell full-length isoforms in human and mouse with long-read sequencing. Genome Biol 22: 310.
#https://github.com/LuyiTian/FLAMES - v0.1
#GFFread - Pertea G, Pertea M. 2020. GFF Utilities: GffRead and GffCompare. F1000Res 9. http://dx.doi.org/10.12688/f1000research.23297.2.
#https://github.com/gpertea/gffread - v0.12.7
	mkdir -p $InDir/$Sample/FLAMES/
	cd $InDir/$Sample/FLAMES/
	echo "Starting FLAMES"

	mkdir -p $InDir/${Sample}_fastq
	cp $InDir/${Sample}.fastq $InDir/${Sample}_fastq/${Sample}.fastq

	/usr/bin/time -v ~/apps/FLAMES/python/bulk_long_pipeline.py --gff3 $reference_gtf --genomefa $reference_fa --outdir ./ --config_file ~/apps/FLAMES/examples/SIRV/data/SIRV_config.json --fq_dir $InDir/${Sample}_fastq

	#filter filtered.gtf on transcripts assessed as "True" in isoform_FSM_annotation.csv
	awk -F"," '$4=="True" {print "transcript:"$1}' isoform_FSM_annotation.csv > True.lst
	gffread --ids True.lst isoform_annotated.filtered.gff3 > isoform_annotated.filtered_True.gff3

	echo "FLAMES terminates"
}

Mandalorion(){
	mkdir -p $InDir/$Sample/Mandalorion/
        cd $InDir/$Sample/Mandalorion/
        echo "Starting Mandalorion"

	/usr/bin/time -v python3 ~/apps/Mandalorion/Mando.py -t $threads -p ./ -g $reference_gtf -G $reference_fa -f $InDir/${Sample}.fastq

	echo "Mandalorion terminates"
}

RATTLE(){
	mkdir -p $InDir/$Sample/rattle/
	cd $InDir/$Sample/rattle/
	echo "Starting RATTLE"

	#Filter fastq to have only reads woth length > 150 (bad alloc error otherwise)
	awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>150{print a"\n"b"\n"c"\n"$0;}' $InDir/${Sample}.fastq > $InDir/${Sample}_sup150.fastq

	if [[ $Sample =~ RNA ]] ; then
        	echo "RNA mode ON"
		rattle cluster -i $InDir/${Sample}_sup150.fastq -o ./ -t $threads --iso --rna
	else
		echo "RNA mode OFF"
		rattle cluster -i InDir/${Sample}_sup150.fastq -o ./ -t $threads --iso
	fi
        mkdir -p clusters
        rattle extract_clusters -i $InDir/${Sample}_sup150.fastq -c ./clusters.out -o ./clusters --fastq
        rattle correct -i $InDir/${Sample}_sup150.fastq -c ./clusters.out -o ./ -t $threads
        rattle polish -i ./consensi.fq -o ./ -t $threads

	echo "RATTLE terminates"
}

Isonclust(){
	mkdir -p $InDir/$Sample/isONclust/
	cd $InDir/$Sample/isONclust/
	echo "Starting isONclust"

	isONclust --ont --fastq $InDir/${Sample}.fastq --outfolder ./ --t $threads
        #making fastq file
        isONclust write_fastq --clusters ./final_clusters.tsv --fastq $InDir/${Sample}.fastq --outfolder ./output_fastq --N 1
        #Convert tsv output into a fastq
        awk -F"\t" '{print "@cluster"$1" "$2"\t"$3"\t+\t"$4}' final_cluster_origins.tsv | sed 's/\t/\n/g' - > final_cluster_origins.fastq

	echo "isONclust terminates"
}


Isonclust2(){
	mkdir -p $InDir/$Sample/isONclust2/
	mkdir -p $InDir/$Sample/isONclust2/batches/
	mkdir -p $InDir/$Sample/isONclust2/cluster/
	mkdir -p $InDir/$Sample/isONclust2/results/
	cd $InDir/$Sample/isONclust2/

	if [[ $Sample =~ RNA ]]; then
		echo "RNA set"
		#isONclust2 can't take U in the input fastq ; conversion into T is necessary
		if ! [ -f $InDir/${Sample}_Tconv.fastq ]; then

			sed '/^[^>]/s/U/T/g' $InDir/${Sample}.fastq > $InDir/${Sample}_Tconv.fastq
		fi
		inputfastq=$InDir/${Sample}_Tconv.fastq
	else
		inputfastq=$InDir/${Sample}.fastq
	fi

	#Create batches (50000kb per batch)
	isONclust2 sort -v -o ./batches $inputfastq
	#Initial clustering
	for file in $InDir/$Sample/isONclust2/batches/batches/*.cer ; do
		FILE_NAME=${file##*/}
		isONclust2 cluster -v -x sahlin -l $file -o $InDir/$Sample/isONclust2/cluster/$FILE_NAME
	done
	#Cluster again by merging initial clusters 2 by 2
	first=($InDir/$Sample/isONclust2/cluster/*.cer)
	first_file=${first[0]}
	i=0
	for file in $InDir/$Sample/isONclust2/cluster/*.cer ; do
		if [[ $file != ${first[0]} ]] ; then
			next=$((i+1))
			isONclust2 cluster -v -x sahlin -l $first_file -r $file -o $InDir/$Sample/isONclust2/cluster/z${i}_${next}.cer
			first_file=$InDir/$Sample/isONclust2/cluster/z${i}_${next}.cer
			i=$next
		fi
	done
	#Dump final results
	var=$(ls $InDir/$Sample/isONclust2/cluster/ | sort -V | tail -n 1)
	echo $var
	isONclust2 dump -v -i $InDir/$Sample/isONclust2/batches/sorted_reads_idx.cer -o $InDir/$Sample/isONclust2/results $InDir/$Sample/isONclust2/cluster/$var
	echo "isONclust2 terminates"
}

RNAbloom(){
	mkdir -p $InDir/$Sample/RNAbloom/
        cd $InDir/$Sample/RNAbloom/
        echo "Starting RNAbloom"

	java -jar ~/apps/RNA-Bloom_v1.4.3/RNA-Bloom.jar -long $InDir/${Sample}.fastq -fpr 0.01 -t $threads -outdir ./

	echo "RNAbloom terminates"
}

RNAbloom2(){
	mkdir -p $InDir/$Sample/RNAbloom2/
        cd $InDir/$Sample/RNAbloom2/
        echo "Starting RNA bloom 2"
        rnabloom -long $InDir/${Sample}.fastq -fpr 0.01 -t $threads -outdir ./

	echo "RNABloom2 terminates"
}

GFF(){
	cd $InDir/$Sample/Final-Assemblies
	mkdir -p $InDir/$Sample/Final-Assemblies/processed
	echo "GFFcompare starts"
	for file in *.f* ; do
		~/apps/minimap2/minimap2-2.24/minimap2 -t $threads -ax splice --secondary=no $reference_fa $file | samtools view -@ 6 -b - | samtools sort -@ 6 - > ${file%*.f*}.bam
        	spliced_bam2gff -M ${file%*.f*}.bam > ${file%*.f*}_converted.gff
        	sed 's/;0/-0/g' ${file%*.f*}_converted.gff | sed 's/;16/-16/g' > ${file%*.f*}_converted-corrected.gff
		rm ${file%*.f*}.bam
		rm ${file%*.f*}_converted.gff
	done
	for file in *.g* ; do
		echo $file
		~/apps/gffcompare/gffcompare -R -r $reference_gtf -o processed/$file $file
		cp ${file}.${file}.tmap ../TMAPs/${file%*.g*}.tmap
	done
	mv *map processed/
	mv *corrected* processed/
	cd ../TMAPs
	for file in *_converted-corrected.tmap ; do
		mv $file ${file%*_converted-corrected.tmap}.tmap
	done
	echo "GFFcompare terminates"
}

SQ3(){
	cd $InDir/$Sample/Final-Assemblies
	echo $InDir/$Sample/Final-Assemblies
        mkdir -p $InDir/$Sample/Final-Assemblies/processed
	mkdir -p $InDir/$Sample/SQ3
        echo "SQANTI3 starts"

	export PYTHONPATH=$PYTHONPATH:~/apps/SQANTI3-5.2/cDNA_Cupcake/sequence/
	export PYTHONPATH=$PYTHONPATH:~/apps/SQANTI3-5.2/cDNA_Cupcake/

	for file in *.g* ; do
		echo $file
		python ~/apps/SQANTI3-5.2/sqanti3_qc.py --force_id_ignore --aligner_choice minimap2 -d processed/ --skipORF --report html $file $reference_gtf  $reference_fa
	done
	for file in *.f* ; do
		echo $file
		python ~/apps/SQANTI3-5.2/sqanti3_qc.py --force_id_ignore --aligner_choice minimap2 -d processed/ --skipORF --fasta --report html $file $reference_gtf $reference_fa
	done
	cp processed/*classification* ../SQ3/
	rm *renamed*
	echo "SQANTI3 terminates"
}


TPFPFN(){
	mkdir -p $InDir/$Sample/TPFP
	cd $InDir/$Sample/TMAPs
	if [[ -f "$InDir/$Sample/GFF-PrecisionSensitivity-update.csv" ]]; then rm $InDir/$Sample/GFF-PrecisionSensitivity.csv; fi
	subset1=${Sample#*b}
	if [[ $subset1 =~ k ]] ; then dividor=1000 ; else dividor=1 ; fi
	subset=${subset1%*k}
	LOD=$(awk "BEGIN {print $subset / $dividor}")

	if [[ $Sample =~ SIRV ]] ; then
		LODsirv=$(awk "BEGIN {print 0.108555 / $LOD}")
        	awk -v LOD="$LODsirv" '$2 <= LOD {print $1}' ~/references/SIRV-list_abund.tsv > $InDir/$Sample/TPFP/Ttoremove.tsv
        	awk -v LOD="$LODsirv" '$2 > LOD {print $1}' ~/references/SIRV-list_abund.tsv > $InDir/$Sample/TPFP/Ttokeep.tsv
	else
		LODsequin=$(awk "BEGIN {print 0.0151965 / $LOD}")
		awk -v LOD="$LODsequin" '$3 <= LOD {print $1}' ~/references/rnasequin_isoforms_2.5.tsv > $InDir/$Sample/TPFP/Ttoremove.tsv
		awk -v LOD="$LODsequin" '$3 > LOD {print $1}' ~/references/rnasequin_isoforms_2.5.tsv | grep -v "NAME" > $InDir/$Sample/TPFP/Ttokeep.tsv
	fi

	for file in *.tmap ; do
		TP=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttokeep.tsv ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l )
		FP1=$( join -1 2 -2 1 -v1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttoremove.tsv ) | awk '$3!="=" && $3!="c" && $3!="u" {print $1}' | grep -v "ref_id" | wc -l )
		FP2=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
		echo $file,$TP,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $InDir/$Sample/TPFP/Ttokeep.tsv ) | wc -l ) | sed 's/.tmap//g' >> $InDir/$Sample/GFF-PrecisionSensitivity.csv
	done

	cd $InDir/$Sample/SQ3
	if [[ -f "$InDir/$Sample/SQ3-PrecisionSensitivity-update.csv" ]]; then rm $InDir/$Sample/SQ3-PrecisionSensitivity.csv; fi
	for file in *_classification.txt ; do
		TP=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttokeep.tsv ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $14!="intron_retention" ) {print $1}' | sort -u | wc -l )
		FP1=$( join -1 8 -2 1 -v1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttoremove.tsv ) | awk '$7!="full-splice_match" && ( $7!="incomplete-splice_match" || $15=="intron_retention" ) && $7!="intergenic" {print $1}' | grep -v "structural_category" | wc -l )
		FP2=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
		echo $file,$TP,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $InDir/$Sample/TPFP/Ttokeep.tsv ) | wc -l ) | sed 's/_classification.txt//g' >> $InDir/$Sample/SQ3-PrecisionSensitivity.csv
	done
}

## executing functions in order:
main(){
	Sample=PCS109_chrIS_mixA_cDNA_sub150k #Fill in sample name
	InDir=/opt/benchmarking/data/sub150k #Fill directory with input .fastq files
	molType=DNA  #<DNA> or <RNA>
	threads=8

	if [[ $Sample =~ SIRV ]]; then
		reference_gtf=~/references/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf
		reference_gtf_talon=~/references/SIRV_ERCC_longSIRV_multi-fasta_20210507_reformatted.gtf
		reference_fa=~/references/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta
		chr="SIRV"
	else
		reference_gtf=~/references/rnasequin_annotation_2.4.gtf
        	reference_gtf_talon=~/references/rnasequin_annotation_2.4.gtf
        	reference_fa=~/references/chrIS.fa
        	chr="chrIS"
	fi

        export Sample=$Sample
	export mainInputDir=$InDir
	export reference_gtf=$reference_gtf
	export reference_gtf_talon=$reference_gtf_talon
        export reference_fa=$reference_fa
	export chr=$chr
	export moltype=$moltype
	export threads=$threads

	mkdir -p $Sample
	mkdir -p ${Sample}/Final-Assemblies
	mkdir -p ${Sample}/TMAPs

	preprocess

	baseline
	cp $InDir/$Sample/control/${Sample}_converted-corrected.gff $InDir/$Sample/Final-Assemblies/Control.gff

#Guided + De novo
	stringtie2
	cp $InDir/$Sample/stringtie2/${Sample}_strg2def.gtf $InDir/$Sample/Final-Assemblies/Stringtie2.gtf
	cp $InDir/$Sample/stringtie2/${Sample}_strg2def-noRef.gtf $InDir/$Sample/Final-Assemblies/Stringtie2-noRef.gtf

	source /home/sagmel/mambaforge/bin/activate flair #conda activate flair
	Flair
	cp $InDir/$Sample/flair/${Sample}_flair.isoforms.gtf $InDir/$Sample/Final-Assemblies/FLAIR.gtf
	cp $InDir/$Sample/flair-noRef/${Sample}_flair-noRef.isoforms.fa $InDir/$Sample/Final-Assemblies/FLAIR-noRef.fa
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

	source /home/sagmel/mambaforge/bin/activate isoquant #conda activate isoquant
	isoQuant
	cp $InDir/$Sample/isoQuant/OUT/OUT.transcript_models.gtf $InDir/$Sample/Final-Assemblies/isoQuant.gtf

	isoQuant-noRef
	cp $InDir/$Sample/isoQuant-noRef/OUT/OUT.transcript_models.gtf $InDir/$Sample/Final-Assemblies/isoQuant-noRef.gtf
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

#Guided only
	source /home/sagmel/mambaforge/bin/activate TALON5 #conda activate TALON5
	TALON
	cp $InDir/$Sample/talon/${Sample}_gtf_talon.gtf $InDir/$Sample/Final-Assemblies/TALON.gtf

	TALON_reco
	cp $InDir/$Sample/talon_reco/${Sample}_gtf_talon.gtf $InDir/$Sample/Final-Assemblies/TALON_reco.gtf
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

	source /home/sagmel/mambaforge/bin/activate FLAMES #conda activate FLAMES
	FLAMES
	cp $InDir/$Sample/FLAMES/isoform_annotated.filtered_True.gff3 $InDir/$Sample/Final-Assemblies/FLAMES.gff3
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

	Mandalorion
	cp $InDir/$Sample/Mandalorion/Isoforms.filtered.clean.gtf $InDir/$Sample/Final-Assemblies/Mandalorion.gtf

#Ab initio
	#run <bulker activate seqtools/seqtools> before running RATTLE
	RATTLE
	cp $InDir/$Sample/rattle/transcriptome.fq $InDir/$Sample/Final-Assemblies/RATTLE.fastq
	#close bulker before running anything else <exit>

	source /home/sagmel/mambaforge/bin/activate isonclust #conda activate isonclust
	Isonclust
	cp $InDir/$Sample/isONclust/final_cluster_origins.fastq $InDir/$Sample/Final-Assemblies/isONclust.fastq
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

	source /home/sagmel/mambaforge/bin/activate isonform #conda activate isonform
	IsonPipeline
	cp $InDir/$Sample/isONpipeline/isoforms/transcriptome.fasta $InDir/$Sample/Final-Assemblies/isONpipeline.fa
	source /home/sagmel/mambaforge/bin/deactivate #conda deactivate

	Isonclust2
	cp $InDir/$Sample/isONclust2/results/cluster_cons.fq $InDir/$Sample/Final-Assemblies/isONclust2.fastq

	RNAbloom
	cp $InDir/$Sample/RNAbloom/rnabloom.transcripts.fa $InDir/$Sample/Final-Assemblies/RNAbloom.fa

	RNAbloom2
	cp $InDir/$Sample/RNAbloom2/rnabloom.transcripts.fa $InDir/$Sample/Final-Assemblies/RNAbloom2.fa


#	bambu   #Doesn't work within the script --> run manually on the side

#Post-process
	GFF
	source /home/sagmel/mambaforge/bin/activate SQANTI3.env
	SQ3
	source /home/sagmel/mambaforge/bin/deactivate

	TPFPFN

	echo "All Done"
}

main
