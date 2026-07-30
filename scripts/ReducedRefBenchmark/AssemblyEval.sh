#!/bin/bash
#Script to run compute TP;FP;FN on discovery assessment datasets
#Fill in the 4 lines of information in the main

LRAA_handling(){
                cd $InDir/$Subset/$Sample/Final-Assemblies
                mkdir -p $InDir/$Subset/$Sample/Final-Assemblies/LRAA_OG_gtfs/
                echo $InDir/$Subset/$Sample/Final-Assemblies
                echo "awk changing LRAA gtfs"
                for lraa in LRAA_ref.gtf LRAA_denovo.gtf; do
                        if [[ -f "$lraa" ]]; then
                                cleaned="${lraa%.gtf}_AWKed.gtf"
                                echo "  Cleaning $lraa → $cleaned"

                                awk '
                                BEGIN{ OFS="\t" }

                                {
                                        attr="";
                                        for (i=9; i<=NF; i++) attr = attr $i " ";

                                        gsub(/" +/, "\"; ", attr);
                                        gsub(/; ;/, ";", attr);

                                        if (attr !~ /;$/) attr = attr ";";
                                        gsub(/;;+/, ";", attr);

                                        sub(/^ +/, "", attr);
                                        sub(/ +$/, "", attr);

                                        # remove final trailing semicolon
                                        sub(/;$/, "", attr);

                                        $9 = attr;
                                        NF = 9;

                                        print
                                }' "$lraa" > "$cleaned"
                        fi
                        mv $InDir/$Subset/$Sample/Final-Assemblies/$lraa $InDir/$Subset/$Sample/Final-Assemblies/LRAA_OG_gtfs/
                done
                echo "awk changing LRAA gtfs completed"
}


GFF(){
        cd $InDir/$Subset/$Sample/Final-Assemblies
        mkdir -p $InDir/$Subset/$Sample/Final-Assemblies/processed
		mkdir -p $InDir/$Subset/$Sample/TMAPs
        echo "GFFcompare starts"
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
        cd $InDir/$Subset/$Sample/Final-Assemblies
        echo $InDir/$Subset/$Sample/Final-Assemblies
        mkdir -p $InDir/$Subset/$Sample/Final-Assemblies/processed
        mkdir -p $InDir/$Subset/$Sample/SQ3
        echo "SQANTI3 starts"

        export PYTHONPATH=$PYTHONPATH:~/apps/SQANTI3-5.2/cDNA_Cupcake/sequence/
        export PYTHONPATH=$PYTHONPATH:~/apps/SQANTI3-5.2/cDNA_Cupcake/

        for file in *.g* ; do
                echo $file
                python ~/apps/SQANTI3-5.2/sqanti3_qc.py --force_id_ignore --aligner_choice minimap2 -d processed/ --skipORF --report html $file $reference_gtf $reference_fa
        done
        cp processed/*classification* ../SQ3/
        rm *renamed*
        echo "SQANTI3 terminates"
}

SplitSections(){
        cd $InDir/$Subset/$Sample/Final-Assemblies
        for file in *.gtf ; do
                if [[ $file == "Mandalorion.gtf" ]] ; then
                        paste <( awk -F"\t" '$3=="transcript" {print $9}' $InDir/$Subset/$Sample/Final-Assemblies/$file | awk -F"; " '{print $1}' | sed 's/transcript_id "//g' | sed 's/ //g' | sed 's/"//g' | sed 's/;//g' ) <( awk -F"\t" '$3=="transcript" {print $4"\t"$5}' $InDir/$Subset/$Sample/Final-Assemblies/$file ) > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst
                        awk -F"\t" '$3<=3522628 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section1.lst
                        awk -F"\t" '$2>=3522628 && $3<=7045256 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section2.lst
                        awk -F"\t" '$2>=7045256 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section3.lst
                else
                        paste <( awk -F"\t" '$3=="transcript" {print $9}' $InDir/$Subset/$Sample/Final-Assemblies/$file | awk -F"; " '{print $2}' | sed 's/transcript_id "//g' | sed 's/ //g' | sed 's/"//g' | sed 's/;//g' ) <( awk -F"\t" '$3=="transcript" {print $4"\t"$5}' $InDir/$Subset/$Sample/Final-Assemblies/$file ) > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst
                        awk -F"\t" '$3<=3522628 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section1.lst
                        awk -F"\t" '$2>=3522628 && $3<=7045256 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section2.lst
                        awk -F"\t" '$2>=7045256 {print $1}' $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}.lst > $InDir/$Subset/$Sample/Final-Assemblies/${file%*.gtf}_section3.lst
                fi
        done

        cd $InDir/$Subset/$Sample/TMAPs
        for file in *.tmap ; do
                join -1 1 -2 5 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*.tmap}_section1.lst ) <( sort -k5,5 $file ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > ${file%*.tmap}-${section1}.tmap
                join -1 1 -2 5 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*.tmap}_section2.lst ) <( sort -k5,5 $file ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > ${file%*.tmap}-${section2}.tmap
                join -1 1 -2 5 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*.tmap}_section3.lst ) <( sort -k5,5 $file ) | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' > ${file%*.tmap}-${section3}.tmap
        done

        cd $InDir/$Subset/$Sample/SQ3
        for file in *classification.txt ; do
                join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*_classification.txt}_section1.lst ) <( sort -k1,1 $file ) | sed 's/ /\t/g' > ${file%*.txt}-${section1}.txt
                join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*_classification.txt}_section2.lst ) <( sort -k1,1 $file ) | sed 's/ /\t/g' > ${file%*.txt}-${section2}.txt
                join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/Final-Assemblies/${file%*_classification.txt}_section3.lst ) <( sort -k1,1 $file ) | sed 's/ /\t/g' > ${file%*.txt}-${section3}.txt
        done
}

TPFPFN(){
        rm -f $InDir/$Subset/$Sample/GFF-TPnovelFPFN.csv
        cd $InDir/$Subset/$Sample/TMAPs
        for file in *full.tmap ; do echo $Subset,$file,"full",$FullSection,$( awk '$3=="=" || $3=="c" {print $2}' $file | sort -u | wc -l ),0,$( awk '$3!="=" && $3!="c" && $3!="u" {print $2}' $file | grep -v "ref_id" | wc -l ),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $FullRef ) | grep -v "NAME" | wc -l ) | sed 's/-full.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN.csv; done
        for file in *null.tmap ; do echo $Subset,$file,"null",$NullSection,0,$( awk '$3=="=" || $3=="c" {print $2}' $file | sort -u | wc -l ),$( awk '$3!="=" && $3!="c" && $3!="u" {print $2}' $file | grep -v "ref_id" | wc -l ),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $NullRef ) | grep -v "NAME" | wc -l ) | sed 's/-null.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN.csv; done
        for file in *partial.tmap ; do echo $Subset,$file,"partial",$PartialSection,$( join -1 2 -2 1 <( sort -k2,2 $file) <( sort -k1,1 $PartialRefGiven ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l ),$( join -1 2 -2 1 <( sort -k2,2 $file) <( sort -k1,1 $PartialRefnotGiven ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l ),$( awk '$3!="=" && $3!="c" && $3!="u" {print $2}' $file | grep -v "ref_id" | wc -l ),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( cat $PartialRefGiven $PartialRefnotGiven | sort -k1,1 ) | grep -v "NAME" | wc -l ) | sed 's/-partial.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN.csv; done

        rm -f $InDir/$Subset/$Sample/SQ3-TPnovelFPFN.csv
        cd $InDir/$Subset/$Sample/SQ3
        for file in *_classification-full.txt ; do echo $Subset,$file,"full",$FullSection,$( awk '$6=="full-splice_match" || ( $6=="incomplete-splice_match" && $14!="intron_retention" ) {print $8}' $file | sort -u | wc -l ),0,$( awk '$6!="full-splice_match" && ( $6!="incomplete-splice_match" || $15=="intron_retention" ) && $6!="intergenic" {print $8}' $file | grep -v "associated_transcript" | wc -l ),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $FullRef ) | grep -v "NAME" | wc -l ) | sed 's/_classification-full.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN.csv; done
        for file in *_classification-null.txt ; do echo $Subset,$file,"null",$NullSection,0,$( awk '$6=="full-splice_match" || ( $6=="incomplete-splice_match" && $14!="intron_retention" ) {print $8}' $file | sort -u | wc -l ),$( awk '$6!="full-splice_match" && ( $6!="incomplete-splice_match" || $15=="intron_retention" ) && $6!="intergenic" {print $8}' $file | grep -v "associated_transcript" | wc -l ),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $NullRef ) | grep -v "NAME" | wc -l ) | sed 's/_classification-null.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN.csv; done
        for file in *_classification-partial.txt ; do echo $Subset,$file,"partial",$PartialSection,$( join -1 8 -2 1 <( sort -k8,8 $file) <( sort -k1,1 $PartialRefGiven ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $15!="intron_retention" ) {print $1}' | sort -u | wc -l ),$( join -1 8 -2 1 <( sort -k8,8 $file) <( sort -k1,1 $PartialRefnotGiven ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $15!="intron_retention" ) {print $1}' | sort -u | wc -l ),$( awk '$6!="full-splice_match" && ( $6!="incomplete-splice_match" || $15=="intron_retention" ) && $6!="intergenic" {print $8}' $file | grep -v "associated_transcript" | wc -l ),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( cat $PartialRefGiven $PartialRefnotGiven | sort -k1,1 ) | grep -v "NAME" | wc -l ) | sed 's/_classification-partial.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN.csv; done
}
TPFPupdate(){
        mkdir -p $InDir/$Subset/$Sample/TPFP
        cd $InDir/$Subset/$Sample/TMAPs
        rm -f $InDir/$Subset/$Sample/GFF-TPnovelFPFN-update.csv
        subset1=${Sample#*b}
        if [[ $subset1 =~ k ]] ; then dividor=1000 ; else dividor=1 ; fi
        subset2=${subset1%*k}
        LOD=$(awk "BEGIN {print $subset2 / $dividor}")

        LODsequin=$(awk "BEGIN {print 0.0151965 / $LOD}")
        awk -v LOD="$LODsequin" '$3 <= LOD {print $1}' /opt/benchmarking/data/references/rnasequin_isoforms_2.5.tsv > $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv
        awk -v LOD="$LODsequin" '$3 > LOD {print $1}' /opt/benchmarking/data/references/rnasequin_isoforms_2.5.tsv | grep -v "NAME" > $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv

        join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) <( sort -k1,1 $FullRef ) | sed 's/ /\t/g' > $InDir/$Subset/$Sample/TPFP/TtokeepFullRef.tsv
        join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) <( sort -k1,1 $NullRef ) | sed 's/ /\t/g' > $InDir/$Subset/$Sample/TPFP/TtokeepNullRef.tsv
        join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) <( cat $PartialRefGiven $PartialRefnotGiven | sort -k1,1 ) | sed 's/ /\t/g' > $InDir/$Subset/$Sample/TPFP/TtokeepPartialRef.tsv
        join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) <( sort -k1,1 $PartialRefGiven ) | sed 's/ /\t/g' > $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefGiven.tsv
        join -1 1 -2 1 <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) <( sort -k1,1 $PartialRefnotGiven ) | sed 's/ /\t/g' > $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefnotGiven.tsv

        for file in *full.tmap ; do
                TP=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l )
                FP1=$( join -1 2 -2 1 -v1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$3!="=" && $3!="c" && $3!="u" {print $1}' | grep -v "ref_id" | wc -l )
                FP2=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"full",$FullSection,$TP,0,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepFullRef.tsv ) | wc -l ) | sed 's/-full.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN-update.csv
        done

        for file in *null.tmap ; do
                novel=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l )
                FP1=$( join -1 2 -2 1 -v1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$3!="=" && $3!="c" && $3!="u" {print $1}' | grep -v "ref_id" | wc -l )
                FP2=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"null",$NullSection,0,$novel,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepNullRef.tsv ) | wc -l ) | sed 's/-null.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN-update.csv
        done

        for file in *partial.tmap ; do
                TP=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefGiven.tsv ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l )
                novel=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefnotGiven.tsv ) | awk '$3=="=" || $3=="c" {print $1}' | sort -u | wc -l )
                FP1=$( join -1 2 -2 1 -v1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$3!="=" && $3!="c" && $3!="u" {print $1}' | grep -v "ref_id" | wc -l )
                FP2=$( join -1 2 -2 1 <( sort -k2,2 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"partial",$PartialSection,$TP,$novel,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 2 -2 1 -v 2  <( sort -k2,2 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRef.tsv ) | wc -l ) | sed 's/-partial.tmap//g' >> $InDir/$Subset/$Sample/GFF-TPnovelFPFN-update.csv
        done

        cd $InDir/$Subset/$Sample/SQ3
        rm -f $InDir/$Subset/$Sample/SQ3-TPnovelFPFN-update.csv
        for file in *_classification-full.txt ; do
                TP=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $14!="intron_retention" ) {print $1}' | sort -u | wc -l )
                FP1=$( join -1 8 -2 1 -v1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$7!="full-splice_match" && ( $7!="incomplete-splice_match" || $15=="intron_retention" ) && $7!="intergenic" {print $1}' | grep -v "associated_transcript" | wc -l )
                FP2=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"full",$FullSection,$TP,0,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepFullRef.tsv ) | wc -l ) | sed 's/_classification-full.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN-update.csv
        done

        for file in *_classification-null.txt ; do
                novel=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttokeep.tsv ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $14!="intron_retention" ) {print $1}' | sort -u | wc -l )
                FP1=$( join -1 8 -2 1 -v1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$7!="full-splice_match" && ( $7!="incomplete-splice_match" || $15=="intron_retention" ) && $7!="intergenic" {print $1}' | grep -v "associated_transcript" | wc -l )
                FP2=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"null",$NullSection,0,$novel,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepNullRef.tsv ) | wc -l ) | sed 's/_classification-null.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN-update.csv
        done

        for file in *_classification-partial.txt ; do
                TP=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefGiven.tsv ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $15!="intron_retention" ) {print $1}' | sort -u | wc -l )
                novel=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRefnotGiven.tsv ) | awk '$7=="full-splice_match" || ( $7=="incomplete-splice_match" && $15!="intron_retention" ) {print $1}' | sort -u | wc -l )
                FP1=$( join -1 8 -2 1 -v1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | awk '$7!="full-splice_match" && ( $7!="incomplete-splice_match" || $15=="intron_retention" ) && $7!="intergenic" {print $1}' | grep -v "associated_transcript" | wc -l )
                FP2=$( join -1 8 -2 1 <( sort -k8,8 $file ) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/Ttoremove.tsv ) | wc -l )
                echo $Subset,$file,"partial",$PartialSection,$TP,$novel,$(awk "BEGIN {print $FP1 + $FP2}"),$( join -1 8 -2 1 -v 2  <( sort -k8,8 $file) <( sort -k1,1 $InDir/$Subset/$Sample/TPFP/TtokeepPartialRef.tsv ) | wc -l ) | sed 's/_classification-partial.txt//g' >> $InDir/$Subset/$Sample/SQ3-TPnovelFPFN-update.csv
        done
}

## executing functions in order:
main(){

        while read GTF_file
        do
                NOW_T="$(date +'%T')"
                Subset=chrIS_${GTF_file}
                Sample=LSK109_chrIS_mixA_cDNA_sub150k #Fill in sample name
                InDir=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/Outputs
                molType=DNA  #<DNA> or <RNA>
                threads=8

                if [[ $Sample =~ SIRV ]]; then
                        reference_gtf=/opt/benchmarking/data/references/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf
                        reference_gtf_talon=/opt/benchmarking/data/references/SIRV_ERCC_longSIRV_multi-fasta_20210507_reformatted.gtf
                        reference_fa=/opt/benchmarking/data/references/SIRV_ERCC_longSIRV_multi-fasta_20210507.fasta
                        chr="SIRV"
                else
                        reference_gtf=/opt/benchmarking/data/references/rnasequin_annotation_2.4.gtf
                        reference_gtf_talon=/opt/benchmarking/data/references/rnasequin_annotation_2.4.gtf
                        reference_fa=/opt/benchmarking/data/references/chrIS.fa
                        chr="chrIS"
                fi

                export Sample=$Sample
                export InDir=$InDir
                export reference_gtf=$reference_gtf
                export reference_gtf_talon=$reference_gtf_talon
                export reference_fa=$reference_fa
                export chr=$chr
                export moltype=$moltype
                export threads=$threads


                if [[ $Subset == "chrIS_full_half_null" ]] ; then
                        FullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section1_full.lst
                        FullSection=Section1
                        section1=full
                        PartialRefGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section2_partial.lst
                        PartialRefnotGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section2_partial-notgiven.lst
                        PartialSection=Section2
                        section2=partial
                        NullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section3_null.lst
                        NullSection=Section3
                        section3=null
                elif [[ $Subset == "chrIS_half_null_full" ]] ; then
                        FullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section3_full.lst
                        FullSection=Section3
                        section3=full
                        PartialRefGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section1_partial.lst
                        PartialRefnotGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section1_partial-notgiven.lst
                        PartialSection=Section1
                        section1=partial
                        NullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section2_null.lst
                        NullSection=Section2
                        section2=null
                elif [[ $Subset == "chrIS_null_full_half" ]] ; then
                        FullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section2_full.lst
                        FullSection=Section2
                        section2=full
                        PartialRefGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section3_partial.lst
                        PartialRefnotGiven=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section3_partial-notgiven.lst
                        PartialSection=Section3
                        section3=partial
                        NullRef=/opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/sections_gtfs/section1_null.lst
                        NullSection=Section1
                        section1=null
                fi

                export Subset=$Subset
                export FullRef=$FullRef
                export FullSection=$FullSection
                export section3=$section3
                export PartialRefGiven=$PartialRefGiven
                export PartialRefnotGiven=$PartialRefnotGiven
                export PartialSection=$PartialSection
                export section2=$section2
                export NullRef=$NullRef
                export NullSection=$NullSection
                export section1=$section1

                echo "##########################################################"
                echo "GTF:"
                echo $GTF_file
                echo "Sample:$Sample"
                echo "time: ${NOW_T}"
				
				GFF
				LRAA_handling
				source /home/sagmel/mambaforge/bin/activate SQANTI3.env
				SQ3
				source /home/sagmel/mambaforge/bin/deactivate
				echo "GFF and SQ3 done"
				
                SplitSections
                echo "Split Sections function done"
                echo "time: ${NOW_T}"
                TPFPFN
                echo "TPFPFN function done"
                echo "time: ${NOW_T}"
                TPFPupdate
                echo "TPFPupdate function done"
                echo "time: ${NOW_T}"

        done < /opt/benchmarking/data/20251118_Benchmarking/novel_discovery_revision4/Outputs/gtfs_list.txt
        echo "All Done"
        echo "time: ${NOW_T}"

}

main