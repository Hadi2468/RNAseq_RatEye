#!/bin/bash
# This is the final script for RSEM analysis (42 samples with one run, and 11 samples with two runs)

clear 
sample_number=1
for file_run1_forward in /home/h/Desktop/RSEM/data/Rat_HS_eyes/C202SC19070808_1/*_1.fq.gz
	do
	filename=${file_run1_forward##*/}
	echo; echo "********** Sample_Number "$sample_number": ${filename:0:12}"; echo "**********"
	date
	file_run1_reverse=`echo $file_run1_forward |sed "s/_1.fq.gz/_2.fq.gz/"`
	echo "File_Run1_Forward:" $file_run1_forward
	echo "File_Run1_Reverse:" $file_run1_reverse
	if [ -e /home/h/Desktop/RSEM/data/Rat_HS_eyes/C202SC19070808_2/$filename ]; then
		file_run2_forward=/home/h/Desktop/RSEM/data/Rat_HS_eyes/C202SC19070808_2/$filename
		file_run2_reverse=`echo $file_run2_forward |sed "s/_1.fq.gz/_2.fq.gz/"`
		echo "File_Run2_Forward:" $file_run2_forward
		echo "File_Run2_Reverse:" $file_run2_reverse
		file_two_runs_forward="/home/h/Desktop/RSEM/data/Rat_HS_eyes/data_merged/tworuns_"$filename
		file_two_runs_reverse=`echo $file_two_runs_forward |sed "s/_1.fq.gz/_2.fq.gz/"`
		# Combine two data for the same rat obtained in two different runs
		cat $file_run1_forward $file_run2_forward > $file_two_runs_forward
		cat $file_run1_reverse $file_run2_reverse > $file_two_runs_reverse
		echo "File_TwoRuns_Forward:" $file_two_runs_forward
		echo "File_TwoRuns_Reverse:" $file_two_runs_reverse
		software/RSEM-1.3/rsem-calculate-expression -p 12 --paired-end --bowtie2 --bowtie2-path software/bowtie2-2.2 --estimate-rspd --append-names --sort-bam-by-coordinate --output-genome-bam $file_two_runs_forward $file_two_runs_reverse ref/rat_ref export/${filename:0:12}
	else
		software/RSEM-1.3/rsem-calculate-expression -p 12 --paired-end --bowtie2 --bowtie2-path software/bowtie2-2.2 --estimate-rspd --append-names --sort-bam-by-coordinate --output-genome-bam $file_run1_forward $file_run1_reverse ref/rat_ref export/${filename:0:12}
	date
	fi
	((sample_number++))
done

