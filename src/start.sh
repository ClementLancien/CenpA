#!/usr/bin/env bash
#==> MAKEFILE ??????

#exit as soon as any line in the script fails
set -e 

#exit as soon as any variable is not init
set -o nounset #equals set -u

. /data/tmp/clancien/Script/all.sh

#__get_bed_files /data/tmp/clancien/Aligned/CENH3_OVEREXP.1000.csv SRR633614 SRR633612
#__get_bed_files /data/tmp/clancien/Aligned/CENH3_OVEREXP.1000.csv SRR633615 SRR633613

# __overlap /data/tmp/clancien/Result/Bed/SRR633612.bed /data/tmp/clancien/Aligned/SRR633612.sorted.bam
# __overlap /data/tmp/clancien/Result/Bed/SRR633613.bed /data/tmp/clancien/Aligned/SRR633613.sorted.bam
# __overlap /data/tmp/clancien/Result/Bed/SRR633614.bed /data/tmp/clancien/Aligned/SRR633614.sorted.bam
# __overlap /data/tmp/clancien/Result/Bed/SRR633615.bed /data/tmp/clancien/Aligned/SRR633615.sorted.bam

# __sort_by_name $OVERLAP/SRR633612.overlap.bam
# __sort_by_name $OVERLAP/SRR633613.overlap.bam
# __sort_by_name $OVERLAP/SRR633614.overlap.bam
# __sort_by_name $OVERLAP/SRR633615.overlap.bam

# __bam_to_fastq $BAM/SRR633612.overlap.sorted.bam
# __bam_to_fastq $BAM/SRR633613.overlap.sorted.bam
# __bam_to_fastq $BAM/SRR633614.overlap.sorted.bam
# __bam_to_fastq $BAM/SRR633615.overlap.sorted.bam

# __count_read $FASTQ/SRR633612.R1.fq
# __count_read $FASTQ/SRR633613.R1.fq
# __count_read $FASTQ/SRR633614.R1.fq
# __count_read $FASTQ/SRR633615.R1.fq

#__merge_count_read $FASTQ/countread.txt $FASTQ/SRR633612.countread.txt $FASTQ/SRR633613.countread.txt $FASTQ/SRR633614.countread.txt $FASTQ/SRR633615.countread.txt

#__downsample_auto $FASTQ/SRR633614.R1.fq $FASTQ/SRR633614.R2.fq $FASTQ/SRR633612.R1.fq $FASTQ/SRR633612.R2.fq $FASTQ/countread.txt 10 3
#__downsample_auto $FASTQ/SRR633615.R1.fq $FASTQ/SRR633615.R2.fq $FASTQ/SRR633613.R1.fq $FASTQ/SRR633613.R2.fq $FASTQ/countread.txt 10 3

# __align_all SRR633614 SRR633612 10 3
# __align_all SRR633615 SRR633613 10 3

#__overlap_cytoband SRR633614 SRR633612 10 3
#__overlap_cytoband SRR633615 SRR633613 10 3

#__merge_overlap_cytoband SRR633614 SRR633612 10 3 1
#__merge_overlap_cytoband SRR633615 SRR633613 10 2 0

# __overlap_features SRR633614 SRR633612 10 3
# __overlap_features SRR633615 SRR633613 10 3

#__merge_overlap_features SRR633614 SRR633612 10 3 1
#__merge_overlap_features SRR633615 SRR633613 10 3 0

__show_options()
{
    printf "$formatPrintBlue" "Usage: $0"
    printf "$formatPrintGreen" "(1) - getbedfiles"
    printf "$formatPrintGreen" "(2) - overlap"
    printf "$formatPrintGreen" "(3) - sortbyname"
    printf "$formatPrintGreen" "(4) - bamtofastq"
    printf "$formatPrintYellow" "(5) - countread"
    printf "$formatPrintYellow" "(6) - mergecountread"
    printf "$formatPrintGreen" "(7) - downsampleauto"
    printf "$formatPrintGreen" "(8) - alignall"
    printf "$formatPrintRed" "(9) - overlapcytoband"
    printf "$formatPrintYellow" "(10) - mergeoverlapcytoband"
    printf "$formatPrintRed" "(11) - overlapfeatures"
    printf "$formatPrintYellow" "(12) - mergeoverlapfeatures"
    printf "$formatPrintYellow" "(13) - getreadnamefastq"
    printf "$formatPrintYellow" "(14) - getreadnamefastqauto"
    printf "$formatPrintGreen" "(15) - overlap_cytoband_count"
    printf "$formatPrintRed" "(16) - extractheaderfastafile"
    printf "$formatPrintBlue" "(17) - replaceheadergenomefasta"
    printf "$formatPrintBlue" "(18) - test"
    printf "$formatPrintRed" "(19) - kallistoindex"
    printf "$formatPrintRed" "(20) - kallistoquant "
    printf "$formatPrintBlue" "(21) - confusionmatrixcytoband"
    printf "$formatPrintBlue" "(22) - mergeconfusionmatrixcytoband"
}

if [ $# -eq 0 ]; then
	__show_options
else
	case "$1" in
	        1)
	                __get_bed_files /data/tmp/clancien/Aligned/CENH3_OVEREXP.1000.csv SRR633614 SRR633612
					__get_bed_files /data/tmp/clancien/Aligned/CENH3_OVEREXP.1000.csv SRR633615 SRR633613
	                ;;
	        2)
					__overlap /data/tmp/clancien/Result/Bed/SRR633612.bed /data/tmp/clancien/Aligned/SRR633612.sorted.bam
					__overlap /data/tmp/clancien/Result/Bed/SRR633613.bed /data/tmp/clancien/Aligned/SRR633613.sorted.bam
					__overlap /data/tmp/clancien/Result/Bed/SRR633614.bed /data/tmp/clancien/Aligned/SRR633614.sorted.bam
					__overlap /data/tmp/clancien/Result/Bed/SRR633615.bed /data/tmp/clancien/Aligned/SRR633615.sorted.bam
					;;
			3)
					__sort_by_name $OVERLAP/SRR633612.overlap.bam
					__sort_by_name $OVERLAP/SRR633613.overlap.bam
					__sort_by_name $OVERLAP/SRR633614.overlap.bam
					__sort_by_name $OVERLAP/SRR633615.overlap.bam
					;;
			4)
					__bam_to_fastq $BAM/SRR633612.overlap.sorted.bam
					__bam_to_fastq $BAM/SRR633613.overlap.sorted.bam
					__bam_to_fastq $BAM/SRR633614.overlap.sorted.bam
					__bam_to_fastq $BAM/SRR633615.overlap.sorted.bam
					;;
			5)
					__count_read $FASTQ/SRR633612.R1.fq
					__count_read $FASTQ/SRR633613.R1.fq
					__count_read $FASTQ/SRR633614.R1.fq
					__count_read $FASTQ/SRR633615.R1.fq
					;;
			6)
					__merge_count_read $FASTQ/countread.txt $FASTQ/SRR633612.countread.txt $FASTQ/SRR633613.countread.txt $FASTQ/SRR633614.countread.txt $FASTQ/SRR633615.countread.txt
					;;
			7)
					__downsample_auto $FASTQ/SRR633614.R1.fq $FASTQ/SRR633614.R2.fq $FASTQ/SRR633612.R1.fq $FASTQ/SRR633612.R2.fq $FASTQ/countread.txt 10 3
					__downsample_auto $FASTQ/SRR633615.R1.fq $FASTQ/SRR633615.R2.fq $FASTQ/SRR633613.R1.fq $FASTQ/SRR633613.R2.fq $FASTQ/countread.txt 10 3
					;;
			8)
					__align_all SRR633614 SRR633612 10 3
					__align_all SRR633615 SRR633613 10 3
					;;
			9)
					__overlap_cytoband SRR633614 SRR633612 10 3
					__overlap_cytoband SRR633615 SRR633613 10 3
					;;
			10)
					__merge_overlap_cytoband SRR633614 SRR633612 10 3 1
					__merge_overlap_cytoband SRR633615 SRR633613 10 3 0
					;;
			11)
					__overlap_features SRR633614 SRR633612 10 3
					__overlap_features SRR633615 SRR633613 10 3
					;;
			12)
					__merge_overlap_features SRR633614 SRR633612 10 3 1
					__merge_overlap_features SRR633615 SRR633613 10 3 0
					;;
			13)
					__get_read_name_id_fastq_file /data/tmp/clancien/Result/Downsampling/SRR633614.SRR633612.R1.1.900000.29913.32169.fq /data/tmp/clancien/Result/Downsampling/SRR633614.SRR633612.R1.1.900000.29913.32169.read.name
					;;
			14)
					__get_read_name_id_fastq_file_auto /data/tmp/clancien/Result/Downsampling
	       			;;
	       	15)
					__overlap_cytoband_count SRR633614 SRR633612 10 3
					__overlap_cytoband_count SRR633615 SRR633613 10 3
					;;
			16)
					__extract_header_fasta_file /data/tmp/clancien/Genomes/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz
					;;
			17)
					__replace_header_genome_fasta /data/tmp/clancien/Genomes/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa /data/tmp/clancien/Genomes/Homo_sapiens.GRCh38.dna_rm.primary_assembly.header.fa
					;;
			18)
					__test /data/tmp/clancien/Genomes/sm/Homo_sapiens.GRCh38.dna_sm.toplevel.fa /data/tmp/clancien/Cytoband/cytoBand.txt /data/tmp/clancien/Genomes/human_genome_split_cytoband.fa
					;;
			19)
					__kallisto_index /data/tmp/clancien/Genomes/human_genome_split_cytoband.fa /data/tmp/clancien/cytoband/cytoband_genome.idx
					;;
			20)
					__kallisto_quant SRR633614 SRR633612 10 3
					#__kallisto_quant SRR633615 SRR633613 10 3
					;;
			21)
					#__confusion_matrix_cytoband SRR633614 SRR633612 10 3
					__confusion_matrix_cytoband SRR633615 SRR633613 10 3
					;;
			22)
					__merge_confusion_matrix_cytoband SRR633614 SRR633612 10 3 1
					#__merge_confusion_matrix_cytoband SRR633615 SRR633613 10 3 0
					;;
	        *)
					__show_options
	                ;;
	esac
fi