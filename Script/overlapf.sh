#!/usr/bin/env bash

# exit when any command fails
set -e
set -o nounset


readonly CWD=$(pwd)
readonly BASE_DIR=$(dirname "$CWD")

readonly RAW="${BASE_DIR}/Aligned"
readonly MAPPEDREAD="${BASE_DIR}/MappedRead"
readonly BED="${BASE_DIR}/Result/Bed"
readonly BAM="${BASE_DIR}/Result/Bam"
readonly FASTQ="${BASE_DIR}/Result/Fastq"
readonly OVERLAP="${BASE_DIR}/Result/Overlap"

readonly DOWNSAMPLING="${BASE_DIR}/Result/Downsampling"
readonly INDEX="${BASE_DIR}/Genomes/GRCh38/dna_sm.primary_assembly/Homo_sapiens.GRCh38"
readonly ALIGN="${BASE_DIR}/Result/Align"

readonly SAMTOOLS="/bioinfo/local/build/samtools/samtools-1.3/bin/samtools"
#readonly SAMTOOLS="/bioinfo/local/build/Centos/samtools/samtools-1.3/bin/samtools"
readonly BEDTOOLS="/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin/bedtools"
readonly SEQTK="/data/tmp/clancien/Tools/seqtk/seqtk"
readonly BOWTIE2="$BASE_DIR/bin/bowtie2-2.2.9/bowtie2"


readonly day=$(date +%d_%m_%y)

function _overlapf() 
{
	#Check if qsub stdout and stderr folders exist else create these folders
	if [ ! -d $OVERLAP/stdout/ ]; then mkdir -p $OVERLAP/stdout/; fi
	if [ ! -d $OVERLAP/stderr/ ]; then mkdir -p $OVERLAP/stderr/; fi

	#INIT qsub options
	opt="-o $OVERLAP/stdout/$day -e $OVERLAP/stdout/$day"
	ovlp_opt="$opt -q batch -l nodes=1:ppn=1,mem=150mb,walltime=01:00:00"

}
