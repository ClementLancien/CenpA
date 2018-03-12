#!/usr/bin/env bash
#==> MAKEFILE ??????

#exit as soon as any line in the script fails
set -e 

#exit as soon as any variable is not init
set -o nounset #equals set -u

###############################
#  	     Variables            #
###############################
### SETUP PROJECT PATHS ###
readonly SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
readonly BASE_DIR="$(dirname "$SRC_DIR")"
readonly RAW="${BASE_DIR}/Aligned"
readonly CYTOBAND="${BASE_DIR}/Cytoband"
readonly FEATURES="${BASE_DIR}/Features"
readonly INDEX="${BASE_DIR}/Genomes/GRCh38/dna_sm.primary_assembly/Homo_sapiens.GRCh38"
readonly BED="${BASE_DIR}/Result/Bed"
readonly BAM="${BASE_DIR}/Result/Bam"
readonly FASTQ="${BASE_DIR}/Result/Fastq"
readonly OVERLAP="${BASE_DIR}/Result/Overlap"
readonly DOWNSAMPLING="${BASE_DIR}/Result/Downsampling"
readonly ALIGN="${BASE_DIR}/Result/Align"
readonly STDERR="${BASE_DIR}/stderr"
readonly STDOUT="${BASE_DIR}/stdout"
readonly LOGFILE="${BASE_DIR}/report.log"
readonly GENOME="${BASE_DIR}/Genomes"
readonly KALLISTO_FOLDER="${BASE_DIR}/Result/Kallisto"
readonly INDEX_KALLISTO="${BASE_DIR}/Genomes/cytoband/cytoband_genome.fa"
readonly THOUSANDGENOMES="${BASE_DIR}/Result/1000Genomes"
readonly TMP="${BASE_DIR}/Result/tmp"
### SETUP TOOLS PATHS ### 
readonly SAMTOOLS="/bioinfo/local/build/samtools/samtools-1.3/bin/samtools"
readonly BEDTOOLS="/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin/bedtools"
readonly SEQTK="/data/tmp/clancien/bin/seqtk/seqtk"
readonly BOWTIE2="${BASE_DIR}/bin/bowtie2-2.2.9/bowtie2"
readonly KALLISTO="/data/tmp/clancien/bin/kallisto_linux-v0.43.0/kallisto"

### SETUP FILE PATH ###
readonly cytoband_centromeric="${CYTOBAND}/centromeric.bed"
readonly cytoband_noncentromeric="${CYTOBAND}/noncentromeric.bed"
readonly features_centromeric="${FEATURES}/repeats.bed"
#readonly features_noncentromeric="${FEATURES}/noncentromeric.bed"

### SETUP FORMAT COLOR PRINTF ###
readonly formatPrintYellow="\e[33m\t%s\n\e[m"
readonly formatPrintGreen="\e[32m\t%s\n\e[m"
readonly formatPrintBlue="\e[34m\t%s\n\e[m"
readonly formatPrintRed="\e[31m\t%s\n\e[m"

### Date format dd-mm-yyyy ###
readonly NOW=$(date +%d_%m_%y)

### Time format hh_mm_ss ###
readonly TIME_FORMAT='%H_%M_%S'

### SETUP QSUB OPTIONS ###
readonly overlap_opt="-q batch -l nodes=1:ppn=1,mem=150mb,walltime=01:00:00"
readonly sort_opt="-q batch -l nodes=1:ppn=8,mem=20gb,walltime=02:00:00"
readonly rm_opt="-q batch -l nodes=1,mem=2gb,walltime=01:00:00"
readonly bam_to_fastq_opt="-q batch -l nodes=1:ppn=2,mem=15mb,walltime=00:45:00"
readonly count_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:30:00"
readonly downsample_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:45:00"
readonly aln_opt="-q batch -l nodes=1:ppn=8,mem=16gb,walltime=24:00:00"
readonly index_opt="-q batch -l nodes=1,mem=2gb,walltime=04:00:00"
readonly overlap_cytoband_centromeric_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:20:00"
readonly overlap_cytoband_noncentromeric_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=02:00:00"
readonly overlap_features_centromeric_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:20:00"
#readonly overlap_features_noncentromeric_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=02:00:00"
readonly fastq_read_name_opt="-q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:20:00"
readonly kallisto_index_opt="-q batch -l nodes=1:ppn=1,mem=500gb,walltime=01:00:00"
readonly kallisto_quant_opt="-q batch -l nodes=1:ppn=8,mem=12gb,walltime=02:00:00"

###############################
#  	     CHECK PATHS          #
###############################

### Check stdout and stderr folder ###
if [ ! -d $STDOUT/$NOW/ ]; then mkdir -p $STDOUT/$NOW/; fi
if [ ! -d $STDERR/$NOW/ ]; then mkdir -p $STDERR/$NOW/; fi

###############################
#  	  HELPER DEFINED FUNCTION #
###############################

###############################
#  	  USER DEFINED FUNCTIONS  #
###############################

#
# Purpose: Create bed file.
#
__get_bed_file()
{
	# $1 : CSV File (binned) - name file
	# $2 : CenH3 enrichment - sample name
	# $3 : WT - sample name
	# $4 : output : centromeric bed file - name file
	# $5 : output : non centromeric bed file - name file

	tail -n+2 $1 | awk -v col1="$2" -v col2="$3" -v out1="$4" -v out2="$5" -F, '
	{
		if($1!="Y")
		{
			if($col2==0)
			{
				if($col1 > 2. )
				{
					print $1"\t"$2*1000"\t"$2*1000+1000 > out1 # centromeric file
				}
				else
				{
					print $1"\t"$2*1000"\t"$2*1000+1000 > out2 # noncentromeric file
				}
			}
			else
			{
				ratio=$col1/$col2
				if(ratio>2.)
				{
					print $1"\t"$2*1000"\t"$2*1000+1000 > out1 # centromeric file
				}
				else
				{
					print $1"\t"$2*1000"\t"$2*1000+1000 > out2 # noncentromeric file
				}
			}
		}
		else
		{
			print $1"\t"$2*1000"\t"$2*1000+1000 > out2 # noncentromeric file
		}
	}'

}
#
# Purpose: Return index of elt in an array.
#
__get_index()
{
	local value=$1
	shift
	local my_array=("$@")

	for i in "${!my_array[@]}"
	do
	   	if [[ "${my_array[$i]}" = "$value" ]]
	   	then
	       	echo ${i} ;
	       	break;
	   	fi
	done
}
#
# Purpose: Create bed file.
#
__get_bed_files()
{
	# $1 : CSV File (binned) - name file
	# $2 : CenH3 enrichment - sample name
	# $3 : WT - sample name

	#Check if csv is found and if bed folder if found
	if [ ! -d "$BED" ]; then mkdir -p "$BED"; fi
	if [ ! -f "$1" ]; then echo "$1 not found"; exit; fi

	# Store headers in an array
	declare -a headers
	local j=0
	# (head -n 1 $1 | tr "," "\n")
	for i in $(cat $1 | head -n 1 | tr "," "\n") 
	do
		headers[$j]=$i
		j=$((j+1))
	done

	#Get index in headers from sample -> Awk col index
	#Add +1 because awk index starts 1
	local centromeric=$(( $(__get_index "$2" "${headers[@]}") + 1 ))
	local noncentromeric=$(( $(__get_index "$3" "${headers[@]}") + 1))

	local path_centromeric="${BED}/${2}.bed"
	local path_noncentromeric="${BED}/${3}.bed"

	if [ -f $path_centromeric ]; then rm $path_centromeric; fi
	if [ -f $path_noncentromeric ]; then rm $path_noncentromeric; fi

	__get_bed_file $1 "$centromeric" "$noncentromeric" "$path_centromeric" "$path_noncentromeric"	
}
#
# Purpose: Create bam file from bed and bam files.
#
__overlap ()
{
	# $1 : bed file
	# $2 : bam file

	if [ ! -f "$1" ]; then echo "$1 not found"; exit; fi;
	if [ ! -f "$2" ]; then echo "$2 not found"; exit; fi;

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/overlap_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_${tTIME}_${UUID} $overlap_opt"
	
	#GET file name
	local name=${1##*/} 
	local sample=${name%.bed}

	#GET file path
	local bam_out="${OVERLAP}/${sample}.overlap.bam"

	#GET job
	local ovlp="$SAMTOOLS view -L $1 $2 -b > $bam_out"

	#SUBMIT job
	local ovlp_jid=$(echo "$ovlp" | qsub $opt -N ${sample}.overlap)

	#LOG
	echo "$NOW $tTIME : overlap" >> $LOGFILE 2>&1
	echo "$ovlp | qsub $opt -N ${sample}.overlap " >> $LOGFILE 2>&1

}
#
# Purpose: Sort bam file by name.
#
__sort_by_name ()
{
	if [ -f $1 ]; then echo "$1 not found"; exit; fi

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/sortbyname_${tTIME}_${UUID} -e ${STDERR}/${NOW}/sortbyname_${tTIME}_${UUID}"
	local opt_s="$opt $sort_opt"
	local opt_r="$opt $rm_opt"

	#GET file name
	local name=${1##*/}
	local sample=${name%.overlap.bam}

	#GET file path
	local bam_out="${BAM}/${sample}.overlap.sorted.bam"

	#GET job
	local _sort="$SAMTOOLS sort -m 2G -@ 8 -n $1 > $bam_out"
	local _rm="rm $1"

	#SUBMIT jobs
	local _sort_jid=$(echo "$_sort" | qsub $opt_s -N sort.${sample})
	local _rm_jid=$(echo "$_rm" | qsub $opt_r -N rm.${sample} -W depend=afterok:$_sort_jid)	

	#LOG
	echo "$NOW $tTIME : sort_by_name" >> $LOGFILE 2>&1
	echo "$_sort | qsub $opt_s -N sort.${sample}" >> $LOGFILE 2>&1

	echo "$NOW $tTIME : sort_by_name" >> $LOGFILE 2>&1
	echo "$_rm | qsub $opt_r -N rm.${sample} -W depend=afterok:$_sort_jid" >> $LOGFILE 2>&1
}
#
# Purpose: Create 2 Fastq files from bam file.
#
__bam_to_fastq ()
{
	if [ -f $1 ]; then echo "$1 not found"; exit; fi

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/bam_to_fastq_${tTIME}_${UUID} -e ${STDERR}/${NOW}/bam_to_fastq_${tTIME}_${UUID}"
	local opt_bam_to_fastq="$opt $bam_to_fastq_opt"

	#GET file name
	local name=${1##*/}
	local sample=${name%.overlap.sorted.bam}

	#GET path out
	local fastq1="${FASTQ}/${sample}.R1.fq"
	local fastq2="${FASTQ}/${sample}.R2.fq"

	#GET job
	local bam_to_fastq="$BEDTOOLS bamtofastq -i $1 -fq $fastq1 -fq2 $fastq2"

	#SUBMIT cmd
	local bam_to_fastq_jid="$(echo "$bam_to_fastq" | qsub $opt_bam_to_fastq -N fastq.${sample})"

	#LOG
	echo "$NOW $tTIME : bam_to_fastq" >> $LOGFILE 2>&1
	echo "$bam_to_fastq | qsub $opt_bam_to_fastq -N fastq.${sample}" >> $LOGFILE 2>&1		
}
#
# Purpose : Count number of read in fastq file. Create one file per sample and write the number of read.
#
__count_read ()
{
	if [ -f $1 ]; then echo "$1 not found"; exit; fi

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/count_read${tTIME}_${UUID} -e ${STDERR}/${NOW}/count_read${tTIME}_${UUID}"
	local opt_count="$opt $count_opt"

	#GET file name
	local name=${1##*/}
	local sample=${name%.R*}

	#GET path out
	local file_out="${FASTQ}/${sample}.countread.txt"

	#GET job
	local count="wc -l < $1 > $file_out"

	#SUBMIT job
	local count_jid="$(echo "$count" | qsub $opt_count -N count.${sample})"

	#LOG
	echo "$NOW $tTIME : count_read" >> $LOGFILE 2>&1
	echo "$count | qsub $opt_count -N fastq.${sample}" >> $LOGFILE 2>&1		
}
#
# Purpose : Merge in one file all the sample.countread.txt
#
__merge_count_read ()
{
	local file_out=$1; shift
	if [ -f $file_out ]; then rm $file_out; fi
	while [ $# -gt 0 ]
	do
		#GET file name
		local name=${1##*/}
		local sample=${name%.countread.txt}
		local number=$(head -n 1 $1)
		number=$((number / 4))

		#CMD
		(echo $sample; echo $number) | paste -d',' - - >> $file_out
		rm $1

		shift;
	done
}
#
# Purpose: Function to downsample fastq file.
#
__downsample ()
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/downsample_${tTIME}_${UUID} -e ${STDERR}/${NOW}/downsample_${tTIME}_${UUID}"
	local opt_downsample="$opt $downsample_opt"
	local opt_rm="$opt $rm_opt"

	#GET file name
	local centromeric_name=${1##*/}
	local centromeric_sample=${centromeric_name%.R*}

	local noncentromeric_name=${3##*/}
	local noncentromeric_sample=${noncentromeric_name%.R*}

	#GET paramaters
	#centromeric_sample=$1
	local centromeric_sample_seed_percent=$5
	local centromeric_read_number=$7
	local centromeric_fastq_size_read_number=$(( 4 * $centromeric_read_number ))

	#noncentromeric_sample=$2
	local noncentromeric_sample_seed_percent=$6
	local noncentromeric_read_number=$8
	local noncentromeric_fastq_size_read_number=$(( 4 * $noncentromeric_read_number ))

	local repeat_number=$9

	#GET random seed (INT) for SEQTK
	local centromeric_random_seed=$(( RANDOM ))  #$(( RANDOM % 1000000 ))
	local noncentromeric_random_seed=$(( RANDOM )) #$(( RANDOM % 1000000 ))

	#GET path in
	local centromeric_sample_R1=$1 #"${FASTQ}/${centromeric_sample}.R1.fq"
	local centromeric_sample_R2=$2 #"${FASTQ}/${centromeric_sample}.R2.fq"

	local noncentromeric_sample_R1=$3 #"${FASTQ}/${noncentromeric_sample}.R1.fq"
	local noncentromeric_sample_R2=$4 # "${FASTQ}/${noncentromeric_sample}.R2.fq"

	#GET path out
	local centromeric_sample_R1_tmp_out="${DOWNSAMPLING}/${centromeric_sample}.R1.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.fq"
	local centromeric_sample_R2_tmp_out="${DOWNSAMPLING}/${centromeric_sample}.R2.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.fq"
	
	local noncentromeric_sample_R1_tmp_out="${DOWNSAMPLING}/${noncentromeric_sample}.R1.${repeat_number}.${noncentromeric_read_number}.${noncentromeric_random_seed}.fq"
	local noncentromeric_sample_R2_tmp_out="${DOWNSAMPLING}/${noncentromeric_sample}.R2.${repeat_number}.${noncentromeric_read_number}.${noncentromeric_random_seed}.fq"

	local cat_centromeric_sample_R1_tmp_out_noncentromeric_sample_R1_tmp_out="${DOWNSAMPLING}/${centromeric_sample}.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.${noncentromeric_random_seed}.fq"
	local cat_centromeric_sample_R2_tmp_out_noncentromeric_sample_R2_tmp_out="${DOWNSAMPLING}/${centromeric_sample}.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.${noncentromeric_random_seed}.fq"

	#GET job
	local centromeric_downsample_R1="$SEQTK sample -s $centromeric_random_seed $centromeric_sample_R1 $centromeric_sample_seed_percent | head -n $centromeric_fastq_size_read_number > $centromeric_sample_R1_tmp_out"
	local centromeric_downsample_R2="$SEQTK sample -s $centromeric_random_seed $centromeric_sample_R2 $centromeric_sample_seed_percent | head -n $centromeric_fastq_size_read_number > $centromeric_sample_R2_tmp_out"
	
	local noncentromeric_downsample_R1="$SEQTK sample -s $noncentromeric_random_seed $noncentromeric_sample_R1 $noncentromeric_sample_seed_percent | head -n $noncentromeric_fastq_size_read_number > $noncentromeric_sample_R1_tmp_out"
	local noncentromeric_downsample_R2="$SEQTK sample -s $noncentromeric_random_seed $noncentromeric_sample_R2 $noncentromeric_sample_seed_percent | head -n $noncentromeric_fastq_size_read_number > $noncentromeric_sample_R2_tmp_out"

	local cat_centromeric_downsample_R1_noncentromeric_downsample_R1="cat $noncentromeric_sample_R1_tmp_out $centromeric_sample_R1_tmp_out > $cat_centromeric_sample_R1_tmp_out_noncentromeric_sample_R1_tmp_out"
	local cat_centromeric_downsample_R2_noncentromeric_downsample_R2="cat $noncentromeric_sample_R2_tmp_out $centromeric_sample_R2_tmp_out > $cat_centromeric_sample_R2_tmp_out_noncentromeric_sample_R2_tmp_out"

	local rm_tmp_R1="rm $centromeric_sample_R1_tmp_out $noncentromeric_sample_R1_tmp_out"
	local rm_tmp_R2="rm $centromeric_sample_R2_tmp_out $noncentromeric_sample_R2_tmp_out"

	#SUBMIT job
	local centromeric_downsample_R1_jid="$(echo "$centromeric_downsample_R1" | qsub $opt_downsample -N downsample.${centromeric_sample}.R1.${repeat_number}.${centromeric_read_number})"
	local centromeric_downsample_R2_jid="$(echo "$centromeric_downsample_R2" | qsub $opt_downsample -N downsample.${centromeric_sample}.R2.${repeat_number}.${centromeric_read_number})"
	
	local noncentromeric_downsample_R1_jid="$(echo "$noncentromeric_downsample_R1" | qsub $opt_downsample -N downsample.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number})"
	local noncentromeric_downsample_R2_jid="$(echo "$noncentromeric_downsample_R2" | qsub $opt_downsample -N downsample.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number})"

	local cat_centromeric_downsample_R1_noncentromeric_downsample_R1_jid="$(echo "$cat_centromeric_downsample_R1_noncentromeric_downsample_R1" | qsub $opt_downsample -N cat.${centromeric_sample}.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R1_jid:$noncentromeric_downsample_R1_jid )"
	local cat_centromeric_downsample_R2_noncentromeric_downsample_R2_jid="$(echo "$cat_centromeric_downsample_R2_noncentromeric_downsample_R2" | qsub $opt_downsample -N cat.${centromeric_sample}.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R2_jid:$noncentromeric_downsample_R2_jid)"

	local rm_tmp_R1_jid="$(echo "$rm_tmp_R1" | qsub $opt_rm -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R1.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R1_noncentromeric_downsample_R1_jid)"
	local rm_tmp_R2_jid="$(echo "$rm_tmp_R2" | qsub $opt_rm -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R2.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R2_noncentromeric_downsample_R2_jid)"

	#LOG
	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$centromeric_downsample_R1 | qsub $opt_downsample -N downsample.${centromeric_sample}.R1.${repeat_number}.${centromeric_read_number}" >> $LOGFILE 2>&1
	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$centromeric_downsample_R2 | qsub $opt_downsample -N downsample.${centromeric_sample}.R2.${repeat_number}.${centromeric_read_number}" >> $LOGFILE 2>&1		

	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$noncentromeric_downsample_R1 | qsub $opt_downsample -N downsample.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number}" >> $LOGFILE 2>&1		
	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$noncentromeric_downsample_R2 | qsub $opt_downsample -N downsample.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number}" >> $LOGFILE 2>&1		

	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$cat_centromeric_downsample_R1_noncentromeric_downsample_R1 | qsub $opt_downsample -N cat.${centromeric_sample}.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R1_jid:$noncentromeric_downsample_R1_jid" >> $LOGFILE 2>&1		
	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$cat_centromeric_downsample_R2_noncentromeric_downsample_R2 | qsub $opt_downsample -N cat.${centromeric_sample}.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R2_jid:$noncentromeric_downsample_R2_jid" >> $LOGFILE 2>&1		

	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$rm_tmp_R1 | qsub $opt_rm -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R1.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R1_noncentromeric_downsample_R1_jid" >> $LOGFILE 2>&1		
	echo "$NOW $tTIME : downsample" >> $LOGFILE 2>&1
	echo "$rm_tmp_R2 | qsub $opt_rm -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R2.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R2_noncentromeric_downsample_R2_jid" >> $LOGFILE 2>&1
}
#
# Purpose : Return number of total read from sample ($1) in countread file ($2).
#
__get_read_number()
{
	awk -v sample_name="$1" -F, '
		{
			if($1==sample_name)
				{
					print $2
					exit
				}
		}
	' $2
}
#
# Purpose : Return seed percent for downsample (number of read / number of total read).
#
__get_seed_percent ()
{
	echo "$(echo $1/$2 | bc -l)"
}
#
# Purpose : Function to downsample fastq file with different size and number of repeats.
#
__downsample_auto ()
{
	if [ ! -f $1 ]; then echo "$1 not found"; exit; fi
	if [ ! -f $2 ]; then echo "$2 not found"; exit; fi
	if [ ! -f $3 ]; then echo "$3 not found"; exit; fi
	if [ ! -f $4 ]; then echo "$4 not found"; exit; fi

	#GET file name
	local centromeric_name=${1##*/}
	local centromeric_sample=${centromeric_name%.R*}

	local noncentromeric_name=${3##*/}
	local noncentromeric_sample=${noncentromeric_name%.R*}

	#GET total read number
	local centromeric_read_number=$(__get_read_number $centromeric_sample $5)
	local noncentromeric_read_number=$(__get_read_number $noncentromeric_sample $5)

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$6*900000;sample_centromeric_size+=900000))
	do
		#Get seed numbers
		local sample_noncentromeric_size=$(( 90000000 - $sample_centromeric_size ))
		local centromeric_seed_percent=$(__get_seed_percent $sample_centromeric_size $centromeric_read_number)
		local noncentromeric_seed_percent=$(__get_seed_percent $sample_noncentromeric_size $noncentromeric_read_number)

		#Loop for 
		for repeat in $(seq 1 $7)
		do				
			__downsample $1 $2 $3 $4 $centromeric_seed_percent $noncentromeric_seed_percent $sample_centromeric_size $sample_noncentromeric_size $repeat
		done
	done

}
#
# Purpose : 
#
__align ()
{
	if [ ! -f "$1" ]; then echo "$1 not found"; exit; fi
	if [ ! -f "$2" ]; then echo "$2 not found"; exit; fi

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/align_${tTIME}_${UUID} -e ${STDERR}/${NOW}/align_${tTIME}_${UUID}"
	local opt_aln="$opt $aln_opt"
	local opt_sort="$opt $sort_opt"
	local opt_index="$opt $index_opt"
	local opt_rm="$opt $rm_opt"

	#GET Parameters
	R1=$1
	R2=$2

	centromeric_sample=$3
	noncentromeric_sample=$4

	repeat=$5
	sample_centromeric_size=$6

	#GET path out
	aln_out=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${sample_centromeric_size}.bam
	sort_out=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${sample_centromeric_size}.sorted.bam

	#GET cmd
	_aln="$BOWTIE2 --threads 8 --very-sensitive -x $INDEX -1 ${R1} -2 ${R2} | $SAMTOOLS view -b - > $aln_out"
	
	_sort="$SAMTOOLS sort -m 2G -@ 8 $aln_out $sort_out"
	
	_index="$SAMTOOLS index -b $sort_out"
	
	_rm="rm $aln_out"

	#SUBMIT jobs
	aln_jid="$(echo $_aln | qsub $opt_aln -N aln.${1}.${2}.${repeat}.${sample_centromeric_size})"
	
	sort_jid="$(echo $_sort | qsub $opt_sort -N sort.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$aln_jid)"
	
	index_jid="$(echo $_index | qsub $opt_index -N index.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"
	
	rm_jid="$(echo $_rm | qsub $opt_rm -N rm.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"

	#LOG
	echo "$NOW $tTIME : align" >> $LOGFILE 2>&1
	echo "$_aln | qsub $opt_aln -N aln.${1}.${2}.${repeat}.${sample_centromeric_size}" >> $LOGFILE 2>&1		
	
	echo "$NOW $tTIME : align" >> $LOGFILE 2>&1
	echo "$_sort | qsub $opt_sort -N sort.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$aln_jid " >> $LOGFILE 2>&1		

	echo "$NOW $tTIME : align" >> $LOGFILE 2>&1
	echo "$_index | qsub $opt_index -N index.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid" >> $LOGFILE 2>&1		
	
	echo "$NOW $tTIME : align" >> $LOGFILE 2>&1
	echo "$_rm | qsub $opt_rm -N rm.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)" >> $LOGFILE 2>&1		

}
#
# Purpose :
#
__align_all ()
{
	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			local seed1="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 6)"
			local seed2="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 7)"
			
			local R1_path=$DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq
			local R2_path=$DOWNSAMPLING/${1}.${2}.R2.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq
			
			__align $R1_path $R2_path $1 $2 $repeat $sample_centromeric_size
		done
	done
}
#
# Purpose : Count number of alignment to compute ratio.
#
__number_of_alignement()
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/number_of_alignement_${tTIME}_${UUID} -e ${STDERR}/${NOW}/number_of_alignement_${tTIME}_${UUID}"
	local opt_count="$opt $count_opt" 

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			local bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"
			local file_out="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.count.txt"
			
			#GET job
			local count="$SAMTOOLS view $bam_in | wc -l > $file_out"

			#SUBMIT job
			local count_jid="$(echo "$count" | qsub $opt_count -N "$file_out")"

			#LOG
			echo "$NOW $tTIME : number_of_alignement" >> $LOGFILE 2>&1
			echo "$count | qsub $opt_count -N $file_out" >> $LOGFILE 2>&1		
		done
	done
}
#
# Purpose : Merge all count file into one.
#
__merge_number_of_alignement()
{	
	local file_out="$ALIGN/align.count.txt"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			local file_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.count.txt"
			local number=$(head -n 1 $file_in)

			$(echo "${1}.${2}.subsample.${repeat}.${sample_centromeric_size},${number}" >> "$file_out")
			rm $file_in
		done
	done
}
#
# Purpose : Count Observed ratio (centromere/noncentromere) with cytoband
#
__overlap_cytoband() #centromeric #noncentromeric ##percent #repeat 
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/overlap_cytoband_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_cytoband_${tTIME}_${UUID}"
	local opt_overlap_cytoband_centromeric="$opt $overlap_cytoband_centromeric_opt" 
 	local opt_overlap_cytoband_noncentromeric="$opt $overlap_cytoband_noncentromeric_opt"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			local bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			local ovlp_centromeric="$SAMTOOLS view -L $cytoband_centromeric $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.cytoband"
			local ovlp_noncentromeric="$SAMTOOLS view -L $cytoband_noncentromeric $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.cytoband"

			#SUBMIT job
			local ovlp_centromeric_jid="$(echo "$ovlp_centromeric" | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} )"
			local ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"
		
			#LOG
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_centromeric | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_noncentromeric | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
		done
	done
}
#
# Purpose :
#
__merge_overlap_cytoband()
{
	centromeric_sample=$1
	noncentromeric_sample=$2
	total=$(($3*900000))
	file_out=$ALIGN/cytoband.count.txt

	if [ $5 -eq 1 ] ; then
		if [ -f $file_out ]; then rm $file_out; fi;
		echo "centromeric_sample,noncentromeric_sample,read_number,subsample,centromeric_number,noncentromeric_number,repeat,observed_ratio,theorical_ratio" >> $file_out
	fi

	for((i=900000;i<=$total;i+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			centromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.centromeric.cytoband
			noncentromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.noncentromeric.cytoband

			centromeric_number=$(head -n 1 $centromeric_path)
			noncentromeric_number=$(head -n 1 $noncentromeric_path)

			observed_ratio=$(echo $centromeric_number/$noncentromeric_number | bc -l)
			theorical_ratio=$(echo $i/90000000 | bc -l)

			echo ${centromeric_sample},${noncentromeric_sample},${i},${repeat},${centromeric_number},${noncentromeric_number},${repeat},${observed_ratio},${theorical_ratio} >> $file_out
		done
	done
}
# __overlap_cytoband() #centromeric #noncentromeric ##percent #repeat 
# {
# 	#TIMESTAMPING
# 	local tTIME=$(date +"${TIME_FORMAT}")
# 	local UUID=$(cat /proc/sys/kernel/random/uuid)

# 	#INIT qsub option
# 	local opt="-o ${STDOUT}/${NOW}/overlap_cytoband_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_cytoband_${tTIME}_${UUID}"
# 	local opt_overlap_cytoband_centromeric="$opt $overlap_cytoband_centromeric_opt" 
#  	local opt_overlap_cytoband_noncentromeric="$opt $overlap_cytoband_noncentromeric_opt"

# 	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
# 	do
# 		for repeat in $(seq 1 $4)
# 		do
# 			#GET PATH IN
# 			local bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

# 			#GET cmd
# 			local ovlp_centromeric="$SAMTOOLS view -L $cytoband_centromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.cytoband"
# 			local ovlp_noncentromeric="$SAMTOOLS view -L $cytoband_noncentromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.cytoband"

# 			#SUBMIT job
# 			local ovlp_centromeric_jid="$(echo "$ovlp_centromeric" | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} )"
# 			local ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"
		
# 			#LOG
# 			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
# 			echo "$ovlp_centromeric | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
# 			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
# 			echo "$ovlp_noncentromeric | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
# 		done
# 	done
# }
# __merge_overlap_cytoband()
# {
# 	centromeric_sample=$1
# 	noncentromeric_sample=$2
# 	total=$(($3*900000))
# 	file_out=$ALIGN/cytoband.count.txt

# 	if [ $4 -eq 1 ] ; then
# 		if [ -f $file_out ]; then rm $file_out; fi;
# 		echo "centromeric_sample,noncentromeric_sample,read_number,subsample,centromeric_number,noncentromeric_number,repeat,observed_ratio,theorical_ratio" >> $file_out
# 	fi

# 	for((i=900000;i<=$total;i+=900000))
# 	do
# 		for repeat in $(seq 1 $4)
# 		do
# 			centromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.centromeric.cytoband
# 			noncentromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.noncentromeric.cytoband

# 			centromeric_number=$(head -n 1 $centromeric_path)
# 			noncentromeric_number=$(head -n 1 $noncentromeric_path)

# 			observed_ratio=$(echo $centromeric_number/$noncentromeric_number | bc -l)
# 			theorical_ratio=$(echo $i/90000000 | bc -l)

# 			echo ${centromeric_sample},${noncentromeric_sample},${i},${repeat},${centromeric_number},${noncentromeric_number},${repeat},${observed_ratio},${theorical_ratio} >> $file_out
# 		done
# 	done
# }
__confusion_matrix_cytoband()
{
	# $1 centromere
	#
	# Use of : comm -12 < (sort file1 | uniq) < (sort file2 | uniq)
	#

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/confusion_matrix_cytoband_${tTIME}_${UUID} -e ${STDERR}/${NOW}/confusion_matrix_cytoband_${tTIME}_${UUID}"
	local opt_count="$opt $count_opt" 

	centromeric_sample=$1
	noncentromeric_sample=$2

	for((i=900000;i<=$3*900000;i+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			local centromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.centromeric.cytoband.bam
			local noncentromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.noncentromeric.cytoband.bam
			#echo "^$1"
			c_c="$SAMTOOLS view $centromeric_path | grep "^$1" | wc -l > $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_c"
			c_nc="$SAMTOOLS view $centromeric_path | grep "^$2" | wc -l >  $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_nc"

			#echo "${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c $c_c $c_nc" >> $output

			nc_c="$SAMTOOLS view $noncentromeric_path | grep "^$1" | wc -l > /data/tmp/clancien/Result/tmp/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_c"
			nc_nc="$SAMTOOLS view $noncentromeric_path | grep "^$2" | wc -l > /data/tmp/clancien/Result/tmp/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_nc"

			c_c_job_id=$(echo $c_c | qsub $opt_count -N ${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_c)
			c_nc_job_id=$(echo $c_nc | qsub $opt_count -N ${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_nc)
			nc_c_job_id=$(echo $nc_c | qsub $opt_count -N ${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_c)
			nc_nc_job_id=$(echo $nc_nc | qsub $opt_count -N ${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_nc)
			#echo "${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc $nc_c $nc_nc" >> $output
		done
	done
}
__merge_confusion_matrix_cytoband()
{
	# $1 centromere
	#
	# Use of : comm -12 < (sort file1 | uniq) < (sort file2 | uniq)
	#

	centromeric_sample=$1
	noncentromeric_sample=$2

	output=$ALIGN/confusion_matrix_cytoband
	if [ $5 -eq 1 ] ; then
		if [ -f $output ]; then rm $output; fi;
	fi
	centromeric_sample=$1
	noncentromeric_sample=$2

	for((i=900000;i<=$3*900000;i+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			local c_c=$( head -n 1 $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_c)
			local c_nc=$( head -n 1 $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c_nc)
			local nc_c=$( head -n 1 $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_c)
			local nc_nc=$( head -n 1 $TMP/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc_nc)
			
			echo "${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.c,$c_c,$c_nc" >> $output
			echo "${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.nc,$nc_c,$nc_nc" >> $output
		done
	done
}
#
# Purpose
#
__overlap_features() #centromeric #noncentromeric ##percent #repeat 
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/overlap_features_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_features_${tTIME}_${UUID}"
	local opt_overlap_features_centromeric="$opt $overlap_features_centromeric_opt"
	#local opt_overlap_features_noncentromeric="$opt $overlap_features_noncentromeric_opt"
	local opt_overlap_features_noncentromeric="$opt $count_opt"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			#ovlp_centromeric="$SAMTOOLS view -L $features_centromeric $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.features"
			#ovlp_noncentromeric="$SAMTOOLS view -L $features_noncentromeric $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.features"
			##ovlp_noncentromeric="$SAMTOOLS view $bam_in | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.features "

			ovlp_centromeric="$SAMTOOLS view -L $features_centromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.features.bam"
			ovlp_noncentromeric="$SAMTOOLS view $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.features.bam"

			#SUBMIT job
			ovlp_centromeric_jid="$(echo "$ovlp_centromeric" | qsub $opt_overlap_features_centromeric -N features.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} )"
			#ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $features_noncentromeric -N features.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"
			ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $opt_overlap_features_noncentromeric -N features.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"

			#LOG
			echo "$NOW $tTIME : overlap_features" >> $LOGFILE 2>&1
			echo "$ovlp_centromeric | qsub $opt_overlap_features_centromeric -N features.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size}" >> $LOGFILE 2>&1
			echo "$NOW $tTIME : overlap_features" >> $LOGFILE 2>&1
			echo "$ovlp_noncentromeric | qsub $opt_overlap_features_noncentromeric -N features.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
		done
	done
}
#
# Purpose :
#
__merge_overlap_features()
{
	local centromeric_sample=$1
	local noncentromeric_sample=$2
	local total=$(($3*900000))
	local file_out=$ALIGN/features.count.txt

	if [ $5 -eq 1 ] ; then
		if [ -f $file_out ]; then rm $file_out; fi;
		echo "centromeric_sample,noncentromeric_sample,read_number,subsample,centromeric_number,noncentromeric_number,repeat,observed_ratio,theorical_ratio" >> $file_out
	fi

	for((i=900000;i<=$total;i+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			local centromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.centromeric.features
			local noncentromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.noncentromeric.features

			local centromeric_number=$(head -n 1 $centromeric_path)
			local noncentromeric_number=$(head -n 1 $noncentromeric_path)
			noncentromeric_number=$(( noncentromeric_number - centromeric_number ))

			local observed_ratio=$(echo $centromeric_number/$noncentromeric_number | bc -l)
			local theorical_ratio=$(echo $i/90000000 | bc -l)

			echo ${centromeric_sample},${noncentromeric_sample},${i},${repeat},${centromeric_number},${noncentromeric_number},${repeat},${observed_ratio},${theorical_ratio} >> $file_out
		done
	done
}
#
# Purpose : Get Read name from Fastq file
#
__get_read_name_id_fastq_file()
{
	# $1 input
	# $2 output
	awk -v output=$2 -F'[@/]' 'NR%4 == 1 { print $2 > output }' $1
}
#
# Purpose : Get Auto getreadnameidfastqfile
#
__get_read_name_id_fastq_file_auto()
{
	for fastq in ${1}/*.fq
	do
		#TIMESTAMPING
		local tTIME=$(date +"${TIME_FORMAT}")
		local UUID=$(cat /proc/sys/kernel/random/uuid)

		#INIT qsub option
		local opt="-o ${STDOUT}/${NOW}/get_read_name_id_fastq_file_auto_${tTIME}_${UUID} -e ${STDERR}/${NOW}/get_read_name_id_fastq_file_auto_${tTIME}_${UUID}"
		local opt_fastq_read_name="$opt $fastq_read_name_opt" 
		
		#GET job

		local job="awk -v output=${fastq%.fq}.read.name -F'[@/]' 'NR%4 == 1 { print \$2 > output }' $fastq "

		#SUBMIT job
		local job_id="$(echo "$job" | qsub $opt_fastq_read_name -N ${fastq##*/})"

		#LOG
		echo "$NOW $tTIME : get_read_name_id_fastq_file_auto" >> $LOGFILE 2>&1
		echo "awk -v output=${fastq%.fq}.read.name -F'[@/]' 'NR%4 == 1 { print \$2 > ${fastq%.fq}.read.name }' $fastq | qsub $opt_fastq_read_name -N ${fastq##*/} " >> $LOGFILE 2>&1
	done
}
#
# Purpose :
#
__get_read_overlap_cytoband()
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/overlap_cytoband_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_cytoband_${tTIME}_${UUID}"
	local opt_overlap_cytoband_centromeric="$opt $overlap_cytoband_centromeric_opt" 
 	local opt_overlap_cytoband_noncentromeric="$opt $overlap_cytoband_noncentromeric_opt"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			local bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			local ovlp_centromeric="$SAMTOOLS view -L $cytoband_centromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.cytoband"
			local ovlp_noncentromeric="$SAMTOOLS view -L $cytoband_noncentromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.cytoband"

			#SUBMIT job
			local ovlp_centromeric_jid="$(echo "$ovlp_centromeric" | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} )"
			local ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"

			#LOG
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_centromeric | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_noncentromeric | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
		done
	done
}
#
# Purpose : 
#
__overlap_cytoband_count() #centromeric #noncentromeric ##percent #repeat 
{
	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/overlap_cytoband_count_${tTIME}_${UUID} -e ${STDERR}/${NOW}/overlap_cytoband_count_${tTIME}_${UUID}"
	local opt_overlap_cytoband_centromeric="$opt $overlap_cytoband_centromeric_opt" 
 	local opt_overlap_cytoband_noncentromeric="$opt $overlap_cytoband_noncentromeric_opt"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			local bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			local ovlp_centromeric="$SAMTOOLS view -L $cytoband_centromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric.cytoband.bam"
			local ovlp_noncentromeric="$SAMTOOLS view -L $cytoband_noncentromeric $bam_in -b > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric.cytoband.bam"

			#SUBMIT job
			local ovlp_centromeric_jid="$(echo "$ovlp_centromeric" | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} )"
			local ovlp_noncentromeric_jid="$(echo "$ovlp_noncentromeric" | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"

			#LOG
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_centromeric | qsub $opt_overlap_cytoband_centromeric -N cytoband.centromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
			echo "$NOW $tTIME : overlap_cytoband" >> $LOGFILE 2>&1
			echo "$ovlp_noncentromeric | qsub $opt_overlap_cytoband_noncentromeric -N cytoband.noncentromeric.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} " >> $LOGFILE 2>&1
		done
	done
}
#
# Purpose: Get genome order by chr and band
#
__extract_header_fasta_file()
{
	#Check if compress file
	if [ ${1:${#1}<3?0:-3} == ".gz" ]; then
		grep "^>" <(zcat $1) > $2
	else
		grep "^>" $1 > $2
	fi
}
#
# Purpose : Create new fasta file by replacing what is after blank by nothing(header).
#
__replace_header_genome_fasta()
{
	sed 's/\s.*$//' $1 > $2
}
#
# Purpose : Get fasta reference sort by chromosome and cytogenetic band.
#
__test()
{

	if [ -f $3 ]; then rm $3; fi

	for header in $(ls | grep "^>" $1)
	do
		#echo $header
		
		if [ ${header:0:1} == ">" ]; then
			local chromosome=${header:1}
			#echo ${chromosome:0:3}
		
			if [ ${chromosome:0:3} != "CHR" ]; then

				local chromosome=$(sed -r 's/\./v/' <<< $chromosome)

				if [ $chromosome == "MT" ]; then 
					chromosome=${chromosome:0:1}; 
				fi

				local i=0
				while read line
				do
					local chr=$(echo $line | cut -d ' ' -f 1)

					if [ "chr${chromosome}" == "$chr" ] || [ ${#chromosome} -gt 2 ] && [[  "$chr" =~ $chromosome ]]; then
						local array=($line)

						str=$(sed -r 's/ /\_/g' <<< $line)

						$($SAMTOOLS faidx $1 "${header:1}:${array[1]}-${array[2]}" > $GENOME/${str})
						#echo $str
						#sed -i "1s/.*/>$str/" $GENOME/${str}
						sed -i "1s/.*/>$str/" $GENOME/${str}
						cat $GENOME/${str} >> $GENOME/${array[0]}.tmp
						rm $GENOME/${str}

					fi
				done < $2

				cat $GENOME/*.tmp >> $3
				rm $GENOME/*.tmp
			fi
		# printf "$formatPrintRed" "END"
		#exit
		fi
	done
}
#
# Purpose : Create kallisto index
#
__kallisto_index()
{
	#$1 fasta name
	#$2 idx file

	#TIMESTAMPING
	local tTIME=$(date +"${TIME_FORMAT}")
	local UUID=$(cat /proc/sys/kernel/random/uuid)

	#INIT qsub option
	local opt="-o ${STDOUT}/${NOW}/kalisto_index_${tTIME}_${UUID} -e ${STDERR}/${NOW}/kalisto_index_${tTIME}_${UUID}"
	local opt_kallisto_index="$opt $kallisto_index_opt" 

	#GET job
	job="$KALLISTO index -i $2 $1"

	#SUBMIT job
	job_id=$(echo "$job" | qsub $opt_kallisto_index -N kallisto_index )

	#LOG
	echo "$NOW $tTIME : kalisto_index" >> $LOGFILE 2>&1
	echo "$job | qsub $opt_kallisto_index -N kallisto_index " >> $LOGFILE 2>&1
}
#
# Purpose : Runs the quantification kallisto algorithm for each downsample
#
__kallisto_quant()
{
	#$1 centromeric
	#$2 noncentromeric
	#$3 max size*900000
	#$4 number of repeat

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#TIMESTAMPING
			local tTIME=$(date +"${TIME_FORMAT}")
			local UUID=$(cat /proc/sys/kernel/random/uuid)

			#INIT qsub option
			local opt="-o ${STDOUT}/${NOW}/kallisto_quant_${tTIME}_${UUID} -e ${STDERR}/${NOW}/kallisto_quant_${tTIME}_${UUID}"
			local opt_kallisto_quant="$opt $kallisto_quant_opt"

			#GET Path in
			local seed1="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 6)"
			local seed2="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 7)"
			
			local R1_path=$DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq
			local R2_path=$DOWNSAMPLING/${1}.${2}.R2.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq

			#GET Path out
			local prefix=$KALLISTO_FOLDER/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}

			#GET job
			local quant_rf="$kallisto quant --rf-stranded --bias -t 8 -i $INDEX_KALLISTO $R1_path $R2_path -o $prefix.rf"
    		local quant_fr="$kallisto quant --fr-stranded --bias -t 8 -i $INDEX_KALLISTO $R1_path $R2_path -o $prefix.fr"
   			local quant_un="$kallisto quant --bias -t 8 -i $INDEX_KALLISTO $R1_path $R2_path -o $prefix.un"

   			#SUBMIT job
   			local quant_rf_jid="$(echo $quant_rf | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.rf )"
			local quant_fr_jid="$(echo $quant_fr | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.fr )"
			local quant_un_jid="$(echo $quant_un | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.un )"
		
			#LOG
			echo "$NOW $tTIME : kallisto_quant" >> $LOGFILE 2>&1
			echo "$quant_rf | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.rf " >> $LOGFILE 2>&1
			echo "$NOW $tTIME : kallisto_quant" >> $LOGFILE 2>&1
			echo "$quant_fr | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.fr " >> $LOGFILE 2>&1		
			echo "$NOW $tTIME : kallisto_quant" >> $LOGFILE 2>&1
			echo "$quant_un | qsub $opt_kallisto_quant -N $${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.un " >> $LOGFILE 2>&1
			exit
		done
	done
}
#
# Purpose
#
_1000_genome()
{
	echo "t"
}

###################################
# Main Script Logic Starts Here   #
###################################

# _menu_get_bed_files()
# {
# 	printf "$formatPrintBlue" "some information"
# 	read -p 'CSV File: ' csvfile
# 	read -p 'CenH3 sample: ' centromeric_sample
# 	read -p 'WT sample ' noncentromeric_sample
# 	__get_bed_files $csvfile $centromeric_sample $noncentromeric_sample
# }

# __main_menu() 
# {
# 	printf "$formatPrintYellow" "Welcome to the menu of blabla"
	
# 	local option=0
# 	until [ "$option" = "4" ]; do
# 		printf "%s\n" "  1.) Get Bed Files"
# 		printf "%s\n" "  2.) Remove user"
# 		printf "%s\n" "  3.) Update user"
# 		printf "%s\n" "  4.) Quit"

# 		read -p 'Select option number: ' option

# 		case $option in
# 	    	1 ) 
# 				_menu_get_bed_files
# 				#printf "%s\n" "Choice 1"
# 				;;
# 	    	2 ) 
# 				printf "%s\n" "Choice 2" 
# 				;;
# 	    	3 ) 
# 				printf "%s\n" "Choice 3"
# 				;;
# 	    	4 ) 
# 				printf "%s\n" "Exit"; 
# 				exit
# 				;;
# 	    	* ) 
# 				printf "$formatPrintRed" "Invalid Option"
# 				printf "%s\n" "Please enter 1, 2, 3, or 4"
# 				;;
# 		esac
# 	done
# }
# __main_menu

# case "$1" in
#         getbedfiles)
#                 __get_bed_files $2 $3 $4
#                 ;;
#         overlap)
# 				__overlap
# 				;;
# 		sortbyname)
# 				__sort_by_name
# 				;;
#         *)
#                 echo "Usage: $0 {getbedfiles|overlap|sortbyname}"
# esac