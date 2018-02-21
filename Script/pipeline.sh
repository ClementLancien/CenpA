#!/usr/bin/env bash

# exit when any command fails
set -e
set -o nounset

#####################################################################################################################################################
# 				____/\\\\\\\\\______/\\\________/\\\__/\\\\\_____/\\\_ 			__/\\\\____________/\\\\__/\\\\\\\\\\\\\\\_ 				        #        
# 				 __/\\\///////\\\___\/\\\_______\/\\\_\/\\\\\\___\/\\\_    		 _\/\\\\\\________/\\\\\\_\/\\\///////////__ 				        #
# 				  _\/\\\_____\/\\\___\/\\\_______\/\\\_\/\\\/\\\__\/\\\_      	  _\/\\\//\\\____/\\\//\\\_\/\\\_____________ 				        #
# 				   _\/\\\\\\\\\\\/____\/\\\_______\/\\\_\/\\\//\\\_\/\\\_     	   _\/\\\\///\\\/\\\/_\/\\\_\/\\\\\\\\\\\_____     					#
# 				    _\/\\\//////\\\____\/\\\_______\/\\\_\/\\\\//\\\\/\\\_    		_\/\\\__\///\\\/___\/\\\_\/\\\///////______    					#
# 				     _\/\\\____\//\\\___\/\\\_______\/\\\_\/\\\_\//\\\/\\\_   		 _\/\\\____\///_____\/\\\_\/\\\_____________ 				    #
# 				      _\/\\\_____\//\\\__\//\\\______/\\\__\/\\\__\//\\\\\\_  		  _\/\\\_____________\/\\\_\/\\\_____________ 				    #
# 				       _\/\\\______\//\\\__\///\\\\\\\\\/___\/\\\___\//\\\\\_ 		   _\/\\\_____________\/\\\_\/\\\\\\\\\\\\\\\_ 				    #
# 				        _\///________\///_____\/////////_____\///_____\/////__  	    _\///______________\///__\///////////////__ 				#
#####################################################################################################################################################

readonly CWD=$(pwd)
readonly BASE_DIR=$(dirname "$CWD")


readonly RAW="${BASE_DIR}/Aligned"
readonly MAPPEDREAD="${BASE_DIR}/MappedRead"
readonly BED="${BASE_DIR}/Result/Bed"
readonly BAM="${BASE_DIR}/Result/Bam"
readonly FASTQ="${BASE_DIR}/Result/Fastq"
readonly OVERLAP="${BASE_DIR}/Result/Overlap"
readonly CYTOBAND="${BASE_DIR}/Cytoband"
# readonly FEATURES

readonly DOWNSAMPLING="${BASE_DIR}/Result/Downsampling"
readonly INDEX="${BASE_DIR}/Genomes/GRCh38/dna_sm.primary_assembly/Homo_sapiens.GRCh38"
readonly ALIGN="${BASE_DIR}/Result/Align"

readonly SAMTOOLS="/bioinfo/local/build/samtools/samtools-1.3/bin/samtools"
#readonly SAMTOOLS="/bioinfo/local/build/Centos/samtools/samtools-1.3/bin/samtools"
readonly BEDTOOLS="/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin/bedtools"
readonly SEQTK="/data/tmp/clancien/Tools/seqtk/seqtk"
readonly BOWTIE2="$BASE_DIR/bin/bowtie2-2.2.9/bowtie2"
readonly day=$(date +%d_%m_%y)

readonly merge_count_out="$BAM/countread.txt"


# opt="-o $BAM/stdout/$day -e $BAM/stdout/$day"
# sub_opt="$opt -q batch -l nodes=1:ppn=1,mem=150mb,walltime=01:00:00"



# readonly SEQTK="/data/tmp/clancien/Tools/seqtk/seqtk"

# random_seed=echo $(( RANDOM % 100000 ))


# R1="$BAM/SRR633614.R1.fq"
# R2="$BAM/SRR633614.R2.fq"

# job_id1="$(echo "$SEQTK sample -s $random_seed $R1 0.015 | head -n 4000000 > $BAM/SRR633614.R1.random.sub.fq" | qsub $sub_opt -N sub.R1)"
# job_id2="$(echo "$SEQTK sample -s $random_seed $R2 0.015 | head -n 4000000 > $BAM/SRR633614.R2.random.sub.fq" | qsub $sub_opt -N sub.R2)"

# exit

#seqtk sample -s100 read1.fq 0.02 > sub1.fq
#seqtk sample -s100 read2.fq 0.02 > sub2.fq



#_sort_jid=$(echo "$_sort" | qsub $sort_opt -N sort.${sample})




function __get_bed_file()
{
	# $1 CSV File
	# $2 output : centromeric bed file
	# $3 output : non centromeric bed file

	tail -n+2 $1 | awk -v col1="$2" -v col2="$3" -v out1="$4" -v out2="$5"  -F, '

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

function __get_index()
{
	value=$1
	shift
	my_array=("$@")
	for i in "${!my_array[@]}"
	do
	   	if [[ "${my_array[$i]}" = "$value" ]]
	   	then
	       	echo ${i} ;
	   	fi
	done

}

function __get_bed_files()
{

	#Check if csv is found and if bed folder if found
	#bed_folder_path="${BASE_DIR/Result/Bed"

	if [ ! -f $1 ]; then echo "CSV File not found";	fi
	if [ ! -d "$BED" ]; then mkdir -p "$BED"; fi

	# Arrays headers contains first line of csv file
	declare -a headers
	local j=0

	for i in $(cat $1 | head -n 1 | tr "," "\n")
	do
		headers[$j]=$i
		let j=j+1
	done
	#echo "${!headers[@]}"
	#echo "${headers[@]}"


	#Get index in headers from sample -> Awk col index
	#Add +1 because awk index starts 1
	centromeric=$(( $(__get_index "$2" "${headers[@]}") + 1 ))
	noncentromeric=$(( $(__get_index "$3" "${headers[@]}") + 1))


	path_centromeric="$BED/$2.bed"
	path_noncentromeric="$BED/$3.bed"

	__get_bed_file $1 "$centromeric" "$noncentromeric" "$path_centromeric" "$path_noncentromeric"
}

# bed_raw=$BASE_DIR/Aligned/CENH3_OVEREXP.1000.csv
# echo "here"
# __get_bed_files $bed_raw SRR633614 SRR633612
# echo "here"
# __get_bed_files $bed_raw SRR633615 SRR633613
# exit

#bed_raw=$BASE_DIR/Aligned/CENH3_OVEREXP.1000.csv

#chromosome,start,SRR633612,SRR633613,SRR633614,SRR633615,SRR633616,SRR633617,SRR633618,SRR633619,SRR891239,SRR891240
#0		   ,1    ,2        ,3        ,4        ,5        ,6 
#NA        ,NA   ,NC       ,NC       ,C        ,C        ,NA
#bed_raw centromeric_sample noncentromeric_sample




# function __menu()
# {
# 	title="Make text"
# 	prompt="Pick an option:"
# 	options=("Get Bed Files" "B" "C")

# 	echo "$title"
# 	PS3="$prompt "
# 	select opt in "${options[@]}" "Quit"; do 

# 	    case "$REPLY" in

# 	    1 ) echo "You picked $opt which is option $REPLY";;
# 	    2 ) echo "You picked $opt which is option $REPLY";;
# 	    3 ) echo "You picked $opt which is option $REPLY";;

# 	    $(( ${#options[@]}+1 )) ) echo "Goodbye!"; break;;
# 	    *) echo "Invalid option. Try another one.";continue;;

# 	    esac

# 	done

# 	while opt=$(zenity --title="$title" --text="$prompt" --list \
# 	                   --column="Options" "${options[@]}"); do

# 	    case "$opt" in
# 	    "${options[0]}" ) zenity --info --text="You picked $opt, option 1";;
# 	    "${options[1]}" ) zenity --info --text="You picked $opt, option 2";;
# 	    "${options[2]}" ) zenity --info --text="You picked $opt, option 3";;
# 	    *) zenity --error --text="Invalid option. Try another one.";;
# 	    esac

# 	done


# }

#__menu

#############
#     1
#############

#__get_bed_files $bed_raw SRR633614 SRR633612
#__get_bed_files $bed_raw SRR633615 SRR633613
#exit

#############
#     2
#############

function __overlap ()
{

	
	#Check if qsub stdout and stderr folders exist else create these folders
	if [ ! -d $BAM/stdout/ ]; then mkdir -p $BAM/stdout/; fi
	if [ ! -d $BAM/stderr/ ]; then mkdir -p $BAM/stderr/; fi

	#INIT qsub options
	opt="-o $BAM/stdout/$day -e $BAM/stdout/$day"
	ovlp_opt="$opt -q batch -l nodes=1:ppn=1,mem=150mb,walltime=01:00:00"

	for bed_in in $BED/*.bed;
	do
		#GET file name
		name=${bed_in##*/}
		sample=${name%.bed}

		#GET file path
		bam_in="$RAW/$sample.sorted.bam"
		bam_out="$BAM/$sample.overlap.bam"

		#GET cmd
		ovlp="$SAMTOOLS view -L $bed_in $bam_in -b > $bam_out"

		#SUBMIT job
		ovlp_jid=$(echo "$ovlp" | qsub $ovlp_opt -N ${sample}.overlap)

	done

}



# __overlap

# exit

#############
#     3
#############



function __sort_by_name ()
{
	#Check if qsub stdout and stderr folders exist else create these folders
	if [ ! -d $BAM/stdout/ ]; then mkdir -p $BAM/stdout/; fi
	if [ ! -d $BAM/stderr/ ]; then mkdir -p $BAM/stderr/; fi

	#INIT qsub options
	opt="-o $BAM/stdout/$day -e $BAM/stderr/$day"
	sort_opt="$opt -q batch -l nodes=1:ppn=8,mem=20gb,walltime=02:00:00"
	rm_opt="$opt -l nodes=1,mem=2gb,walltime=01:00:00"

	for bam_in in $BAM/*.overlap.bam
	do
		#GET file name
		name=${bam_in##*/}
		sample=${name%.overlap.bam}

		#GET file path
		bam_out="$BAM/${sample}.overlap.sorted.bam"

		#GET cmd
		_sort="$SAMTOOLS sort -m 2G -@ 8 -n $bam_in  > $bam_out"
		_rm="rm $bam_in"

		# echo $_sort
		# echo $_rm

		#SUBMIT jobs
		_sort_jid=$(echo "$_sort" | qsub $sort_opt -N sort.${sample})
		_rm_jid=$(echo "$_rm" | qsub $rm_opt -N rm.${sample} -W depend=afterok:$_sort_jid)	
		
	done

}

# __sort_by_name
# exit




function _bam_to_fastq ()
{
	#GET qsub options
	opt="-o $FASTQ/stdout/$day -e $FASTQ/stderr/$day -q batch"
	bam_to_fastq_opt="$opt -l nodes=1:ppn=2,mem=15mb,walltime=00:30:00"

	if [ ! -d $FASTQ/stdout/ ]; then mkdir -p $FASTQ/stdout/; fi
	if [ ! -d $FASTQ/stderr/ ]; then mkdir -p $FASTQ/stderr/; fi

	#GET parameters
	sample=$1

	#GET path in
	bam_in="$BAM/${sample}.overlap.sorted.bam"

	#GET path out
	fastq1="$FASTQ/${sample}.R1.fq"
	fastq2="$FASTQ/${sample}.R2.fq"

	#GET cmd
	bam_to_fastq="$BEDTOOLS bamtofastq -i $bam_in -fq $fastq1 -fq2 $fastq2"

	#SUBMIT cmd
	bam_to_fastq_jid="$(echo "$bam_to_fastq" | qsub $bam_to_fastq_opt -N fastq.${sample})"


	# $(echo "/bioinfo/local/build/BEDTools/BEDTools_2.25.0/bin/bedtools 
	# bamtofastq -i /data/tmp/clancien/Result/Bam/SRR633612.overlap.sorted.bam 
	# -fq /data/tmp/clancien/Result/Bam/SRR633612.R1.fq 
	# -fq2 /data/tmp/clancien/Result/Bam/SRR633612.R2.fq" 
	# | qsub -o /data/tmp/clancien/Result/Bam/out -e /data/tmp/clancien/Result/Bam/err -q batch 
	# -l nodes=1:ppn=8,mem=15mb,walltime=00:30:00 -N fastq)



}
# _bam_to_fastq SRR633612
# _bam_to_fastq SRR633613
# _bam_to_fastq SRR633614
# _bam_to_fastq SRR633615
# exit



function __count_read ()
{
	#GET qsub options
	opt="-o $FASTQ/stdout/$day -e $FASTQ/stderr/$day"
	count_opt="$opt -q batch -l nodes=1:ppn=8,mem=15mb,walltime=00:30:00"

	for fastq in $FASTQ/*R1*
	do
		name=${fastq##*/}
		sample=${name%.R*}

		#GET path out
		COUNTREAD="$FASTQ/${sample}.countread.txt"

		#GET cmd
		count="wc -l $fastq > $COUNTREAD"

		#SUBMIT job

		count_jid="$(echo "$count" | qsub $count_opt -N count.${sample})"


	done
}

# __count_read
# exit

function __merge_count_read ()
{

	COUNTREAD="$FASTQ/countread.txt"

	for count_in in $FASTQ/*.countread.txt
	do
		#GET file name
		name=${count_in##*/}
		sample=${name%.countread.txt}

		#GET cmd
		(echo $sample; cat $count_in) | paste -d',' - - >> $COUNTREAD
		rm $count_in
		#"$(echo "$merge_count")"
		#"$(echo "$rm_count_in")"
	done

}

# __merge_count_read
# echo "end"
# exit

function __downsample ()
{
	#GET qsub options
	opt="-o $DOWNSAMPLING/stdout/$day -e $DOWNSAMPLING/stderr/$day -q batch"
	downsample_opt="$opt -l nodes=1:ppn=1,mem=15mb,walltime=00:45:00"
	rm_opt="$opt -l nodes=1:ppn=1,mem=10mb,walltime=00:10:00"

	if [ ! -d $DOWNSAMPLING/stdout/ ]; then mkdir -p $DOWNSAMPLING/stdout/; fi
	if [ ! -d $DOWNSAMPLING/stderr/ ]; then mkdir -p $DOWNSAMPLING/stderr/; fi

	#GET paramaters
	centromeric_sample=$1
	centromeric_sample_seed_percent=$3
	centromeric_read_number=$5
	centromeric_fastq_size_read_number=$((4*$centromeric_read_number))

	noncentromeric_sample=$2
	noncentromeric_sample_seed_percent=$4
	noncentromeric_read_number=$6
	noncentromeric_fastq_size_read_number=$((4*$noncentromeric_read_number))

	repeat_number=$7

	#GET var
	centromeric_random_seed=$(( RANDOM % 1000000 ))
	noncentromeric_random_seed=$(( RANDOM % 1000000 ))

	#GET path in
	centromeric_sample_R1="$FASTQ/${centromeric_sample}.R1.fq"
	centromeric_sample_R2="$FASTQ/${centromeric_sample}.R2.fq"

	noncentromeric_sample_R1="$FASTQ/${noncentromeric_sample}.R1.fq"
	noncentromeric_sample_R2="$FASTQ/${noncentromeric_sample}.R2.fq"

	#GET path out
	centromeric_sample_R1_tmp_out="$DOWNSAMPLING/${centromeric_sample}.R1.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.fq"
	centromeric_sample_R2_tmp_out="$DOWNSAMPLING/${centromeric_sample}.R2.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.fq"
	
	noncentromeric_sample_R1_tmp_out="$DOWNSAMPLING/${noncentromeric_sample}.R1.${repeat_number}.${noncentromeric_read_number}.${noncentromeric_random_seed}.fq"
	noncentromeric_sample_R2_tmp_out="$DOWNSAMPLING/${noncentromeric_sample}.R2.${repeat_number}.${noncentromeric_read_number}.${noncentromeric_random_seed}.fq"

	cat_centromeric_sample_R1_tmp_out_noncentromeric_sample_R1_tmp_out="$DOWNSAMPLING/${centromeric_sample}.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.${noncentromeric_random_seed}.fq"
	cat_centromeric_sample_R2_tmp_out_noncentromeric_sample_R2_tmp_out="$DOWNSAMPLING/${centromeric_sample}.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number}.${centromeric_random_seed}.${noncentromeric_random_seed}.fq"

	#GET cmd
	centromeric_downsample_R1="$SEQTK sample -s $centromeric_random_seed $centromeric_sample_R1 $centromeric_sample_seed_percent | head -n $centromeric_fastq_size_read_number > $centromeric_sample_R1_tmp_out"
	centromeric_downsample_R2="$SEQTK sample -s $centromeric_random_seed $centromeric_sample_R2 $centromeric_sample_seed_percent | head -n $centromeric_fastq_size_read_number > $centromeric_sample_R2_tmp_out"
	
	noncentromeric_downsample_R1="$SEQTK sample -s $noncentromeric_random_seed $noncentromeric_sample_R1 $noncentromeric_sample_seed_percent | head -n $noncentromeric_fastq_size_read_number > $noncentromeric_sample_R1_tmp_out"
	noncentromeric_downsample_R2="$SEQTK sample -s $noncentromeric_random_seed $noncentromeric_sample_R2 $noncentromeric_sample_seed_percent | head -n $noncentromeric_fastq_size_read_number > $noncentromeric_sample_R2_tmp_out"

	cat_centromeric_downsample_R1_noncentromeric_downsample_R1="cat $noncentromeric_sample_R1_tmp_out $centromeric_sample_R1_tmp_out > $cat_centromeric_sample_R1_tmp_out_noncentromeric_sample_R1_tmp_out"
	cat_centromeric_downsample_R2_noncentromeric_downsample_R2="cat $noncentromeric_sample_R2_tmp_out $centromeric_sample_R2_tmp_out  > $cat_centromeric_sample_R2_tmp_out_noncentromeric_sample_R2_tmp_out"

	rm_tmp_R1="rm $centromeric_sample_R1_tmp_out $noncentromeric_sample_R1_tmp_out"
	rm_tmp_R2="rm $centromeric_sample_R2_tmp_out $noncentromeric_sample_R2_tmp_out"






	#SUBMIT job
	centromeric_downsample_R1_jid="$(echo "$centromeric_downsample_R1" | qsub $downsample_opt -N downsample.${centromeric_sample}.R1.${repeat_number}.${centromeric_read_number})"
	centromeric_downsample_R2_jid="$(echo "$centromeric_downsample_R2" | qsub $downsample_opt -N downsample.${centromeric_sample}.R2.${repeat_number}.${centromeric_read_number})"

	noncentromeric_downsample_R1_jid="$(echo "$noncentromeric_downsample_R1" | qsub $downsample_opt -N downsample.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number})"
	noncentromeric_downsample_R2_jid="$(echo "$noncentromeric_downsample_R2" | qsub $downsample_opt -N downsample.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number})"

	cat_centromeric_downsample_R1_noncentromeric_downsample_R1_jid="$(echo "$cat_centromeric_downsample_R1_noncentromeric_downsample_R1" | qsub $downsample_opt -N cat.${centromeric_sample}.${noncentromeric_sample}.R1.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R1_jid:$noncentromeric_downsample_R1_jid )"
	cat_centromeric_downsample_R2_noncentromeric_downsample_R2_jid="$(echo "$cat_centromeric_downsample_R2_noncentromeric_downsample_R2" | qsub $downsample_opt -N cat.${centromeric_sample}.${noncentromeric_sample}.R2.${repeat_number}.${centromeric_read_number} -W depend=afterok:$centromeric_downsample_R2_jid:$noncentromeric_downsample_R2_jid)"

	rm_tmp_R1_jid="$(echo "$rm_tmp_R1" | qsub $rm_opt -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R1.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R1_noncentromeric_downsample_R1_jid)"
	rm_tmp_R2_jid="$(echo "$rm_tmp_R2" | qsub $rm_opt -N rm.${centromeric_sample}.${noncentromeric_sample}.tmp.R2.${centromeric_read_number} -W depend=afterok:$cat_centromeric_downsample_R2_noncentromeric_downsample_R2_jid)"

	# opt="-o $BAM/stdout/$day -e $BAM/stdout/$day"
	# sub_opt="$opt -q batch -l nodes=1:ppn=1,mem=150mb,walltime=01:00:00"



# readonly SEQTK="/data/tmp/clancien/Tools/seqtk/seqtk"

# random_seed=echo $(( RANDOM % 100000 ))


# R1="$BAM/SRR633614.R1.fq"
# R2="$BAM/SRR633614.R2.fq"

# job_id1="$(echo "$SEQTK sample -s $random_seed $R1 0.015 | head -n 4000000 > $BAM/SRR633614.R1.random.sub.fq" | qsub $sub_opt -N sub.R1)"
# job_id2="$(echo "$SEQTK sample -s $random_seed $R2 0.015 | head -n 4000000 > $BAM/SRR633614.R2.random.sub.fq" | qsub $sub_opt -N sub.R2)"


}





# __downsample SRR633614 SRR633612 0.15 0.96 900000 89100000 1
# exit


function __get_read_number ()
{
	awk -v sample_name="$1" -F, '
		{
			if($1==sample_name)
				{
					print $2/4
				}
		}
	' $FASTQ/countread.txt
}


function __get_seed_percent ()
{

	echo "$(echo $1/$2 | bc -l) "

}
function __downsample_auto ()
{


	if [ ! -d "$DOWNSAMPLING" ]; then mkdir -p "$DOWNSAMPLING"; fi
	#j=0
	centromeric_sample=$1
	noncentromeric_sample=$2
	centromeric_read_number=$(__get_read_number $centromeric_sample)
	noncentromeric_read_number=$(__get_read_number $noncentromeric_sample)

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		#echo $sample_centromeric_size
		sample_noncentromeric_size=$(( 90000000 - $sample_centromeric_size ))
		centromeric_seed_percent=$(__get_seed_percent $sample_centromeric_size $noncentromeric_read_number)
		noncentromeric_seed_percent=$(__get_seed_percent $sample_noncentromeric_size $noncentromeric_read_number)

		# echo "$centromeric_seed_percent"
		# echo "$noncentromeric_seed_percent"
		# echo "-----------------"



		for repeat in $(seq 1 $4)
		do

			__downsample $centromeric_sample $noncentromeric_sample $centromeric_seed_percent $noncentromeric_seed_percent $sample_centromeric_size $sample_noncentromeric_size $repeat
			#exit

			# # __downsampling SRR633614 SRR633612 1000000 99000000 1
			# __downsampling $centromeric_sample $noncentromeric_sample $sample_centromeric_size $sample_noncentromeric_size $repeat

			# __downsample SRR633614 SRR633612 0.15 0.96 1000000 89000000
			# #let j=j+7
		done

	done
}


# __downsample_auto SRR633614 SRR633612 10 3
# __downsample_auto SRR633615 SRR633613 10 3
# exit



function __align ()
{
	#GET qsub options
	opt="-o $ALIGN/stdout/$day -e $ALIGN/stderr/$day -q batch"

	aln_opt="$opt -l nodes=1:ppn=8,mem=16gb,walltime=24:00:00"

	sort_opt="$opt -l nodes=1:ppn=8,mem=20gb,walltime=04:00:00"

	index_opt="$opt -l nodes=1,mem=2gb,walltime=04:00:00"

	rm_opt="$opt -l nodes=1,mem=2gb,walltime=01:00:00"

	#CHECK if qsub stdout and stderr folders exist else create them
	if [ ! -d "$ALIGN/stdout/$day" ]; then mkdir -p "$ALIGN/stdout/$day"; fi
	if [ ! -d "$ALIGN/stderr/$day" ]; then mkdir -p "$ALIGN/stderr/$day"; fi

	#GET Parameters
	
	R1=$1
	R2=$2

	centromeric_sample=$3
	noncentromeric_sample=$4

	repeat=$5
	sample_centromeric_size=$6

	aln_out=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${sample_centromeric_size}.bam
	sort_out=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${sample_centromeric_size}.sorted.bam

	#GET cmd
	_aln="$BOWTIE2 --threads 8 --very-sensitive -x $INDEX -1 ${R1} -2 ${R2} | $SAMTOOLS view -b - > $aln_out"
	_sort="$SAMTOOLS sort -m 2G -@ 8 $aln_out $sort_out"
	_index="$SAMTOOLS index -b $sort_out"
	_rm="rm aln_out"

	#SUBMIT jobs
	aln_jid="$(echo $_aln | qsub $aln_opt -N aln.${1}.${2}.${repeat}.${sample_centromeric_size})"
	sort_jid="$(echo $_sort | qsub $sort_opt -N sort.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$aln_jid)"
	index_jid="$(echo $_index | qsub $index_opt -N index.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"
	rm_jid="$(echo $_rm | qsub $rm_opt -N rm.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"


}

function __align_all ()
{
	if [ ! -d "$ALIGN" ]; then mkdir -p "$ALIGN"; fi


	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		#echo $sample_centromeric_size
		#echo $sample_noncentromeric_size
		for repeat in $(seq 1 $4)
		do
			seed1="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 6)"
			seed2="$(ls $DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.*.fq | cut -d'.' -f 7)"

			R1_path=$DOWNSAMPLING/${1}.${2}.R1.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq
			R2_path=$DOWNSAMPLING/${1}.${2}.R2.${repeat}.${sample_centromeric_size}.${seed1}.${seed2}.fq

			__align $R1_path $R2_path $1 $2 $repeat $sample_centromeric_size
			#echo $R1
			# __downsampling SRR633614 SRR633612 1000000 99000000 1
			#__align SRR633614.SRR633612.subsample.1000000.R1.fastq SRR633614.SRR633612.subsample.1000000.R2.fastq

			#__align $1 $2 $sample_centromeric_size $sample_noncentromeric_size $repeat
		done

	done
}

# __align_all SRR633614 SRR633612 10 3
# __align_all SRR633615 SRR633613 10 3
#exit


function __number_of_alignement()
{
	#INIT qsub options
	opt="-o $ALIGN/stdout/day -e $ALIGN/stdout/day"
	count_opt="$opt -q batch -l nodes=1:ppn=8,mem=15mb,walltime=00:30:00"

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"
			file_out="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.count.txt"
			
			#GET job
			count="$SAMTOOLS view $bam_in | wc -l > $file_out"

			#SUBMIT job
			count_jid="$(echo "$count" | qsub $opt -N "$file_out")"
		done
	done
}

# __number_of_alignement SRR633614 SRR633612 10 3
# __number_of_alignement SRR633615 SRR633613 10 3
# exit

function __merge_number_of_alignement()
{
	file_out="$ALIGN/align.count.txt"
	$(echo "filename,number_of_read" > "$file_out")

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			file_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.count.txt"
			number=$(head -n 1 $file_in)

			$(echo "${file_in},${number}" >> "$file_out")
		done
	done

}

__merge_number_of_alignement SRR633614 SRR633612 10 3
#__merge_number_of_alignement SRR633615 SRR633613 10 3
exit


function _overlap_cytoband() #centromeric #noncentromeric ##percent #repeat 
{
	#Check if qsub stdout and stderr folders exist else create these folders
	if [ ! -d $OVERLAP/stdout/ ]; then mkdir -p $OVERLAP/stdout/; fi
	if [ ! -d $OVERLAP/stderr/ ]; then mkdir -p $OVERLAP/stderr/; fi

	#INIT qsub options
	opt="-o $OVERLAP/stdout/$day -e $OVERLAP/stdout/$day"
	ovlpf_opt_centromeric="$opt -q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:20:00"
	ovlpf_opt_noncentromeric="$opt -q batch -l nodes=1:ppn=1,mem=15mb,walltime=02:00:00"
	centromeric_bed=$CYTOBAND/centromeric.bed
	noncentromeric_bed=$CYTOBAND/noncentromeric.bed

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		#echo $sample_centromeric_size
		#echo $sample_noncentromeric_size
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			ovlpf_centromeric="$SAMTOOLS view -L $centromeric_bed $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric"
			ovlpf_noncentromeric="$SAMTOOLS view -L $noncentromeric_bed $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric"

			#SUBMIT job
			ovlpf_centromeric_jid="$(echo "$ovlpf_centromeric" | qsub $ovlpf_opt_centromeric -N centromeric.test )"
			ovlpf_noncentromeric_jid="$(echo "$ovlpf_noncentromeric" | qsub $ovlpf_opt_noncentromeric -N noncentromeric.test)"
			#exit # First : /data/tmp/clancien/Result/Align/SRR633614.SRR633612.subsample.1.900000.bam
		done

	done

}
# _overlap_cytoband SRR633615 SRR633613 10 3
# exit


function __merge_overlap_cytoband()
{
	centromeric_sample=$1
	noncentromeric_sample=$2
	total=$(($3*900000))
	file_out=$ALIGN/test2.txt
	echo $file_out
	if [ -f $file_out ]; then rm $file_out; fi;


	echo "centromeric_sample,noncentromeric_sample,read_number,subsample,centromeric_number,noncentromeric_number,repeat,observed_ratio,theorical_ratio" >> $file_out

	for((i=900000;i<=$total;i+=900000))
	do
		for repeat in $(seq 1 $4)
		do
			centromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.centromeric
			noncentromeric_path=$ALIGN/${centromeric_sample}.${noncentromeric_sample}.subsample.${repeat}.${i}.noncentromeric

			centromeric_number=$(head -n 1 $centromeric_path)
			noncentromeric_number=$(head -n 1 $noncentromeric_path)

			observed_ratio=$(echo $centromeric_number/$noncentromeric_number | bc -l)
			theorical_ratio=$(echo $i/90000000 | bc -l)


			echo ${centromeric_sample},${noncentromeric_sample},${i},${repeat},${centromeric_number},${noncentromeric_number},${repeat},${observed_ratio},${theorical_ratio} >> $file_out

		done
	done
}
# __merge_overlap_cytoband SRR633614 SRR633612 10 3
# __merge_overlap_cytoband SRR633615 SRR633613 10 2
# exit


function _overlap_cytoband() #centromeric #noncentromeric ##percent #repeat 
{
	#Check if qsub stdout and stderr folders exist else create these folders
	if [ ! -d $OVERLAP/stdout/ ]; then mkdir -p $OVERLAP/stdout/; fi
	if [ ! -d $OVERLAP/stderr/ ]; then mkdir -p $OVERLAP/stderr/; fi

	#INIT qsub options
	opt="-o $OVERLAP/stdout/$day -e $OVERLAP/stdout/$day"
	ovlpf_opt_centromeric="$opt -q batch -l nodes=1:ppn=1,mem=15mb,walltime=00:20:00"
	ovlpf_opt_noncentromeric="$opt -q batch -l nodes=1:ppn=1,mem=15mb,walltime=02:00:00"
	centromeric_bed=$CYTOBAND/centromeric.bed
	noncentromeric_bed=$CYTOBAND/noncentromeric.bed

	for ((sample_centromeric_size=900000;sample_centromeric_size<=$3*900000;sample_centromeric_size+=900000))
	do
		#echo $sample_centromeric_size
		#echo $sample_noncentromeric_size
		for repeat in $(seq 1 $4)
		do
			#GET PATH IN
			bam_in="$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam"

			#GET cmd
			ovlpf_centromeric="$SAMTOOLS view -L $centromeric_bed $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.centromeric"
			ovlpf_noncentromeric="$SAMTOOLS view -L $noncentromeric_bed $bam_in -b | $SAMTOOLS view | wc -l > $ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.noncentromeric"

			#SUBMIT job
			ovlpf_centromeric_jid="$(echo "$ovlpf_centromeric" | qsub $ovlpf_opt_centromeric -N centromeric.test )"
			ovlpf_noncentromeric_jid="$(echo "$ovlpf_noncentromeric" | qsub $ovlpf_opt_noncentromeric -N noncentromeric.test)"
			#exit # First : /data/tmp/clancien/Result/Align/SRR633614.SRR633612.subsample.1.900000.bam
		done

	done

}







# __align_all SRR633614 SRR633612 10 3
# __align_all SRR633615 SRR633613 10 3

##Only to test
# function __get_mapped_read_only () 
# {
# 	#GET qsub options
# 	option="-o $BAM/stdout/$day -e $BAM/stderr/$day -q batch"
# 	option_get_mapped_read_only=" $option -l nodes=1:ppn=8,mem=15mb,walltime=00:30:00"

# 	#GET parameters
# 	study=$1
	

# 	#GET path in
# 	bam_in=$BAM/${study}.overlap.sorted.bam

# 	#GET path out
# 	bam_out=$BAM/${study}.overlap.sorted.mapped.read.only.bam

# 	#GET cmd
# 	mapped_read_only="$SAMTOOLS view -f 0x02 $bam_in -b > $bam_out"

# 	# -f 0x02 -F 0x04 

# 	#SUBMIT job
# 	mapped_read_only="$(echo "$mapped_read_only" | qsub $option_get_mapped_read_only -N mapped_read_only.${study})"

# 	# $(echo "/bioinfo/local/build/samtools/samtools-1.3/bin/samtools view -f 0x02 /data/tmp/clancien/Result/Bam/SRR633612.overlap.sorted.bam -b >
# 	#  /data/tmp/clancien/Result/Bam/SRR633612.mapped.read.only.bam" |
# 	#   qsub   -N mapped)
# 	# samtools view -f 0x02 -b in.bam > out.aligned.bam


# }

# __get_mapped_read_only SRR633612
# __get_mapped_read_only SRR633613
# __get_mapped_read_only SRR633614
# __get_mapped_read_only SRR633615









# __/\\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\_____/\\\\\\\\\\\____/\\\\\\\\\\\\\\\_        
#  _\///////\\\/////__\/\\\///////////____/\\\/////////\\\_\///////\\\/////__       
#   _______\/\\\_______\/\\\______________\//\\\______\///________\/\\\_______      
#    _______\/\\\_______\/\\\\\\\\\\\_______\////\\\_______________\/\\\_______     
#     _______\/\\\_______\/\\\///////___________\////\\\____________\/\\\_______    
#      _______\/\\\_______\/\\\_____________________\////\\\_________\/\\\_______   
#       _______\/\\\_______\/\\\______________/\\\______\//\\\________\/\\\_______  
#        _______\/\\\_______\/\\\\\\\\\\\\\\\_\///\\\\\\\\\\\/_________\/\\\_______ 
#         _______\///________\///////////////____\///////////___________\///________

# Start test

function _downsample ()
{
	#sample=$1
	#rawdir=$2
	#outdir=$3

	R1="$BAM/SRR633614.R1.fq"
	R2="$BAM/SRR633614.R2.fq"

	tmp="$BAM/SRR633614.tmp.fastq"

	if [ ! -f $tmp ]; then

	    paste <($R1) <($R2) | less
	    exit
	    # awk '{printf("%s", $0); n++; if (n%4 == 0) printf("\n"); else printf("\t\t");}' > $tmp

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

}

# _downsample
# exit
# /data/tmp/clancien/Result/Bam

# SRR633614.R1.fq  SRR633614.R2.fq



# End test


















function __count_read ()
{
	#INIT qsub options
	opt="-o $BAM/stdout/day -e $BAM/stdout/day"
	count_opt="$opt -q batch -l nodes=1:ppn=8,mem=15mb,walltime=00:30:00"

	for sort_in in $BAM/*.overlap.sorted.bam
	do
		#GET file name
		name=${sort_in##*/}
		sample=${name%.overlap.sorted.bam}

		#GET file path
		count_out="$BAM/${sample}.countread.txt"

		#GET cmd
		count_uniq="$SAMTOOLS view $sort_in | cut -f 1 | uniq | wc -l >> $count_out"
		count_all="$SAMTOOLS view $sort_in | wc -l >> $count_out"

		#echo "$count_uniq"
		#echo "$count_all"
		#exit
		#SUBMIT jobs
		count_uniq_jid="$(echo "$count_uniq" | qsub $count_opt -N count_uniq.${sample})"
		count_all_jid="$(echo "$count_all" | qsub $count_opt -N count_all.${sample} -W depend=afterok:$count_uniq_jid)"

	done
}

# __count_read
# exit

function __merge_count_read ()
{

	##ADD header in countread

	#GET cmd
	add_header="echo sample,uniq_read_number,read_number >> $merge_count_out"

	#SUBMIT cmd
	add_header_id="$( eval $add_header )"


	##MERGE all ${sample}.countread.txt in countread.txt then
	##REMOVE all ${sample}.countread.txt

	for count_in in $BAM/*.countread.txt
	do
		#GET file name
		name=${count_in##*/}
		sample=${name%.countread.txt}

		#GET cmd
		merge_count="(echo $sample; cat $count_in) | paste -d',' - - - >> $merge_count_out"
		rm_count_in="rm $count_in"
		
		#SUBMIT cmd
		merge_count_id="$( eval $merge_count )"

		rm_count_in_id="$( eval $rm_count_in )"
		#echo "$(eval $merge_count)"
		#(echo $sample ; $merge_count;)| paste -d',' - - -
		#(echo $sample ; $merge_count;) | paste -d',' - - - >>$merge_count_out
		##{ command1 & command2; }


	done

}

#__merge_count_read
exit

function __get_read_number ()
{
	tail -n+2 $merge_count_out |  awk -v col_name="$1" -v sample_name="$2" -F, '
		{
			if($1==sample_name)
				{
					if(col_name=="uniq_read_number")
						{
							print $2
						}
					else if(col_name=="read_number")
						{
							print $3
						}
					else
						{
							print "Error col_name"
							print "Choose uniq_read_number or read_number"
							exit 1
						}
				}
		}
	'
}

#__get_read_number uniq_read_number SRR633613


function __get_centromeric_seed_percent ()
{
	# #Explain the last line

	# #GET parameters
	# centromeric_sample=$1
	# sample_centromeric_size=$2
	# sample_centromeric_read_number=$(__get_read_number read_number "$centromeric_sample")

	# #echo $sample_read_number

	# #GET cmd
	# centromeric_seed_percent="$(echo $sample_centromeric_size/$sample_centromeric_read_number | bc -l)"

	# #RETURN centromeric_seed_percent value
	# echo "$centromeric_seed_percent"

	echo "$(echo $2/$(__get_read_number read_number "$1") | bc -l)"

}

function __get_noncentromeric_seed_percent ()
{
	# #Explain the last line
	
	# #GET parameters
	# noncentromeric_sample=$1
	# sample_noncentromeric_size=$2
	# sample_noncentromeric_read_number=$(__get_read_number read_number "$noncentromeric_sample")

	# #echo $sample_read_number

	# #GET cmd
	# noncentromeric_seed_percent="$(echo $sample_noncentromeric_size/$sample_noncentromeric_read_number | bc -l)"

	# #RETURN centromeric_seed_percent value
	# echo "$noncentromeric_seed_percent"

	echo "$(echo $2/$(__get_read_number read_number "$1") | bc -l) "
}

function __get_bam_header_length ()
{
	#$1 bam file path
	#$2 search headers only in #$2 first line --> faster
	echo "$($SAMTOOLS view -h $1 | head -n $2 | grep "^@" | wc -l)"

}

function __downsampling ()
{
	## Need 3 paramters :
	# size of final file (i.e. 100 M reads)
	# centromeric read number (i.e. 1M reads)
	# sample number (i.e. repeat 3 times)

	#GET Parameters
	centromeric_sample=$1
	noncentromeric_sample=$2
	sample_centromeric_size=$3
	sample_noncentromeric_size=$4
	sample_number=$5

	#GET seed percent
	centromeric_seed_percent=$(__get_centromeric_seed_percent "$centromeric_sample" "$sample_centromeric_size")	

	noncentromeric_seed_percent=$(__get_noncentromeric_seed_percent "$noncentromeric_sample" "$sample_noncentromeric_size")


	#INIT qsub options
	opt="-o $DOWNSAMPLING/stdout/$day -e $$DOWNSAMPLING/stderr/$day -q batch"
	subsample_opt="$opt -l nodes=1:ppn=1,mem=20mb,walltime=01:00:00"
	cat_opt="$opt -l nodes=1:ppn=1,mem=100mb,walltime=00:15:00"
	rm_opt="$opt -l nodes=1:ppn=1,mem=300mb,walltime=00:15:00"

	#CHECK if qsub stdout and stderr folders exist else create them

	if [ ! -d "$DOWNSAMPLING/stdout/$day" ]; then mkdir -p "$DOWNSAMPLING/stdout/$day"; fi
	if [ ! -d "$DOWNSAMPLING/stderr/$day" ]; then mkdir -p "$DOWNSAMPLING/stderr/$day"; fi

	#GET PATH IN
	centromeric_sample_file=$BAM/${centromeric_sample}.overlap.sorted.bam
	noncentromeric_sample_file=$BAM/${noncentromeric_sample}.overlap.sorted.bam

	#GET PATH OUT
	centromeric_out=$DOWNSAMPLING/${centromeric_sample}.subsample.${sample_number}.${sample_centromeric_size}.bam
	noncentromeric_out=$DOWNSAMPLING/${noncentromeric_sample}.subsample.${sample_number}.${sample_noncentromeric_size}.bam
	cat_out=$DOWNSAMPLING/${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_number}.${sample_centromeric_size}.bam
	fastq_R1_out=$DOWNSAMPLING/${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_number}.${sample_centromeric_size}.R1.fastq
	fastq_R2_out=$DOWNSAMPLING/${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_number}.${sample_centromeric_size}.R2.fastq
	
	#GET Head line number
	centromeric_header_length=$(__get_bam_header_length "$centromeric_sample_file" 10000)
	noncentromeric_header_length=$(__get_bam_header_length "$noncentromeric_sample_file" 10000)

	#GET Total size = size sample + size header
	centromeric_total_size=$(( $sample_centromeric_size + $centromeric_header_length ))
	noncentromeric_total_size=$(( $sample_noncentromeric_size + $noncentromeric_header_length ))

	#GET cmd
	subsample_centromeric="$SAMTOOLS view -h -f 0x02 -s $centromeric_seed_percent $centromeric_sample_file -b | $SAMTOOLS view -h | head -n $centromeric_total_size | $SAMTOOLS view -b > $centromeric_out "
	subsample_noncentromeric="$SAMTOOLS view -h -f 0x02 -s $noncentromeric_seed_percent $noncentromeric_sample_file -b | $SAMTOOLS view -h | head -n $noncentromeric_total_size | $SAMTOOLS view -b > $noncentromeric_out "
	subsample_cat="$SAMTOOLS cat -o $cat_out $centromeric_out $noncentromeric_out"
	subsample_fastq="$SAMTOOLS fastq $cat_out  -1 $fastq_R1_out -2 $fastq_R2_out"
	
	rm_subsample_centromeric="rm $centromeric_out"
	rm_subsample_noncentromeric="rm $noncentromeric_out"
	rm_subsample_cat="rm $cat_out"

	#SUBMIT job
	subsample_centromeric_jid="$(echo $subsample_centromeric | qsub $subsample_opt -N seed.${centromeric_sample}.reads.${sample_centromeric_size} )"
	subsample_noncentromeric_jid="$(echo $subsample_noncentromeric | qsub $subsample_opt -N seed.${noncentromeric_sample}.reads.${sample_noncentromeric_size} )"
	subsample_cat_jid="$(echo $subsample_cat |  qsub $cat_opt -N cat.${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_centromeric_size}  -W depend=afterok:$subsample_centromeric_jid:$subsample_noncentromeric_jid )"
	subsample_fastq_jid="$(echo $subsample_fastq | qsub $cat_opt -N fastq.${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_centromeric_size} -W depend=afterok:$subsample_cat_jid )"

	rm_subsample_centromeric_jid="$(echo $rm_subsample_centromeric | qsub $rm_opt -N rm.${centromeric_sample}.reads.${sample_centromeric_size} -W depend=afterok:$subsample_cat_jid )"
	rm_subsample_noncentromeric_jid="$(echo $rm_subsample_noncentromeric | qsub $rm_opt -N rm.${noncentromeric_sample}.reads.${sample_noncentromeric_size} -W depend=afterok:$subsample_cat_jid )"
	rm_subsample_cat_jid="$(echo $rm_subsample_cat | qsub $rm_opt -N rm.${centromeric_sample}.${noncentromeric_sample}.subsample.${sample_centromeric_size} -W depend=afterok:$subsample_fastq_jid )"

}

# __downsampling SRR633614 SRR633612 1000000 99000000 1
# exit

function __downsampling_auto ()
{

	if [ ! -d "$DOWNSAMPLING" ]; then mkdir -p "$DOWNSAMPLING"; fi
	#j=0
	centromeric_sample=$1
	noncentromeric_sample=$2
	for ((sample_centromeric_size=1000000;sample_centromeric_size<=$3*1000000;sample_centromeric_size+=1000000))
	do
		#echo $sample_centromeric_size
		sample_noncentromeric_size=$(( 100000000 - $sample_centromeric_size ))
		#echo $sample_noncentromeric_size
		for repeat in $(seq 1 $4)
		do
			# __downsampling SRR633614 SRR633612 1000000 99000000 1
			__downsampling $centromeric_sample $noncentromeric_sample $sample_centromeric_size $sample_noncentromeric_size $repeat
			#let j=j+7
		done

	done
	#echo $j
}

# __downsampling_auto SRR633614 SRR633612 10 3
# __downsampling_auto SRR633615 SRR633613 10 3

# exit

function __align ()
{
	#GET qsub options
	opt="-o $ALIGN/stdout/$day -e $ALIGN/stderr/$day -q batch"

	aln_opt="$opt -l nodes=1:ppn=8,mem=16gb,walltime=24:00:00"

	sort_opt="$opt -l nodes=1:ppn=8,mem=20gb,walltime=04:00:00"

	index_opt="$opt -l nodes=1,mem=2gb,walltime=04:00:00"

	rm_opt="$opt -l nodes=1,mem=2gb,walltime=01:00:00"

	#CHECK if qsub stdout and stderr folders exist else create them
	if [ ! -d "$ALIGN/stdout/$day" ]; then mkdir -p "$ALIGN/stdout/$day"; fi
	if [ ! -d "$ALIGN/stderr/$day" ]; then mkdir -p "$ALIGN/stderr/$day"; fi

	#GET Parameters
	R1=$DOWNSAMPLING/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.R1.fastq
	R2=$DOWNSAMPLING/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.R2.fastq

	aln_out=$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.bam
	sort_out=$ALIGN/${1}.${2}.subsample.${repeat}.${sample_centromeric_size}.sorted.bam

	#GET cmd
	_aln="$BOWTIE2 --threads 8 --very-sensitive -x $INDEX -1 ${R1} -2 ${R2} | $SAMTOOLS view -b - > $aln_out"
	_sort="$SAMTOOLS sort -m 2G -@ 8 $aln_out $sort_out"
	_index="$SAMTOOLS index -b $sort_out"
	_rm="rm aln_out"

	#SUBMIT jobs
	aln_jid="$(echo $_aln | qsub $aln_opt -N aln.${1}.${2}.subsample.${repeat}.${sample_centromeric_size})"
	sort_jid="$(echo $_sort | qsub $sort_opt -N sort.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$aln_jid)"
	index_jid="$(echo $_index | qsub $index_opt -N index.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"
	rm_jid="$(echo $_rm | qsub $rm_opt -N rm.${1}.${2}.subsample.${repeat}.${sample_centromeric_size} -W depend=afterok:$sort_jid)"


	# exit
	# R1=$1
	# R2=$2
	# base_name=${R1%.R*}
	# name=$ALIGN/${base_name}.subsample.${3}.${4}
	# echo $name
	# exit
	# aln="$BOWTIE2 --threads 8 --very-sensitive -x $INDEX -1 ${R1} -2 ${R2}"

	

	# alnbam="$aln | $SAMTOOLS view -b - > $name.aln.bam"
	# #echo "$alnbam"
	# sortbam="$SAMTOOLS sort -m 2G -@ 8 $name.aln.bam $name.aln.sorted"
	# indexbam="$SAMTOOLS index -b $name.aln.sorted.bam"
	# rmbam="rm $name.aln.bam"
	# echo "$alnbam"
	# echo "$sortbam"
	# echo "$indexbam"


	# aln_jid=`echo "$alnbam" | qsub $aln_opt -N ${name}.aln`
                   
 #    sort_jid=`echo "$sortbam" | qsub $sort_opt -N ${name}.sort -W depend=afterok:$aln_jid`
        
 #    index_jid=`echo "$indexbam" | qsub $index_opt -N ${name}.index -W depend=afterok:$sort_jid` 

 #    rm_jid=`echo "$rmbam" | qsub $rm_opt -N ${name}.rm -W depend=afterok:$index_jid`




	#SRR633614.SRR633612.subsample.1000000.R2.fastq
	#aln="$bowtie2 --threads 8 --very-sensitive -x $index -1 ${centromericRead[i]} -2 ${centromericRead[i+1]}"

}

function __align_all ()
{
	if [ ! -d "$ALIGN" ]; then mkdir -p "$ALIGN"; fi


	for ((sample_centromeric_size=1000000;sample_centromeric_size<=$3*1000000;sample_centromeric_size+=1000000))
	do
		#echo $sample_centromeric_size
		sample_noncentromeric_size=$(( 100000000 - $sample_centromeric_size ))
		#echo $sample_noncentromeric_size
		for repeat in $(seq 1 $4)
		do

			#echo $R1
			#exit
			# __downsampling SRR633614 SRR633612 1000000 99000000 1
			#__align SRR633614.SRR633612.subsample.1000000.R1.fastq SRR633614.SRR633612.subsample.1000000.R2.fastq

			__align $1 $2 $sample_centromeric_size $sample_noncentromeric_size $repeat
		done

	done
}

# __align_all SRR633614 SRR633612 10 3
# __align_all SRR633615 SRR633613 10 3
exit
#SRR633614.SRR633612.subsample.2.2000000.R1.fastq

#__align SRR633614.SRR633612.subsample.1000000.R1.fastq SRR633614.SRR633612.subsample.1000000.R2.fastq
#create __align_all
#loop over 1 000 000..10 000 000 --> X
#loop over 1..3
#call two times __align_all :
#__align SRR633614.SRR633612.subsample.Y.X.R1.fastq SRR633614.SRR633612.subsample.Y.X.R2.fastq 10 3
#__align SRR633615.SRR633613.subsample.Y.X.R1.fastq SRR633615.SRR633613.subsample.Y.X.R2.fastq 10 3
#countread_file=$BASE_DIR/test4/countread.txt #$1

## concat all count files








__getReadNumber()
{
	#printf "%s\n" "here"

	tail -n+2 $1 | while read line 
	do
		local filename=$(echo $line | cut -d',' -f 1)
		#printf "%s\t%s\n" "$filename" "$2"
		if [ $filename == "$2" ]; then
			echo $(echo $line | cut -d ',' -f 2)
			break
		fi
	done
}

__downsampling()

{
	centromeric_path="$BASE_DIR/test4/SRR633612.noncentromeric.1.overlap.sorted.bam"
	seed_centromeric=$(echo "1000000/$(__getReadNumber $countread_file "SRR633612.noncentromeric.1")" | bc -l)  #$(( $(__getReadNumber $countread_file "SRR633612.noncentromeric.1") / 100000000))

	downsample_centromere="samtools view -s $centromeric_seed -b $centromeric_path > $centromeric_path.downsample.bam"
	job_id=`samtools view -s $centromeric_seed -b $centromeric_path > $centromeric_path.downsample.bam` 

	#seed_noncentromeric=$(__getReadNumber $countread_file "SRR633614.centromeric.1")

	echo $seed_centromeric
	echo $seed_noncentromeric

}

#__downsampling

# SAMTOOLS=/bioinfo/local/build/samtools/samtools-1.3/bin/samtools
# OPTION=" -o $BASE_DIR/Result/Bam/stdout/${DATE} -e $BASE_DIR/Result/Bam/stdout/${DATE}"
# PBS="$OPTION -q batch -l nodes=1:ppn=1,mem=10mb,walltime=00:30:00"

