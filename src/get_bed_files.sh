#!/usr/bin/env bash

#
# Purpose: Create bed file.
#
function __get_bed_file()
{
	# $1 : CSV File (binned) - name file
	# $2 : CenH3 enrichment - sample name
	# $3 : WT - sample name
	# $4 : output : centromeric bed file - name file
	# $5 : output : non centromeric bed file - name file

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

function myfunc()
{
    local  __resultvar=$1
    local  myresult='some value'
    if [[ "$__resultvar" ]]; then
        eval $__resultvar="'$myresult'"
    else
        echo "$myresult"
    fi
}

myfunc result
echo $result
result2=$(myfunc)
echo $result2

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

	# Store headers in an array
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
