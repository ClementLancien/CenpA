#!/usr/bin/env bash
#==> MAKEFILE ??????

#exit as soon as any line in the script fails
set -e

#exit as soon as any variable is not init
set -o nounset #equals set -u

#./get_ALR_RepBase.sh < ../RepBase23.01.fasta/humrep.ref > output.fa
bool=0

while read line
do
	#echo $line
	#echo ${line:0:1}
	if [[ ${line:0:1} = ">" ]]; then
		#echo "here1"
		if [[ $line =~ .*">ALR".* ]] &&	 [[ "$line" =~ .*"Homo sapiens".* ]]; then
			bool=1
			echo $line

		else
			bool=0
		fi
	else
		#echo $boolean
		if [ $bool -eq 1 ]; then
			echo $line
		fi
	fi
done 