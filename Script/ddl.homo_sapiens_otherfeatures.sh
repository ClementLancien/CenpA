#!/usr/bin/env bash

# exit when any command fails
set -e
set -o nounset

#PATHS
readonly CWD=$(pwd)
readonly BASE_DIR=$(dirname "$CWD")
readonly FEATURES=${BASE_DIR}/Features/

#Check if folder exists else create it
if [ ! -d "$FEATURES" ]; then mkdir -p $FEATURES ;fi

#Download - Move to folder - Unzip the file

wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/repeat_consensus.txt.gz -P $FEATURES
#mv repeat_consensus.txt.gz $FEATURES/repeat_consensus.txt.gz
#gzip -d $FEATURES/repeat_consensus.txt.gz

wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/repeat_feature.txt.gz -P $FEATURES
#mv repeat_feature.txt.gz $FEATURES/repeat_feature.txt.gz
#gzip -d $FEATURES/repeat_feature.txt.gz

wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/seq_region.txt.gz -P $FEATURES
#mv seq_region.txt.gz $FEATURES/seq_region.txt.gz
#gzip -d $FEATURES/seq_region.txt.gz

wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/coord_system.txt.gz -P $FEATURES
#mv coord_system.txt.gz $FEATURES/coord_system.txt.gz

wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/seq_region_synonym.txt.gz -P $FEATURES


# #GET parameter
# _link="ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/"
# protein_feature_name="protein_feature.txt.gz"
# repeat_consensus_name="repeat_consensus.txt.gz"
# repeat_feature_name="repeat_feature.txt.gz"
# seq_region_name="seq_region.txt.gz"

# #WGET options
# proxy_option="-e use_proxy=yes -e http_proxy=www-cache:3128 -e https_proxy=www-cache:3128 -e ftp_proxy=www-cache:3128 "

# output_directory="$FEATURES"
# #GET cmd
# #echo ${_link}${protein_feature_name}
# #echo ${output_directory}${protein_feature_name}
# echo "wget ${_link}${protein_feature_name} -P ${FEATURES} ${proxy_option}"
# wget ${_link}${protein_feature_name} -P ${FEATURES} ${proxy_option}
# # && gzip -d ${output_directory}${protein_feature_name}

# #wget -q "${_link}${repeat_consensus_name}"  #&& gzip -d ${output_directory}${repeat_consensus_name}
# #wget -q "${_link}${repeat_feature_name}" #&& gzip -d ${output_directory}${repeat_feature_name}
# #wget -q "${_link}${seq_region_name}"  #&& gzip -d ${output_directory}${seq_region_name}
