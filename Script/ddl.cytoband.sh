#!/usr/bin/env bash

# exit when any command fails
set -e
set -o nounset

#PATHS
readonly CWD=$(pwd)
readonly BASE_DIR=$(dirname "$CWD")
readonly CYTOBAND=${BASE_DIR}/Cytoband

#Check if folder exists else create it
if [ ! -d $CYTOBAND ]; then mkdir -p $CYTOBAND ;fi

#Download - Move to folder - Unzip the file
wget -nc ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg38/database/cytoBand.txt.gz
mv cytoBand.txt.gz $CYTOBAND/cytoBand.txt.gz
gzip -d $CYTOBAND/cytoBand.txt.gz



# GET parameter
# _link="ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg38/database/"
# cytoBand_name="cytoBand.txt.gz"

# #GET cmd
# #echo ${_link}${protein_feature_name}
# #echo ${output_directory}${protein_feature_name}
# #echo "wget ${_link}${cytoBand_name} -P ${CYTOBAND} " 
# wget -Nq ${_link}${cytoBand_name} -P ${CYTOBAND} -e use_proxy=yes -e http_proxy=www-cache:3128 -e https_proxy=www-cache:3128 -e ftp_proxy=www-cache:3128 && gzip -d ${CYTOBAND}${cytoBand_name}

# #wget -q "${_link}${repeat_consensus_name}"  #&& gzip -d ${output_directory}${repeat_consensus_name}
# #wget -q "${_link}${repeat_feature_name}" #&& gzip -d ${output_directory}${repeat_feature_name}
# #wget -q "${_link}${seq_region_name}"  #&& gzip -d ${output_directory}${seq_region_name}


# wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/repeat_consensus.txt.gz

# wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/repeat_feature.txt.gz

# wget -nc ftp://ftp.ensembl.org/pub/release-90/mysql/homo_sapiens_core_90_38/seq_region.txt.gz