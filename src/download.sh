#!/bin/bash

#Download homo sapiens genome

wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz -P ../Data/Reference_Genome/

#Download homo sapiens features

wget ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/protein_feature.txt.gz -P ../Data/Features/
wget ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/repeat_consensus.txt.gz -P ../Data/Features/
wget ftp://ftp.ensembl.org/pub/release-91/mysql/homo_sapiens_otherfeatures_91_38/repeat_feature.txt.gz - P ../Data/Features/
