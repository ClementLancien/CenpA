#!/usr/bin/env python
# coding: utf-8

# Columns definitions :
# https://genome.ucsc.edu/cgi-bin/hgTables
# From :
#

import pandas as pd

cytoband = pd.read_table( 	'../Cytoband/cytoBand.txt',
							header = None,
							sep = '\t',
							dtype = 'str',
							names = [
										'chromosome_name',
										'chromosome_start',
										'chromosome_end',
										'name',
										'giemsa_stain'
									]
							
						)
#replace chr# with @
cytoband['chromosome_name'] = cytoband['chromosome_name'].str.replace('[chr]', '')

#keep centromeric and non centromeric
centromeric = cytoband[cytoband['giemsa_stain'] == 'acen']
noncentromeric = cytoband[cytoband['giemsa_stain'] != 'acen']

#export to bed files
#keep 'chromosome_name','chromosome_start','chromosome_end' columns for bed files
centromeric.to_csv('../Cytoband/centromeric.bed', columns =['chromosome_name','chromosome_start','chromosome_end'], sep='\t', index = False, header=None)
noncentromeric.to_csv('../Cytoband/noncentromeric.bed', columns =['chromosome_name','chromosome_start','chromosome_end'], sep='\t', index = False, header=None)
