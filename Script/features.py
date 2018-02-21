#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re

regions = pd.read_table(    '../Features/seq_region.txt.gz', 
                            header = None, 
                            names = [   'seq_region_id',
                                        'name',
                                        'coord_system_id',
                                        'length'
                                    ]
                        )
regions = regions[
				(regions['coord_system_id'] == 3) | #check coord_system file
				(regions['coord_system_id'] == 4)
			    ]
regions = regions[['seq_region_id', 'name']]
region = regions[['seq_region_id']]
name = regions[['name']]

#regions = regions['seq_region_id'].tolist()
consensus = pd.read_table(  '../Features/repeat_consensus.txt.gz',
                            header = None, 
                            names = [   'repeat_consensus_id',
                                        'repeat_name',
                                        'repeat_class',
                                        'repeat_type',
                                        'repeat_consensus'
                                    ],
                            usecols = ['repeat_consensus_id','repeat_name']
                         )

consensus = consensus[
						consensus['repeat_name'].str.contains('^[A][L][R][/][A][L][P][H][A]', flags=re.IGNORECASE)
					 ]

repeat_consensus_id = consensus['repeat_consensus_id'].values[0]

#print repeat_consensus_id
repeats = pd.read_table(    '../Features/repeat_feature.txt.gz',
                            header = None,
                            names = [   'repeat_feature_id',
                                        'seq_region_id',
                                        'seq_region_start',
                                        'seq_region_end',
                                        'seq_region_strand',
                                        'repeat_start',
                                        'repeat_end',
                                        'repeat_consensus_id',
                                        'analysis_id', 
                                        'score'
                                    ],
                            usecols = ['seq_region_id', 'seq_region_start', 'seq_region_end', 'repeat_consensus_id' ]
                        )
repeats = repeats[
					repeats['repeat_consensus_id'] == repeat_consensus_id
				 ]

result = repeats.merge(regions, on='seq_region_id')
result = result.sort_values(by=['name'])
result.to_csv('../Features/repeats.csv', header=False, index=False, columns=['name', 'seq_region_start', 'seq_region_end'])




