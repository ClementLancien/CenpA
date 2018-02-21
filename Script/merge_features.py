#!/usr/bin/env python
# coding: utf-8

import pandas as pd

regions = pd.read_table(    '../Features/seq_region.txt.gz', 
                            header = None, 
                            names = [   'seq_region_id',
                                        'name',
                                        'coord_system_id',
                                        'length'
                                    ]
                        )


consensus = pd.read_table(  '../Features/repeat_consensus.txt.gz',
                            header = None, 
                            names = [   'repeat_consensus_id',
                                        'repeat_name',
                                        'repeat_class',
                                        'repeat_type',
                                        'repeat_consensus'
                                    ]
                         )

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
                                    ]
                        )

repeats = pd.merge(repeats, regions, on = ['seq_region_id', 'seq_region_id'])

repeats = pd.merge(repeats, consensus, on = ['repeat_consensus_id', 'repeat_consensus_id'])
 
#repeats = pd.merge([(repeats['repeat_class'] != 'dust') & (repeats['repeat_class'] != 'trf')] 


#repeats = repeats[['name', 'seq_region_start', 'seq_region_end', 'repeat_name', 'repeat_class', 'repeat_type']]

#repeats.columns = ['chromosome', 'start', 'end', 'repeat_name', 'repeat_class', 'repeat_type']

repeats.to_csv('../Features/repeats.csv', index = False)