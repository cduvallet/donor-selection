#!/usr/bin/env python
"""
Make the clean Kump 2018 metadata
"""

import pandas as pd

# This file is fine as-is.

fmeta = 'data/raw/kump2018/mapping_file.txt'
fout = 'data/clean/kump2018.metadata.txt'
df = pd.read_csv(fmeta, sep='\t', index_col=0)
df.to_csv(fout, sep='\t')
