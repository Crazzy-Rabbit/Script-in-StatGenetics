#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   add_samplesize.py
@Time    :   2024/07/19 17:15:08
@Author  :   Lulu Shi 
@Mails   :   crazzy_rabbit@163.com
@link    :   https://github.com/Crazzy-Rabbit
'''

import os
import gzip
import click
import pandas as pd
from io import StringIO


@click.command()
@click.option('--infile', '-I',help="input file you need to add samplesize", required=True)
@click.option('--samplesize', "-S", help="sample size you want to add", type=int, required=True)
@click.option('--output', '-O', help="output txt.gz you added samplesize", required=True)
def main(infile, samplesize, output):
    """
    add sample size and change GWAS summary data header to MR-MEGA in file format \n
    such: MARKERNAME EA NEA BETA SE EAF N CHR POS
    """
    if infile.endswith(".gz"):
        with gzip.open(infile, 'rt') as f_in:
            csv_content = f_in.read()   
        df = pd.read_csv(StringIO(csv_content), dtype={'POS':int}, 
                         sep='\\s+', engine='python')  # 使用StringIO处理非字符串内容
        df['nSample'] = samplesize
        df[['SNP', 'ALLELE0', 'ALLELE1', 'BETA', 'SE',
            'A1FREQ', 'nSample', 'CHR', 'BP']].to_csv(f"{f_out}.MRMEGAin", sep='\t', index=False)

        os.system(f"gzip {f_out}.MRMEGAin")
        
    else:
        df = pd.read_csv(infile, dtype={'POS':int}, 
                         sep='\\s+', engine='python')
        df['nSample'] = samplesize
        df[['SNP', 'A1', 'A2', 'BETA', 'SE',
            'A1Frq', 'nSample', 'CHR', 'POS']].to_csv(f"{output}.MRMEGAin", sep='\t', index=False)

if __name__ == '__main__':
    main()
