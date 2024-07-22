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
    if infile.endswith(".gz"):
        with gzip.open(infile, 'rt') as f_in:
            csv_content = f_in.read()
   
        df = pd.read_csv(StringIO(csv_content))  # 使用StringIO处理非字符串内容
        df['nSample'] = samplesize
        with gzip.open(output, 'wt') as f_out:
            df.to_csv(f_out, sep='\t', index=False)
    else:
        with open(infile, 'r') as f_in:
            csv_content = f_in.read()

        df = pd.read_csv(StringIO(csv_content))
        df['nSample'] = samplesize
        with open(output, 'w') as f_out:
            df.to_csv(f_out, sep='\t', index=False)

if __name__ == '__main__':
    main()
