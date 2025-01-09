#! usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :    ReformatGWAScols.py
@Time    :    2025/01/09 17:41:09
@Author  :    Lulu Shi
@Mails   :    crazzy_rabbit@163.com
@line    :    https://github.com/Crazzy-Rabbit
'''
import click
import pandas as pd

@click.command()
@click.option('--gwas', type=str, help="GWAS summary data", required=True)
@click.option('--SNP', type=str, help="SNP colnames", required=True)
@click.option('--A1', type=str, help="A1 colnames for effect allele", required=True)
@click.option('--A2', type=str, help="A2 colnames for another allele", required=True)
@click.option('--freq', type=str, help="frequency colnames for effect allele", required=True)
@click.option('--beta', type=str, help="beta colnames for effect allele", required=True)
@click.option('--se', type=str, help="se colnames", required=True)
@click.option('--pval', type=str, help="p colnames", required=True)
@click.option('--num', type=str, help="sample size colnames or sample size number", required=True)
@click.option('--outprefix', type=str, help="outprefix for your file", required=True)
def process_gwas(gwas, SNP, A1, A2, freq, beta, se, pval, num, outprefix):
    '''
    reformat gwas colnames and select cols that used in downstream analysis
    '''
    # Load GWAS summary data
    gwas_dt = pd.read_csv(gwas, delimiter='\t')

    # Check if num is a numeric value or column name
    if num.isdigit():
        gwas_dt['N'] = int(num)
        columns_map = {SNP: 'SNP', A1: 'A1', A2: 'A2', freq: 'freq', beta: 'beta', se: 'se', pval: 'p', 'N': 'N'}
    elif num in gwas_dt.columns:
        columns_map = {SNP: 'SNP', A1: 'A1', A2: 'A2', freq: 'freq', beta: 'beta', se: 'se', pval: 'p', num: 'N'}
    else:
        raise ValueError(f"Column {num} neither numeric nor found in the data")

    # Rename columns
    gwas_dt.rename(columns=columns_map, inplace=True)
    gwas_dt.to_csv(outprefix, sep='\t', index=False)

if __name__ == '__main__':
    process_gwas()
