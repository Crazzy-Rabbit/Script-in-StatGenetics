#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   run_SMRinPy_filter.py
@Time    :   2024/10/23 20:21:31
@Author  :   Lulu Shi 
@Mails   :   crazzy_rabbit@163.com
@link    :   https://github.com/Crazzy-Rabbit
'''
import os
import argparse
import pandas as pd

def run_smr_qtl(gwas, qtl, outname):
    SMR = "/public/home/lijunpeng/smr-1.3.1-linux-x86_64/smr-1.3.1"
    bfile = "/public/home/lijunpeng/smr-1.3.1-linux-x86_64/g1000_eur/g1000_eur"

    os.system(f'{SMR} --bfile {bfile} \
                      --beqtl-summary {qtl} \
                      --gwas-summary {gwas} \
                      --diff-freq-prop 0.5 \
                      --maf 0.01 \
                      --cis-wind 2000 \
                      --out {outname} \
                      --thread-num 10 > {outname}.log')
def run_filter(smr, outname):
    smr_file = pa.read_csv(smr, sep="\t")
    smr_tsh = "{:.2e}".format(0.05 / len(smr_file))

    flt_file = smr_file[(smr_file['p_SMR'] <= 0.05 / len(smr_file)) & (smr_file['p_HEIDI'] >= 0.01)]
    flt_file.to_csv(f"{outname}_{smr_tsh}smr_0.01heidi.txt", sep="\t", index=False, quoting=False)

def main():
    parser = argparse.ArgumentParser(description="Run SMR and filter results.")
    parser.add_argument("--func", type=str, default="both", help="Function to run: 'run_smr_qtl', 'run_filter', or 'both'")
    parser.add_argument("--gwas", type=str, help="GWAS summary data to run smr")
    parser.add_argument("--qtl", type=str, help="QTL file to run smr")  
    parser.add_argument("--smr", type=str, help="SMR file for filtering")
    parser.add_argument("--out", type=str, help="Output prefix")
    
    args = parser.parse_args()
    
    gwas = args.gwas
    qtl = args.qtl
    out = args.out
    func = args.func
    smr = args.smr

    if func == "run_smr_qtl":
        if not gwas or not smr:
            raise ValueError("When func is 'run_smr_qtl', --gwas and --qtl must be provided.")
        run_smr_qtl(gwas, qtl, out)
    elif func == "run_filter":
        if not smr:
            raise ValueError("When func is 'run_filter', --smr must be provided.")
        run_filter(smr, out)
    else:
        if not gwas or not smr:
            raise ValueError("When func is 'both' or not provided, --gwas and --qtl must be provided.")
        run_smr_qtl(gwas, qtl, out)
        run_filter(f"{out}.smr", out)

if __name__ == '__main__':
    main()
