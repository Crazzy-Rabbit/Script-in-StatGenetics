#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   RunMetal.py
@Time    :   2024/07/19 13:16:56
@Author  :   Lulu Shi 
@Mails   :   crazzy_rabbit@163.com
@link    :   https://github.com/Crazzy-Rabbit
'''

import os
import click

# software dir
Metal = "/public/home/shilulu/software/metal/metal"

def GenerateMetalInFinal(colsnp, eallele, nallele, eallfrq, beta, se, pvalue,  traitlist, metalfile, outfile):
# MARKER  指定SNP ID对应的列标题
# ALLELE  指定test allele和other allele对应的列标题
# PVALUE  指定P值对应的列标题
# EFFECT  指定效应方向对应的列标题
# WEIGHT  指定样本量对应的列标题
# PROCESS 指定gwas结果对应的文件
# OUTFILE 指定输出结果的名称（注意：此为空格分隔的前缀和后缀）
# ANALYZE 表示进行分析
    Final = []
    # 将要meta的每个trait的metal处理信息添加到Final中
    Final.append("SCHEME\tSTDERR \n")
    with open (traitlist, 'r') as f:
        for Trait in f:
            Final.append("# Describe and process the input files \n")
            Final.append(f"MARKER\t{colsnp} \n")
            Final.append(f"ALLELE\t{eallele} {nallele} \n")
            Final.append(f"FREQ\t{eallfrq} \n")
            Final.append(f"EFFECT\t{beta} \n")
            Final.append(f"STDERR\t{se} \n")
            #Final.append(f"WEIGHT\t{nSample} \n")
            Final.append(f"PVAL\t{pvalue} \n")
            Final.append(f"PROCESS\t{Trait} \n\n")
    Final.append(f"OUTFILE {outfile} .tbl \n")
    Final.append("ANALYZE \n")
    Final.append("ANALYZE HETEROGENEITY \n")
    Final.append("QUIT")
    
    with open (metalfile, 'w') as Out:
        for i in Final:
            Out.write(i)

@click.command()
@click.option('--colsnp', help="the column name of SNP", type=str, required=True)
@click.option('--eallele', help="the column name of effect allele", type=str, required=True)
@click.option('--nallele', help="the column name of noneffect allele", type=str, required=True)
@click.option('--eallfrq', help="the column name of frequency of effect allele", type=str, required=True)
@click.option('--beta', help="the column name of effect size", type=str, required=True)
@click.option('--se', help="the column name of standard error", type=str, required=True)
@click.option('--pvalue', help="the column name of P value", type=str, required=True)
#@click.option('--nSample', help="the column name of sample size", type=str, required=True)
@click.option('--traitlist', help="the list file of meta trait, one row one filename", type=str, required=True)
@click.option('--metalfile', help="the filename of metal infile which you want to define, any str is ok, default is metalfile", type=str, default="metalfile")
@click.option('--outfile', help="the outfile name of metal", type=str, required=True)
def main(colsnp, eallele, nallele, eallfrq, beta, se, pvalue,  traitlist, metalfile, outfile):

    GenerateMetalInFinal(colsnp, eallele, nallele, eallfrq, beta, se, pvalue, traitlist, metalfile, outfile)
    os.system(f"srun -N 1 \
                     -J Metal_{metalfile} \
                     -p normal \
                    {Metal} {metalfile}")

if __name__ == '__main__':
    main()
