#### popcorn  trans-ethnic heritability and genetic correlation
```
popcorn 用于估计在不同群体中常见SNP的因果效应量的相关性，同一性状。该方法纳入了研究中所有SNP的关联信息（区别于其他需要事先进行LD-pruning的方法），
而且仅使用GWAS的概括数据，减少了计算负担，并避免了无法获取基因型数据的问题。
```
#### step 1 计算ld
```
popcorn compute -v 2 --bfile1 Popcorn/test/EUR_ref --bfile2 Popcorn/test/EAS_ref Popcorn/test/EUR_EAS_test_ge.cscore

# --gen_effect 计算遗传效应相关，两步都需加上该参数，无则结果为遗传冲击相关
# 软件提供了 eur-eas的
```
#### step 2 计算correlation
```
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --sfile1 $eur.txt.popcornin --sfile2 $eas.txt.popcornin correlation.txt

# summary data col: rsid/SNP a1/A1 a2/A2 af N beta SE 其中af对应a2的freq
# --gen_effect 计算遗传效应相关，两步都需加上该参数
# --sfile1 和 sfile2顺序需与 cscore 中一致
```

#### change file to popcornin
```
python $popcornIn.py infile 1349 outfile.txt
```
popcornIn.py
```
#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   popcornIn.py
@Time    :   2024/7/23 16:05:24
@Author  :   Lulu Shi    
@Mails   :   crazzy_rabbit@163.com
@link    :   https://github.com/Crazzy-Rabbit
'''

import os 
import sys
import gzip
import pandas as pd
from io import StringIO
"""
change BBJ SNP format same as UKB \n
such 1:612688_TCTC_T 1:rs1234145 \n
and change file header as rsid/SNP a1/A1 a2/A2 af N beta SE \n
"""

infile = sys.argv[1]
samplesize = sys.argv[2]
outfile = sys.argv[3]

if infile.endswith('.gz'):
    with gzip.open(infile, 'rt') as f_in:
            csv_content = f_in.read()   
    df = pd.read_csv(StringIO(csv_content), dtype={'BP':int}, sep="\\s+", engine='python')
    
    df.columns = ['SNP', 'CHR', 'BP', 'a2', 'a1', 'af', 'INFO', 'beta', 'SE', 'P']
    df['N'] = samplesize
    df.to_csv(f"{outfile}.popcornin", sep='\t', index=False)
else:
    df = pd.read_csv(infile, dtype={'POS':int}, sep='\\s+', engine='python')
    df.loc[df['SNP'].str.startswith('chr'), 'SNP'] = \
           df['SNP'].str.replace('chr', '', regex=True).str.replace("_", ":") + '_' + df['A1'] + '_' + df['A2']
    
    df.columns = ['SNP', 'CHR', 'BP', 'a2', 'a1', 'af', 'Rsq', 'beta', 'SE', 'P']
    df['N'] = samplesize
    df.to_csv(f"{outfile}.popcornin", sep='\t', index=False)
```
