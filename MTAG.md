#### `mtag` (Multi-Trait Analysis of GWAS)

使用GWAS概括新数据进行多表型联合分析的方法，相比于单表型的GWAS，
`MTAG`可以利用关联表型的信息提升目标表型的检验统计power

核心假设: 对于不同表型，所有的SNP都**共享同一个效应量的方差协方差矩阵**
#### 需要python2.7环境，及几个特定版本的python包
```
数据格式： snpid    chr    bpos    a1    a2    freq   z   pval    n
可用 BETA和SE代替z， --use_beta_se 这个参数有问题，因此作者建议使用 Z，Z=BETA/SE
```
数据可以指定表头名，且多个数据以逗号`,`分开
```
conda activate mtag

python mtag.py --sumstats trait1.txt,trait2.txt --chr_name CHR --snp_name SNP --bpos_name BP \
           --a1_name A1 --a2_name A2 --eaf_name A1Freq --p_name P --stream_stdout \
           --n_name N --z_name Z --n_min 0.0 --out mtag.result
```
