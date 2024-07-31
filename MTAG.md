#### `mtag` (Multi-Trait Analysis of GWAS)

使用GWAS概括新数据进行多表型联合分析的方法，相比于单表型的GWAS，
MTAG可以利用关联表型的信息提升目标表型的检验统计power,
核心假设: 对于不同表型，所有的SNP都共享同一个效应量的方差协方差矩阵
```
数据格式： snpid    chr    bpos    a1    a2    freq   z   pval    n
可用 BETA和SE代替z， --use_beta_se 这个参数有问题，因此作者建议使用 Z，Z=BETA/SE
```
