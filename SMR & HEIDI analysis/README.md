### `smr` & `heidi`分析
#### run SMR and HEIDI
##### 1
```
smr --bfile --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mysmr --thread-num 10
```
`--bfile` SNP基因型数据
`--gwas-summary` summary-level data from GWAS. 格式如下
```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830

每一列含义如下：
SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size.
```
##### 2
···
smr --bld mybld --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mysmr --thread-num 10 
···
