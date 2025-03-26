
### 1.`METAL` 进行meta GWAS
只需要一个配置文件即可，还是比较简单的
他是基于`fixed-effect model`

$$
\hat\beta_{META} = \sum_{c}\frac{\hat\beta_c}{{SE_{\hat\beta_c}^2}}
$$

$$
 SE_{\hat\beta_{META}} = \sqrt{\frac{1}{\sum_{c}\frac{1}{SE_{\hat\beta_c}^2}}}
$$
```
metal METAL.conf
```

### 2.`MR-MEGA`进行跨族裔GWAS荟萃分析
```
https://gwaslab.org/category/gwas-%e5%85%a8%e5%9f%ba%e5%9b%a0%e7%bb%84%e5%85%b3%e8%81%94%e5%88%86%e6%9e%90/13-%e8%8d%9f%e8%90%83%e5%88%86%e6%9e%90-meta-analysis/
```
通过**跨族裔荟萃回归来检测并精准定位复杂表型关联信号**的工具

### 3.`MTAG`进行多个相似表型的meta分析
使用GWAS概括新数据进行多表型联合分析的方法，相比于单表型的GWAS，`MTAG`可以利用关联表型的信息提升目标表型的检验统计power

核心假设: 对于不同表型，所有的SNP都**共享同一个效应量的方差协方差矩阵**
```
需要python2.7环境，及几个特定版本的python包
数据格式： snpid    chr    bpos    a1    a2    freq   z   pval    n
可用 BETA和SE代替z， --use_beta_se 这个参数有问题，因此作者建议使用 Z，Z=BETA/SE
```
数据可以指定表头名，且多个数据以逗号`,`分开
```
conda activate mtag

python mtag.py --sumstats trait1.txt,trait2.txt \
               --stream_stdout \
               --n_min 0.0 \
               --out mtag.result
```
输入文件很多时的简单方法（没办法，太懒了，不想一个一个加）
```
input_str=`ls *.mtag | tr '\n' ',' | sed 's/,$//'`

python mtag.py --sumstats $input_str  \
           --stream_stdout \
           --n_min 0.0 \
           --out mtag.result
```