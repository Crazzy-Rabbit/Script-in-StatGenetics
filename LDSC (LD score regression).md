#### 一、什么是LDSC？
>![image.png](https://upload-images.jianshu.io/upload_images/26461834-c93099c36f067fbf.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
- **LDSC**，LD分数回归，2015年由 `Brendan K Bulik-Sullivan` 提出的方法，旨在从样本量日益增加的`GWAS`结果的`inflation`中辨别 `confounding`（混杂因素）还是 `polygenicity`（多基因效应）。
#### 二、LDSC模型介绍
- LDSC认为，与`causal variant`处于LD的变异（即共享一定的遗传背景），该变异位点的测试统计量会因为其与`causal variant`的LD程度（通常以`r²`表示）而升高，且这种升高是成比例的。
- 与家系相关性（cryptic relatedness）或群体结构（population stratification）导致的测试统计量膨胀不同，这些因素不依赖于LD，而是由于共同遗传背景、遗传漂变等引起的统计量膨胀。这种膨胀不会与LD有相关性。
- 因此，LDSC通过SNP的LD分数构建了一个线性模型，来表征测试统计量的膨胀情况。同时，还能计算该`trait`的遗传力。
>$$
E(\chi^2 | l_j) = \frac{N \cdot h^2 \cdot l_j}M + N \cdot a + 1
$$
其中，$ l_j $为该SNP j 的LD score总和，左边为$ l_j $的$\chi$^2^，N为样本量，M为SNP数量，h^2^为该trait的遗传力，a为混杂因素（confounding）。因此，它本质上是个线性回归模型，该模型有两个未知数 h^2^ 和 a，通过拟合得到最适的 h^2^ 和 a。
![image.png](https://upload-images.jianshu.io/upload_images/26461834-21ca5a907901208b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
#### 三、LDSC分析实践
##### 1、数据格式转换
数据格式需要转换成它要求的`sumstat.gz`格式，使用hapmap3的SNP进行（LDSC提供了）
```
python munge_sumstats.py  --sumstats ${gwas} \
--merge-alleles ${SNPlist} \
--chunksize 500000 \
--a1 A1 \
--a2 A2 \
--out ${gwas}_ldsc

# --a1 effect allele a2 is another allele
```
##### 2、估计遗传力及判断confounding
对于连续性状，只需如下计算
```
python ldsc.py --h2 ${gwas}_ldsc.sumstats.gz \
        --ref-ld-chr ${REF_LD_CHR} \
        --w-ld-chr ${REF_LD_CHR} \
        --out ${gwas}_h2

# --ref-ld-chr  参考的LD score文件
```
对与二元性状，即疾病性状，需将其转换成`libility scale`
```
python ldsc.py --h2 ${gwas}_ldsc.sumstats.gz \
--ref-ld-chr ${REF_LD_CHR} \
--w-ld-chr ${REF_LD_CHR} \
--out ${gwas}_h2 \
--samp-prev 0.297 \
--pop-prev 0.1

# --pop-prev 为患病率
# --samp-prev 该summary data中的患病率
```


