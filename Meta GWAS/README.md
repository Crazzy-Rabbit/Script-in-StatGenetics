
#### `METAL` 进行meta GWAS
```
python RunMetal.py --colsnp SNP --eallele ALLELE1 --nallele ALLELE0 --eallfrq A1FREQ --beta BETA --se SE --pvalue P_BOLT_LMM --traitlist UKB_milk_meta.txt --metalfile UKB_milk_metal --outfile ~/Meta_GWAS/UKB-dietary_meta/UKB_milk_meta
```

`MR-MEGA`进行跨族裔GWAS荟萃分析
```
https://gwaslab.org/category/gwas-%e5%85%a8%e5%9f%ba%e5%9b%a0%e7%bb%84%e5%85%b3%e8%81%94%e5%88%86%e6%9e%90/13-%e8%8d%9f%e8%90%83%e5%88%86%e6%9e%90-meta-analysis/
```
通过**跨族裔荟萃回归来检测并精准定位复杂表型关联信号**的工具
