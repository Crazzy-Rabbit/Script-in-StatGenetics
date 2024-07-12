#### Install
```
### install bcftools >= 1.20
wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar xjvf bcftools-1.20.tar.bz2

cd bcftools-1.20/
/bin/rm -f plugins/{{mocha,beta_binom,genome_rules}.h,{mocha,trio-phase,mochatools,extendFMT}.c}
wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/{{mocha,beta_binom,genome_rules}.h,{mocha,trio-phase,mochatools,extendFMT}.c}
make

### add below info in the last line of .bashrc 
vim ~/.bashrc
export BCFTOOLS_PLUGINS="$HOME/software/bcftools-1.20/plugins:$BCFTOOLS_PLUGINS"
```
#### `MochA` input vcf of WGS
```
##fileformat=VCFv4.2
##INFO=<ID=GC,Number=1,Type=Float,Description="GC ratio content around the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	752566	rs3094315	G	A	.	.	GC=0.3675	GT:AD	1|1:0,31
1	776546	rs12124819	A	G	.	.	GC=0.435	GT:AD	0|1:21,23
1	798959	rs11240777	G	A	.	.	GC=0.4075	GT:AD	0|0:31,0
1	932457	rs1891910	G	A	.	.	GC=0.6425	GT:AD	1|0:18,14
```
#### The vcf of WGS need columns `GC GT AD`, where `GT AD` were contained in the `GATK` out vcffile
#### If your VCF does not include the `GC` field, this can be added with the command
```
bcftools +mochatools --no-version -Ob $vcf -- -t GC -f $ref
```

