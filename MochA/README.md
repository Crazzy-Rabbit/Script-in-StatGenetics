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
#### prepare data
```
vcf="..." # input VCF file with phased GT, LRR, and BAF
pfx="..." # output prefix
thr="..." # number of threads to use
crt="..." # tab delimited file with call rate information (first column sample ID, second column call rate)
sex="..." # tab delimited file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="..." # directory where output files will be generated
mkdir -p $dir
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
#### step 1: If your VCF does not include the `GC` field, this can be added with the command
```
bcftools +mochatools --no-version -Ob $vcf -- -t GC -f $ref
```
#### step 2: Creat a minimal binary VCF and `set to missing all genotypes that have low coverage or low genotyping quality`, as these can cause issues
```
bcftools view --no-version -h $vcf | \
  sed 's/^\(##FORMAT=<ID=AD,Number=\)\./\1R/' | \
  bcftools reheader -h /dev/stdin $vcf | \
  bcftools filter --no-version -Ou -e "FMT/DP<10 | FMT/GQ<20" --set-GT . | \
  bcftools annotate --no-version -Ou -x ID,QUAL,^INFO/GC,^FMT/GT,^FMT/AD | \
  bcftools norm --no-version -Ou -m -any --keep-sum AD | \
  bcftools norm --no-version -o $dir/$pfx.unpheased.bcf -Ob -f $ref --write-index
```
#### step3: Creat a list of `variants that will be excluded from modeling by both eagle and mocha` (Optional)
```
awk -F"\t" '$2<.97 {print $1}' $crt > sample_xcl_list.txt

echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin $dir/$pfx.unpheased.bcf | \
  bcftools view --no-version -Ou -S ^sample_xcl_list.txt | \
  bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
  bcftools view --no-version -Ou -G | \
  bcftools annotate --no-version -o $dir/$pfx.xcl.bcf -Ob --write-index \
    -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.2 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
    -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING

/bin/rm sample_xcl_list.txt
```
#### step4: phase genotypes
first you need split autosomes and chromosome X
```
bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 $dir/$pfx.unpheased.bcf $dir/$pfx.xcl.bcf | \
  bcftools view --no-version -Ou --min-ac 0 --exclude-uncalled | \
  bcftools annotate --no-version -Ou --remove ID,QUAL,INFO,^FMT/GT | \
  bcftools +scatter --no-version -Ob --output $dir --scatter $(echo chr{1..22},X | tr ' ' ',') --prefix $prefix.
```
If used `GRCh37` rather than `GRCh38`, use `--scatter $(echo {{1..22},X} | tr ' ' ',') --prefix $pfx.chr` instead

Phase VCF file `by chromosome` with `SHAPEIT5`
```
for chr in {1..22} X; do
  bcftools index --force $dir/$pfx.chr$chr.bcf
  zcat $map | sed 's/^23/X/' | awk -v chr=$chr '$1==chr {print $2,$3,$4}' > $dir/genetic_map.chr.$chr.txt
  phase_common \
    --thread $thr \
    --input $dir/$pfx.chr.$chr.bcf \
    --reference $panel_pfx${chr}$panel_sfx.bcf \
    --map $dir/genetic_map.chr$chr.txt \
    --region chr$chr \
    --output $dir/$pfx.chr$chr.pgt.bcf
done
```
If used `GRCh37` rather than `GRCh38`, use `--region $chr` instead. If phasing genotypes from WGS data, include the `--sequencing` option
