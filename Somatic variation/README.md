### The pipeline of somatic variation calling using `Mutect2 of GATK4.1.2`
step1: call somatic variation using `Mutect2`
> Note: The germline resource and ref panel of normal variation can download from [here](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/)

```
GATK4="/public/software/apps/GATK/4.1.2.0/gatk"
REF="/public/home/shilulu/reference/hs37d5/hs37d5.fa"
tmpdir="/public/share/wchirdzhq2022/shilulu/tmp"
PON="/public/share/wchirdzhq2022/Wulab_share/somatic/somatic-b37_Mutect2-WGS-panel-b37.vcf"
Germline="/public/share/wchirdzhq2022/Wulab_share/somatic/somatic-b37_af-only-gnomad.raw.sites.vcf"

indir="/public/share/wchirdzhq2022/Wulab_share/somatic"
outdir="/public/share/wchirdzhq2022/Wulab_share/somatic"

${GATK4}  Mutect2 \
      -R ${REF} \
      -I ${indir}/HG002_bwa_sort_rg.bam \
      --panel-of-normals  ${PON}\
      --germline-resource ${Germline} \
      --tmp-dir ${tmpdir} \
      -O HG002_bwa_sort_rg.vcf.gz
```
