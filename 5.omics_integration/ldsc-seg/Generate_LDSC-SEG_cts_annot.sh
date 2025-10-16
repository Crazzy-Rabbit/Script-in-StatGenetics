# split chr for ref genotype file
cd "/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eas"

for i in $(seq 1 22); do
    plink --bfile g1000_eas --chr $i --make-bed --out g1000_eas.${i}
done

# run make_annot.py 为control 和 每个cell type生成
conda activate ldsc
bfile="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eas/g1000_eas"
ldsc="/public/home/shilulu/software/ldsc"
wkdir="/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/EAS"

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  # for chr in {1..22}; do 
  cmd="python ${ldsc}/make_annot.py \
      --gene-set-file ${annot}.GeneSet \
      --gene-coord-file ${ldsc}/ENSG_coord.txt \
      --bimfile ${bfile}.{TASK_ID}.bim \
      --windowsize 100000 \
      --annot-file ${wkdir}/${annot}.{TASK_ID}.annot.gz"
  # done 
  qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22"
done

## Step 2: Computing LD scores with an annot file
# if reference bfile not have cM, you can use --ld-wind-kb to replace
bfile="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eas/g1000_eas"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/share/wchirdzhq2022/Wulab_share/LDSC"
wkdir="/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/EAS"

awk '{if ($1!="SNP") {print $1} }' ${ldsc_dir}/w_hm3.snplist > ${ldsc_dir}/listHM3.txt

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  cmd="python ${ldsc}/ldsc.py --l2 \
    --bfile ${bfile}.{TASK_ID} \
    --print-snps ${ldsc_dir}/listHM3.txt \
    --ld-wind-kb 1000 \
    --annot ${wkdir}/${annot}.{TASK_ID}.annot.gz \
    --thin-annot \
    --out ${wkdir}/${annot}.{TASK_ID}"
  qsubshcom "$cmd" 1 10G ldsc_l2 20:00:00 "-array=1-22"
done

# generate .ldcts file for my sc data
wkdir="/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/EAS"
ls *GeneSet | while read id; do
  tissue=$(basename -- "${id}" | sed 's/^Sc_P8_12_20_//;s/.GeneSet$//')
  cell=${wkdir}/$(basename -- "${id}" "GeneSet")
  control=${wkdir}/"Sc_P8_12_20_control."
  if [ "$tissue" == "control" ]; then
    continue 
  else
    echo "${tissue} ${cell},${control}" >> ${wkdir}/Mouse_cochleae_gene_expr.ldcts
  fi
done

