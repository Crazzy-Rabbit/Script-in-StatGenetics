#=========================================================
# make_annot.py to generate annot.gz file 
# using 100K widow as Hilary K. Finucane et al. Nat Genet
conda activate ldsc

conf_dir="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Eshel/LDSC"

#- make thin annot
bfile="/public/home/shilulu/Wulab/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
coord="/public/home/shilulu/software/ldsc/Gene_coord.txt"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/home/shilulu/Wulab/LDSC"

files=( *.bed )
total=${#files[@]}
batch_size=1
i=0

while [ $i -lt $total ]; do 
    batch=()
    job_ids=()

    for ((j=0; j<batch_size && (i + j) < total; j++)); do 
        file=("${files[$((i + j))]}")
        annot=$(basename -- "${file}" ".bed")
        
        if [ -f "${annot}.1.l2.ldscore.gz" ]; then
            echo "✅ Found ${annot}.1.l2.ldscore.gz, skipping $annot..."
            continue
        fi

        cmd1="python ${ldsc}/make_annot.py \
        --gene-set-file ${annot}.bed \
        --gene-coord-file $coord \
        --bimfile ${bfile}.{TASK_ID}.bim \
        --windowsize 100000 \
        --annot-file ${annot}.{TASK_ID}.annot.gz"

        sub1=$(qsubshcom "$cmd1" 1 10G ldsc_anot 1:00:00 "-array=1-22")

        cmd2="python ${ldsc}/ldsc.py --l2 \
        --bfile ${bfile}.{TASK_ID} \
        --print-snps ${ldsc_dir}/listHM3.txt \
        --ld-wind-cm 1 \
        --annot ${annot}.{TASK_ID}.annot.gz \
        --thin-annot \
        --out ${annot}.{TASK_ID}"
        jobid=$(qsubshcom "$cmd2" 1 10G ldsc_comp 1:00:00 "-array=1-22 -wait=$sub1")
        echo "🚀 Submitted job $jobid for $annot, waiting for it finish..."
        
        job_ids+=("$jobid")
    done 

    echo "⏳ Waiting for jobs: ${job_ids[*]}"
    times=0
    while true; do 
        sleep 30
        all_done=true 
        for jid in "${job_ids[@]}"; do 
            if squeue -j "$jid" | grep -q "$jid"; then 
                all_done=false 
                times=$((times + 30))
                echo "⏳ Jobs in batch still running ${times}s..."
                break
            fi
        done 

        if $all_done; then
            echo "✅ All jobs in batch completed."
            break
        else
            echo "⏳ Still waiting for some jobs in batch to finish..."
        fi 
    done 

    i=$((i + batch_size))
done

#- generate .ldcts file for my sc data
dir=`pwd`
ls *bed | while read id; do
  tissue=$(basename -- ${id} | sed 's/\.bed$//')
  if [[ "$tissue" == "control" ]]; then
    echo "control not to be the tissue or cell type"
  else
    cell=${dir}/$(basename -- ${id} "bed")
    control=${dir}/"control."
    echo "${tissue} ${cell},${control}" >> scRNA_gene_expr.ldcts
  fi
done