#-------------------------------------------------//
# make annot  for control and per cell type 
#-------------------------------------------------//
conda activate ldsc
#-------make annot----//
#! /bin/bash
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"

files=( *.GeneSet )
total=${#files[@]}
batch_size=3
i=0

while [ $i -lt $total ]; do 
    batch=()
    job_ids=()

    for ((j=0; j<batch_size && (i + j) < total; j++)); do 
        file=("${files[$((i + j))]}")
        annot=$(basename -- "${file}" ".GeneSet")

        cmd="python ${ldsc}/make_annot.py \
        --gene-set-file ${annot}.GeneSet \
        --gene-coord-file ${ldsc}/Gene_coord.txt \
        --bimfile ${bfile}.{TASK_ID}.bim \
        --windowsize 100000 \
        --annot-file ${annot}.{TASK_ID}.annot.gz"

        jobid=$(qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22")
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

#-------compute LD----//
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/share/wchirdzhq2022/Wulab_share/LDSC"

files=( *.GeneSet )
total=${#files[@]}
batch_size=3
i=0

while [ $i -lt $total ]; do 
    batch=()
    job_ids=()

    for ((j=0; j<batch_size && (i + j) < total; j++)); do 
        file=("${files[$((i + j))]}")
        annot=$(basename -- "${file}" ".GeneSet")

        cmd="python ${ldsc}/ldsc.py --l2 \
        --bfile ${bfile}.{TASK_ID} \
        --print-snps ${ldsc_dir}/listHM3.txt \
        --ld-wind-cm 1 \
        --annot ${annot}.{TASK_ID}.annot.gz \
        --thin-annot \
        --out ${annot}.{TASK_ID}"
    
        jobid=$(qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22")
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

#-------generate .ldcts file for my sc data----//
dir=`pwd`
ls *GeneSet | while read id; do
    tissue=$(basename -- ${id} | sed 's/\.GeneSet$//')
    if [[ "$tissue" == "control" ]]; then
        echo "control not to be the tissue or cell type"
    else
        cell=${dir}/$(basename -- ${id} "GeneSet")
        control=${dir}/"control."
        echo "${tissue} ${cell},${control}" >> Bone_Teeth_scRNA.ldcts
    fi
done
# rm control file in cell 