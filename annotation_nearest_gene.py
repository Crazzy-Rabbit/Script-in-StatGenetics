# hg37 to hg38
# library(dplyr)
# library(data.table)
# df = fread("ARHL_MVP_AJHG_BBJ_reformatMETAL_addchr.gz")
# df_out = df %>% filter(p < 5e-8)
# fwrite(df_out, file="5e-8_ARHL.txt", sep="\t")

df_out = fread("5e-8_ARHL.txt")
merge_out = data.table()
for (chr in c(1:22)){
    setwd("/public/share/wchirdzhq2022/Wulab_share/dbSNP/GRCh38p7")
    Chr = fread(paste0("chr", chr, ".txt"))[, c(1:3)]
    colnames(Chr) = c("CHR", "POS", "SNP")
    df_merge = merge(df_out, Chr, by = "SNP")[, c("CHR", "POS", "SNP")]
    merge_out = rbind(merge_out, df_merge)
}
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta")
fwrite(merge_out, file="5e-5_ARHL_hg38.txt", sep="\t")


# annotation gene
import os 
import pandas as pd
os.chdir("F:\\百度网盘同步\\BaiduSyncdisk\\24-28川大华西\\2_project_hearing loss\\new_run\\gene_annotation")
os.getcwd()

gwa_df = pd.read_csv('5e-8_ARHL.txt', sep='\t')
glist_df = pd.read_csv('glist-hg19', sep=' ', header=None, names=['CHR', 'START', 'END', 'GENE'])

glist_df['CHR'] = glist_df['CHR'].astype(str)
gwa_df['CHR'] = gwa_df['CHR'].astype(str)

gwa_df['GENE'] = ''

for index, row in gwa_df.iterrows():
    chr = row['CHR']
    pos = row['POS']
    genes = glist_df[glist_df['CHR'] == chr]
    if genes.empty:
        gwa_df.at[index, 'GENE'] = 'No genes on chromosome'
        continue
    match = genes[(genes['START'] <= pos) & (genes['END'] >= pos)]
    if not match.empty:
        gwa_df.at[index, 'GENE'] = match['GENE'].values[0]
    else:
        left_genes = genes[genes['END'] < pos].sort_values(by='END', ascending=False)
        right_genes = genes[genes['START'] > pos].sort_values(by='START', ascending=True)  
        if not left_genes.empty and not right_genes.empty:
            nearest_genes = f"{left_genes['GENE'].values[0]}/{right_genes['GENE'].values[0]}"
        elif not left_genes.empty:
            nearest_genes = left_genes['GENE'].values[0]
        elif not right_genes.empty:
            nearest_genes = right_genes['GENE'].values[0]
        else:
            nearest_genes = 'No neighboring genes found'
        
        gwa_df.at[index, 'GENE'] = nearest_genes


gwa_df.to_csv('5e-8_ARHL_with_genes.txt', sep='\t', index=False)
