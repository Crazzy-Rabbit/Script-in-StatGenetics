#---- liftover hg38 to hg19
library(rtracklayer)
library(GenomicRanges)

regions = unique(as.character(merged_dt$Gene))
regions_str = strsplit(regions, "-")
bed_df = do.call(rbind, lapply(seq_along(regions), function(i) {
    x = regions_str[[i]]
    chr = x[1]
    start = as.numeric(x[2]) - 1
    end = as.numeric(x[3])
    peak_id = regions[i]
    return(data.frame(chr=chr, start=start, end=end, peak_id=peak_id))
}))

gr <- GRanges(seqnames = bed_df$chr, 
            ranges = IRanges(start=bed_df$start, end=bed_df$end),
            peak_id = bed_df$peak_id)
chain <- import.chain("/public/home/shilulu/Wulab/refGenome/hg38/hg38ToHg19.over.chain")

gr_new_list <- liftOver(gr, chain)
gr_new <- unlist(gr_new_list)
gr_clean <- gr_new[elementNROWS(gr_new_list) == 1]
df_clean <- as.data.frame(gr_clean)[, c("seqnames", "start", "end", "peak_id")]
df_clean_unique <- df_clean[!duplicated(df_clean$peak_id), ]
# peak coord file
coords <- df_clean_unique[, c("peak_id", "seqnames", "start", "end")]
colnames(coords) <- c("GENE", "CHR", "START", "END")
coords$CHR <- sub("^chr", "", coords$CHR)
fwrite(coords, "Peak_coord_hg19.txt", sep="\t")
