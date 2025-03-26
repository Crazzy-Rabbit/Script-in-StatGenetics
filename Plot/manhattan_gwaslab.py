import gwaslab as gl
import pandas as pd
import matplotlib.pyplot as plt


# 1. data read
gwas = gl.Sumstats("ARHL_MVP_AJHG_BBJ_reformatMETAL_addchr.gz", 
                   snpid="SNP", 
                   chrom="CHR",
                   pos="POS",
                   ea="A1",
                   nea="A2",
                   eaf="freq",
                   beta="beta",
                   se="SE",
                   p="p",
                   n="N",
                   build="19")

# 2. plot manhattan with annota
df = pd.read_csv("novel_snp_ARHL.txt", sep="\t")
anno_list = df["SNP"].tolist()

gwas.plot_mqq(mode="m", skip=0, sig_line_color="red", fontsize=12, marker_size=(5,5),
              anno="GENENAME", anno_style="expand", anno_set=anno_list, anno_fontsize=12, repel_force=0.01, arm_scale=1,
              xtight=True, ylim=(0,38), chrpad=0.01, xpad=0.05, # cut=40, cut_line_color="white",
              fig_args={"figsize": (18, 5), "dpi": 500},
              save="mqq_plot.png", save_args={"dpi": 500}, check=False, verbose=False)