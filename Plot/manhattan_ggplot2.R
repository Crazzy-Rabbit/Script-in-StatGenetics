library(data.table)
library(ggplot2)
library(dplyr)

dt = fread(infile)[, c("CHR", "BP", "P")]
# 1. 计算每个染色体的长度（BP最大值）
chr_lengths <- dt %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP))

# 2. 累计基因组位置
data <- dt %>%
  arrange(CHR, BP) %>%
  mutate(chr_cumsum = cumsum(c(0, head(chr_lengths$chr_len, -1)))[CHR]) %>%
  mutate(BPcum = BP + chr_cumsum)

# 3. 计算每个染色体在x轴上的中间点，用来显示染色体名称
axis_set <- data %>%
  group_by(CHR) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)

# 4. 绘制曼哈顿图
p = ggplot(data, aes(x=BPcum, y=-log10(P), color=factor(CHR))) +
  geom_point(alpha=0.9) +
  geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed") +
  scale_color_manual(values=rep(c("#1B2C62", "#4695BC"), 22)) + 
  scale_x_continuous(label=axis_set$CHR, 
  	                 breaks=axis_set$center, 
                     expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) + 
  labs(x="Chromosome", y=expression(-log[10](P))) +
  theme_classic() +
  theme(legend.position="none",
  	    axis.text.x = element_text(face="bold", size=20, color="black"),
        axis.text.y = element_text(face="bold", size=20, color="black"),
        axis.title.x = element_text(face="bold", size=20, margin=margin(t=10), color="black"),
        axis.title.y = element_text(face="bold", size=20, margin=margin(t=10), color="black"),
        axis.line = element_line(color="black"), 
        axis.ticks = element_line(color="black"),
        axis.ticks.length = unit(0.3, "cm"))
ggsave("my.png", p, height = 8, width = 23, dpi=300)