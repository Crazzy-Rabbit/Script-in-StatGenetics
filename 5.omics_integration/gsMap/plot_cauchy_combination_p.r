# # # # # plot cauchy combination p value
library(ggplot2)
library(data.table)

data <- fread("gsMap_cauchy_combination_result.txt")
data <- data[order(data$P_Cauchy, decreasing=TRUE),]

# color for dif tissue
custom_colors <- c("Hair cells" = "#DD7694", 
                "Supporting cells" = "#BCD4E7", 
                "Surrounding structures" = "#056E83", 
                "Lateral wall" = "#E9D9BF", 
                "Circulating cells" = "#D4920A", 
                "Glial cells" = "#5AA4AE", 
                "Neurons" = "#65472F")   

# set Name as factor
data$Annotation <- factor(data$Annotation, levels=unique(data$Annotation))

p=ggplot(data, aes(x=Annotation, y=-log10(P_Cauchy), fill=Annotation)) +
    geom_bar(stat="identity", position="dodge", width=0.8, fill="#056E83")+
    scale_y_continuous(expand=c(0,0), limits=c(0, 9))+
    labs(x=NULL, y="-log10(P Cauchy)", title="")+
    coord_flip()+
    theme_bw()+
    theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="none",
        legend.title=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.text = element_blank(),
        strip.background=element_blank())
ggsave("plot cauchy combination p value.png", p, width=3.5, height=5, dpi=500)
