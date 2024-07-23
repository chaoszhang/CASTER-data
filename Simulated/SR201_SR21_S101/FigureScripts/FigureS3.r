setwd("C:/OneDrive - UC San Diego/data/CASTER/additionalSim/score_true")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra)
require(stats)

d = read.csv("caster-w-true-tree.tsv", sep="\t", stringsAsFactors = TRUE)
ggplot(d, aes(y=score_reconstructed-score_true, color=method))+theme_classic()+
  facet_wrap(.~method, scales = "free")+
  geom_boxplot() + 
  scale_color_discrete(name = "") +
  scale_y_continuous(name="") +
  scale_x_discrete(name="Score difference (Reconstructed - True)")
ggsave(file="score_true.pdf", width = 7, height = 5)
