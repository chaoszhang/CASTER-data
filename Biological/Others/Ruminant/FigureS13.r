setwd("C:/OneDrive - UC San Diego/data/CASTER")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra); require(tidyquant); require(MASS)

d = read.csv("Ruminant/sliding.tsv", sep="\t", stringsAsFactors = TRUE)
d$Topology = factor(d$Topology, c("Bongo+Sitatunga", "Bongo+MountainNyala", "Sitatunga+MountainNyala"))
g1 = ggplot(aes(y=Score,x=Pos,color=Topology), data = d) +
  geom_point(aes(alpha = Score > 1, size = Score > 1))+
  geom_ma(aes(x=Pos-2.5), ma_fun = SMA, n = 50, linetype="solid") +
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*100), expand = c(0.01, 0.01))+
  scale_y_continuous(name="") +
  scale_alpha_manual(guide="none",values=c(0.1, 1))+
  scale_size_manual(guide="none",values=c(0.5, 1))+
  theme_classic() +
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))
  
d2 = read.csv("ButterflyPKU/butterfly.sliding.tsv", sep="\t", stringsAsFactors = TRUE)
d2 = d2[d2$Chr != 32,]
r = sum(d2$Score) / nrow(d2[d2$Topology == "inachus+paralekta",])
g2 = ggplot(aes(y=Score/r,x=Pos/1000000,color=Topology), data = d2) +
  geom_point(aes(alpha = Score > 1, size = Score > 1))+
  #geom_point(alpha = 0.3)+
  facet_grid(.~Chr,scales="free", space='free_x')+
  geom_ma(aes(x=Pos/1000000-0.25), ma_fun = SMA, n = 5, linetype="solid") +
  scale_x_continuous(name="Position (million SNPs)", breaks = c(0,1,2,3), labels = c(" 0","1","2 ","3 "), expand = c(0, 0))+
  scale_y_continuous(name="") +
  scale_alpha_manual(guide="none",values=c(0.25, 1))+
  scale_size_manual(guide="none",values=c(0.5, 1))+
  theme_classic() +
  theme(legend.position="bottom", panel.border = element_rect(color = "#AAAAAA", fill=NA, linewidth=0.5), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))

g = arrangeGrob(g1, g2, widths = c(1,2), layout_matrix = rbind(c(1, 2)))
ggsave(file="sliding.pdf", g, width = 18, height = 5.5)
