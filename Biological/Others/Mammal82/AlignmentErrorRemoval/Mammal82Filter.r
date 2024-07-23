setwd("C:/OneDrive - UC San Diego/data/CASTER/Mammal2")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra)
require(stats)

d = read.csv("whale.stat", sep="\t", stringsAsFactors = TRUE)
d$Other_toothed_whales.Minke_whale = d$Other_toothed_whales.Minke_whale / d$QuartetCnt
d$Sperm_whale.Minke_whale = d$Sperm_whale.Minke_whale / d$QuartetCnt
d$Sperm_whale.Other_toothed_whales = d$Sperm_whale.Other_toothed_whales / d$QuartetCnt
d$sum = d$Other_toothed_whales.Minke_whale + d$Sperm_whale.Minke_whale + d$Sperm_whale.Other_toothed_whales
d = d[d$QuartetCnt >= 255 * 200,]

d$absz = abs(d$sum - median(d$sum)) / mad(d$sum)

dmelt = melt(d, c(1))
dmelt = dmelt[dmelt$variable %in% c("Other_toothed_whales.Minke_whale", "Sperm_whale.Minke_whale", "Sperm_whale.Other_toothed_whales"),]
  
pa = ggplot(dmelt, aes(y=variable, x=value, color=variable))+theme_classic()+
  geom_boxplot() + 
  scale_color_discrete(name = "") +
  scale_y_discrete(name="Topology") +
  scale_x_continuous(name="CASTER score")+
  theme(legend.position="right",plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="A"); pa

dmelt2 = melt(d[d$absz < 10,], c(1))
dmelt2 = dmelt2[dmelt2$variable %in% c("Other_toothed_whales.Minke_whale", "Sperm_whale.Minke_whale", "Sperm_whale.Other_toothed_whales"),]

medians <- dmelt2 %>% group_by(variable) %>% summarize(median_value = median(value))
ggplot(dmelt2, aes(x=value, color=variable))+theme_classic()+
  geom_density(n=4096)+coord_cartesian(xlim = c(-2e-6, 2e-6))+
  geom_vline(data = medians, aes(xintercept = median_value, color=variable), linetype = "dashed")

rmv = d[d$absz > 10,c("seq","absz")]

d = read.csv("arctoidea.stat", sep="\t", stringsAsFactors = TRUE)
d$Seals.Giant_panda = d$Seals.Giant_panda / d$QuartetCnt
d$Ferret.Giant_panda = d$Ferret.Giant_panda / d$QuartetCnt
d$Ferret.Seals = d$Ferret.Seals / d$QuartetCnt
d$sum = d$Seals.Giant_panda + d$Ferret.Giant_panda + d$Ferret.Seals
d = d[d$QuartetCnt >= 172 * 200 & d$sum != 0,]

d$absz = abs(d$sum - median(d$sum)) / mad(d$sum)

dmelt = melt(d, c(1))
dmelt = dmelt[dmelt$variable %in% c("Seals.Giant_panda", "Ferret.Giant_panda", "Ferret.Seals"),]

ggplot(dmelt, aes(x=value, color=variable))+theme_classic()+
  geom_density()

pb = ggplot(dmelt, aes(y=variable, x=value, color=variable))+theme_classic()+
  geom_boxplot() + 
  scale_color_discrete(name = "") +
  scale_y_discrete(name="Topology") +
  scale_x_continuous(name="CASTER score")+
  theme(legend.position="right",plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="B"); pb

dmelt2 = melt(d[d$absz < 10,], c(1))
dmelt2 = dmelt2[dmelt2$variable %in% c("Seals.Giant_panda", "Ferret.Giant_panda", "Ferret.Seals"),]

medians <- dmelt2 %>% group_by(variable) %>% summarize(median_value = median(value))
ggplot(dmelt2, aes(x=value, color=variable))+theme_classic()+
  geom_density(n=4096)+coord_cartesian(xlim = c(-1e-5, 1e-5))+
  geom_vline(data = medians, aes(xintercept = median_value, color=variable), linetype = "dashed")

rmv = rbind(rmv, d[d$absz > 10,c("seq","absz")])

d = read.csv("paenungulata.stat", sep="\t", stringsAsFactors = TRUE)
d$Elephant.Rock_hyrax = d$Elephant.Rock_hyrax / d$QuartetCnt
d$Manatee.Rock_hyrax = d$Manatee.Rock_hyrax / d$QuartetCnt
d$Manatee.Elephant = d$Manatee.Elephant / d$QuartetCnt
d$sum = d$Elephant.Rock_hyrax + d$Manatee.Rock_hyrax + d$Manatee.Elephant
d = d[d$QuartetCnt >= 87 * 200 & d$sum != 0,]

d$absz = abs(d$sum - median(d$sum)) / mad(d$sum)

dmelt = melt(d, c(1))
dmelt = dmelt[dmelt$variable %in% c("Elephant.Rock_hyrax", "Manatee.Rock_hyrax", "Manatee.Elephant"),]

ggplot(dmelt, aes(x=value, color=variable))+theme_classic()+
  geom_density()

pc = ggplot(dmelt, aes(y=variable, x=value, color=variable))+theme_classic()+
  geom_boxplot() + 
  scale_color_discrete(name = "") +
  scale_y_discrete(name="Topology") +
  scale_x_continuous(name="CASTER score")+
  theme(legend.position="right",plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="C"); pc

dmelt2 = melt(d[d$absz < 10,], c(1))
dmelt2 = dmelt2[dmelt2$variable %in% c("Elephant.Rock_hyrax", "Manatee.Rock_hyrax", "Manatee.Elephant"),]

medians <- dmelt2 %>% group_by(variable) %>% summarize(median_value = median(value))
ggplot(dmelt2, aes(x=value, color=variable))+theme_classic()+
  geom_density(n=4096)+coord_cartesian(xlim = c(-1e-5, 1e-5))+
  geom_vline(data = medians, aes(xintercept = median_value, color=variable), linetype = "dashed")

rmv = rbind(rmv, d[d$absz > 10,c("seq","absz")])
write.csv(rmv$seq, "z10.list", row.names=FALSE)

pa / pb / pc
ggsave(file="mammal82_outliers.pdf", width = 10, height = 6)
