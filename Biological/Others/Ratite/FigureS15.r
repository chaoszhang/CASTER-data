setwd("C:/OneDrive - UC San Diego/data/CASTER/Ratite")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra)
require(stats)

d = read.csv("res.sliding", sep="\t", stringsAsFactors = TRUE)
d$TDR = d$TD.R / d$QuartetCnt
d$ACR = d$AC.R / d$QuartetCnt
d$ACTD = d$AC.TD / d$QuartetCnt
d$sum = d$TDR + d$ACR + d$ACTD
d = d[d$QuartetCnt >= 1000 * 60 & d$sum != 0,]

dCNEE = d[d$type == "CNEE",]
dIntron = d[d$type == "intron",]
dUCE = d[d$type == "UCE",]

dtemp = dCNEE
ggplot(dtemp, aes(x=sum)) + theme_classic()+
  geom_histogram(aes(y=..density..), bins = 100, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = median(dtemp$sum) - 5 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) - 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 5 * mad(dtemp$sum))

dtemp = dIntron
ggplot(dtemp, aes(x=sum)) + theme_classic()+
  geom_histogram(aes(y=..density..), bins = 100, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = median(dtemp$sum) - 5 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) - 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 5 * mad(dtemp$sum))

dtemp = dUCE
ggplot(dtemp, aes(x=sum)) + theme_classic()+
  geom_histogram(aes(y=..density..), bins = 100, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = median(dtemp$sum) - 5 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) - 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 3 * mad(dtemp$sum)) +
  geom_vline(xintercept = median(dtemp$sum) + 5 * mad(dtemp$sum))

dCNEE$absz = abs(dCNEE$sum - median(dCNEE$sum)) / mad(dCNEE$sum)
dIntron$absz = abs(dIntron$sum - median(dIntron$sum)) / mad(dIntron$sum)
dUCE$absz = abs(dUCE$sum - median(dUCE$sum)) / mad(dUCE$sum)

d = rbind(dCNEE, dIntron, dUCE)

dmelt = melt(d, c(1,2))
dmelt = dmelt[dmelt$variable %in% c("TDR", "ACR", "ACTD"),]
  
ggplot(dmelt, aes(x=value, y=variable, color=variable)) + 
  facet_wrap(type~., ncol = 1, scales = "free", strip.position = "right")+theme_classic()+
  geom_boxplot()+
  scale_y_discrete(name="Topology", labels=c("Tinamiformes+Rheiformes", "Apterygiformes+Casuariiformes+Rheiformes", "Apterygiformes+Casuariiformes+Tinamiformes"))+
  scale_x_continuous(name="CASTER score") +
  theme(legend.position="none", plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="A")
ggsave(file="ratite_outliers.pdf", width = 8, height = 4.5)

dmelt2 = melt(d[d$absz < 3,], c(1,2))
dmelt2 = dmelt2[dmelt2$variable %in% c("TDR", "ACR", "ACTD"),]
medians <- dmelt %>% group_by(variable, type) %>% summarize(median_value = median(value))
ggplot(dmelt2, aes(x=value, color=variable)) + 
  facet_wrap(.~type, scales = "free")+theme_classic()+
  geom_density()+
  geom_vline(data = medians, aes(xintercept = median_value, color=variable), linetype = "dashed")

temp = d[d$absz < 3,]
write.csv(temp, "included.csv", row.names=FALSE)


