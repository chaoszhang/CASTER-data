setwd("C:/OneDrive - UC San Diego/data/CASTER/bird")
require(data.table); require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra); require(purrr)

d = read.csv("sliding.tsv", sep="\t", stringsAsFactors = TRUE)
d2 = d %>% group_by(branch, chrbin, group) %>% summarize(cnt = length(Q), A = sum(A) / sum(Q), B = sum(B) / sum(Q), C = sum(C) / sum(Q))
d2 = d2[d2$cnt > 2,]
dmelt = melt(d2, 1:4)
dmelt$variable = factor(dmelt$variable, levels = c("A", "C", "B"))
dmelt$chr = factor(dmelt$chr, levels = c(1:3, "4*", 5:10, "11-28*","Z*"))
dmelt$value = ifelse(dmelt$value > 0, dmelt$value, 0)
ggplot(aes(y=value,x=group+2.5,fill=variable),data=dmelt)+
  facet_grid(branch~chr,scales="free", space='free_x')+
  geom_hline(yintercept=0.333,linetype=2,size=0.4) +
  geom_hline(yintercept=0.667,linetype=2,size=0.4) +
  geom_bar(position="fill", stat="identity")+ # TODO add color
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(name="", breaks = c(0,0.33,0.67,1)) +
  scale_fill_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#D7E4F0","#FDB96B","#DA382A"),name="") + #TODO
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines")) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="B")
ggsave("avian_bar_chr.pdf", width = 9, height = 5)

d2$chrbin = factor(d2$chrbin, levels = c(1:3, "4*", 5:10, "11-28*","Z*"))
d2$rmv = "Kept"
d2[d2$chrbin %in% c("4*", "11-28*","Z*"),]$rmv = "Fully removed"
d2[interaction(d2$chrbin, d2$group) %in% c("1.0","1.190","10.0","10.15","2.0","2.145","3.0","3.100","5.0","5.55","7.5","7.30","8.5","8.25","9.0"),]$rmv = "Partially removed"
d2[interaction(d2$chrbin, d2$group) %in% c("3.105","6.30","7.0","8.0","9.20"),]$rmv = "Fully removed"
ggplot(aes(y=A+B+C,x=group+2.5,fill=rmv),data=d2)+
  facet_grid(branch~chrbin,scales="free", space='free_x')+
  geom_bar(stat="identity")+ # TODO add color
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(name="CASTER-site score") +
  scale_fill_manual(breaks = c("Kept", "Partially removed", "Fully removed"),values=c("#AAAAAA", "#E06060", "#FF2222"),name="") + 
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines")) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave("avian_score_chr.pdf", width = 9, height = 5)

g = read.csv("sliding_gc.tsv", sep="\t", stringsAsFactors = TRUE)
g2 = g %>% group_by(branch, chrbin, group) %>% summarize(cnt = length(gcA), A = sum(gcA) / sum(gcA + atA), B = sum(gcB) / sum(gcB + atB),
                                                         C = sum(gcC) / sum(gcC + atC), D = sum(gcD) / sum(gcD + atD))
g2 = g2[g2$cnt > 2,]
g2$A = g2$A - g2$D
g2$B = g2$B - g2$D
g2$C = g2$C - g2$D
gmelt = melt(g2, 1:4)
gmelt = gmelt[gmelt$variable != "D",]
gmelt$chr = factor(gmelt$chr, levels = c(1:10, "11-28","Z"))
ggplot(aes(y=value,x=group,color=variable),data=gmelt)+
  facet_grid(branch~chr,scales="free", space='free_x')+
  geom_hline(yintercept=0,linetype=2,size=0.4) +
  geom_line()+
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(name="") +
  scale_color_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#5F93C3","#FDB96B","#DA382A"),name="") + #TODO
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines")) #, strip.text.y = element_blank()
ggsave("avian_gc_chr.pdf", width = 9, height = 5)


(ggplot(aes(y=value,x=group+2.5,fill=variable),data=dmelt)+
    facet_grid(branch~chr,scales="free", space='free_x')+
    geom_hline(yintercept=0.333,linetype=2,size=0.4) +
    geom_hline(yintercept=0.667,linetype=2,size=0.4) +
    geom_bar(position="fill", stat="identity")+ # TODO add color
    theme_classic()+
    scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
    scale_y_continuous(name="", breaks = c(0,0.33,0.67,1)) +
    scale_fill_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#D7E4F0","#FDB96B","#DA382A"),name="") + #TODO
    theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))) + 
(ggplot(aes(y=value,x=group,color=variable),data=gmelt)+
   facet_grid(branch~chr,scales="free", space='free_x')+
   geom_hline(yintercept=0,linetype=2,size=0.4) +
   geom_line()+
   theme_classic()+
   scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
   scale_y_continuous(name="") +
   scale_color_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#5F93C3","#FDB96B","#DA382A"),name="") + #TODO
   theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))
)
