setwd("C:/OneDrive - UC San Diego/data/CASTER")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra)
mycolor = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "SVDQuartets" = "#E7298A",
            "wASTRAL" = "#D95F02")
myshape = c("CASTER-site" = 17, "CASTER-pair" = 19, "RAxML-ng" = 15, "SVDQuartets" = 18,
            "wASTRAL" = 88)
mainMethods = c("CASTER-site", "CASTER-pair", "RAxML-ng", "SVDQuartets", "wASTRAL")

subsample = c(3,4,11,15,16,22,23,26,27,28,29,31,32,36,38,40,42,43,45,50)

d = read.csv("caster.tsv", sep="\t", stringsAsFactors = TRUE)
d$condition = ""
d[d$data == "M201" & d$mut == 5e-8 & d$pop == 1e5 & d$ploidy == 1,]$condition = "Default"
d[d$data == "M201" & d$pop == 1e6,]$condition = "10X population size"
d[d$data == "M201" & d$mut == 5e-9,]$condition = "0.1X mutation rate"
d[d$data == "M201" & d$mut == 5e-7,]$condition = "10X mutation rate"
d[d$data == "M201" & d$ploidy == 2,]$condition = "Diploid (unphased)"
d$condition = factor(d$condition, c("Default", "10X population size", "0.1X mutation rate", "10X mutation rate", "Diploid (unphased)", ""))
d$method = factor(d$method, c("CASTER-site","CASTER-pair","RAxML-ng","RAxML-ng (GTGTR4)","SVDQuartets","wASTRAL","ASTRAL"))

s = d %>% group_by(data, mut, pop, len, ind, ploidy, method, condition) %>% summarise(rf = mean(rf), time = mean(time), core = mean(core))
s$corehour = s$time * s$core / 3600

pa = ggplot(aes(y=rf,x=time*core/3600,color=method,shape=method),data=d[d$data == "M201" & d$method %in% mainMethods & d$condition == "Default",])+
  geom_point(alpha=0.1)+
  geom_point(data=s[s$data == "M201" & s$method != "wASTRAL*" & s$condition == "Default",], size=5)+
  #facet_wrap(V3~.,scales="free")+
  theme_classic()+
  scale_color_manual(values = mycolor, name="")+
  scale_shape_manual(values = myshape, name="")+
  scale_x_log10(name="Running time (core hours)")+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 21, 2)) + 
  #coord_cartesian(ylim = c(0, 6)) +
  theme(legend.position=c(0.2,0.78),panel.border = element_rect(fill=NA), legend.text = element_text(size = 8), 
        legend.title = element_blank(), legend.key.size = unit(2, 'pt'))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="A"); pa

pb = ggplot(aes(y=rf,x=condition,fill=method),data=d[d$data == "M201" & (d$method %in% mainMethods | d$method == "RAxML-ng (GTGTR4)"),])+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "RAxML-ng (GTGTR4)" = "#40A0E3",
                               "SVDQuartets" = "#E7298A", "wASTRAL" = "#D95F02"), name="")+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 21, 2)) + 
  theme(legend.position=c(0.14,0.75),panel.border = element_rect(fill=NA), 
        legend.title = element_blank(), legend.key.size = unit(12, 'pt'))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="B"); pb

spa = ggplot(aes(y=rf,x=condition,fill=method),data=d[d$data == "M201",])+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "RAxML-ng (GTGTR4)" = "#40A0E3",
                               "SVDQuartets" = "#E7298A", "wASTRAL" = "#D95F02", "ASTRAL" = "#E6AB02"), name="")+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 30, 5)) + 
  theme(legend.position=c(0.85, 0.6),panel.border = element_rect(fill=NA),
        legend.title = element_blank(), legend.key.size = unit(16, 'pt'))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"), 
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="A"); spa

spb = ggplot(aes(y=time*core/3600,x=condition,fill=method),data=d[d$data == "M201",])+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "RAxML-ng (GTGTR4)" = "#40A0E3",
                               "SVDQuartets" = "#E7298A", "wASTRAL" = "#D95F02", "ASTRAL" = "#E6AB02"), name="")+
  scale_y_log10(name="Running time (core hours)")+
  scale_x_discrete(name="")+
  theme(legend.position="none",panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="B"); spb

pd = ggplot(aes(y=rf,x=as.factor(ind),fill=method),data=d[d$data == "M21" & d$method%in% mainMethods,])+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = mycolor, name="")+
  scale_x_discrete(name="# individuals per species")+
  scale_y_continuous(name="Species tree error (FN)") + 
  theme(legend.position="None",panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="D"); pd

spd = ggplot(aes(y=rf,x=as.factor(ind),fill=method),data=d[d$data == "M21",])+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "RAxML-ng (GTGTR4)" = "#40A0E3",
                               "SVDQuartets" = "#E7298A", "wASTRAL" = "#D95F02", "ASTRAL" = "#E6AB02"), name="")+
  scale_x_discrete(name="# individuals per species")+
  scale_y_continuous(name="Species tree error (FN)") + 
  theme(legend.position="None",panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="D"); spd

spe = ggplot(aes(y=time*core/3600,x=as.factor(ind),color=method,shape=method,group=method),data=d[d$data == "M21" & d$method %in% mainMethods & d$pop == 1e6,])+
  stat_summary(fun = mean, geom = "point", size = 1.5) + 
  stat_summary(fun = mean, geom = "line") + 
  stat_summary(geom="errorbar",width=0.1)+
  theme_classic()+
  scale_color_manual(values = mycolor, name="")+
  scale_shape_manual(values = myshape, name="")+
  scale_y_log10(name="Running time (core hours)")+
  scale_x_discrete(name="# individuals per species")+
  theme(legend.position=c(0.75,0.75),panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"), legend.title = element_blank(), 
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="E"); spe

t.test(rf ~ condition, data = d[d$condition %in% c("Default","0.1X mutation rate") & d$method == "wASTRAL",], paired = TRUE, alternative = "less")
t.test(rf ~ condition, data = d[d$condition %in% c("Default","10X mutation rate") & d$method == "wASTRAL",], paired = TRUE, alternative = "less")


t.test(rf ~ method, data = d[d$condition == "Default" & d$method %in% c("CASTER-site","wASTRAL") & d$rep %in% subsample,], paired = TRUE, alternative = "less")

t.test(rf ~ method, data = d[d$condition == "Diploid (unphased)" & d$method %in% c("CASTER-pair","RAxML-ng"),], paired = TRUE, alternative = "less")

t.test(rf ~ method, data = d[d$condition == "Default" & d$method %in% c("CASTER-site","CASTER-pair"),], paired = TRUE, alternative = "greater")
t.test(rf ~ method, data = d[d$condition == "Default" & d$method %in% c("RAxML-ng","CASTER-site"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d[d$condition == "10X population size" & d$method %in% c("CASTER-site","CASTER-pair"),], paired = TRUE, alternative = "greater")
t.test(rf ~ method, data = d[d$condition == "10X population size" & d$method %in% c("SVDQuartets","CASTER-site"),], paired = TRUE, alternative = "less")
summary(aov(rf~method*pop+Error(rep),data=d[d$condition %in% c("Default", "10X population size") & d$method != "wASTRAL*",]))
summary(aov((core*time)~method*pop+Error(rep),data=d[d$condition %in% c("Default", "10X population size") & d$method != "wASTRAL*",]))

summary(aov(rf~method*as.factor(mut)+Error(rep),data=d[(d$condition == "Default" | d$mut != 5e-8) & d$method != "wASTRAL*",]))
summary(aov((core*time)~method*as.factor(mut)+Error(rep),data=d[(d$condition == "Default" | d$mut != 5e-8) & d$method != "wASTRAL*",]))

t.test(rf~condition, data = d[d$condition %in% c("Default", "0.1X mutation rate") & d$method == "CASTER-site",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = d[d$condition %in% c("Default", "0.1X mutation rate") & d$method == "CASTER-pair",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = d[d$condition %in% c("Default", "0.1X mutation rate") & d$method == "RAxML-ng",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = d[d$condition %in% c("Default", "0.1X mutation rate") & d$method == "SVDQuartets",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = d[d$condition %in% c("Default", "0.1X mutation rate") & d$method == "wASTRAL",], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d[d$condition == "0.1X mutation rate" & d$method %in% c("RAxML-ng","CASTER-pair"),], paired = TRUE, alternative = "less")

dsub = d[d$rep %in% subsample,]
ssub = dsub %>% group_by(data, mut, pop, len, ind, ploidy, method, condition) %>% summarise(rf = mean(rf), time = mean(time), core = mean(core))
mean(dsub[dsub$mut == 5e-07 & dsub$method == "CASTER-site",]$rf) / mean(dsub[dsub$condition == "Default" & dsub$method == "CASTER-site",]$rf)
mean(dsub[dsub$mut == 5e-07 & dsub$method == "CASTER-pair",]$rf) / mean(dsub[dsub$condition == "Default" & dsub$method == "CASTER-pair",]$rf)
mean(dsub[dsub$mut == 5e-07 & dsub$method == "RAxML-ng",]$rf) / mean(dsub[dsub$condition == "Default" & dsub$method == "RAxML-ng",]$rf)
mean(dsub[dsub$mut == 5e-07 & dsub$method == "SVDQuartets",]$rf) / mean(dsub[dsub$condition == "Default" & dsub$method == "SVDQuartets",]$rf)
mean(dsub[dsub$mut == 5e-07 & dsub$method == "wASTRAL",]$rf) / mean(dsub[dsub$condition == "Default" & dsub$method == "wASTRAL",]$rf)
t.test(rf~condition, data = dsub[(dsub$mut == 5e-07 | dsub$condition == "Default") & dsub$method == "CASTER-site",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = dsub[(dsub$mut == 5e-07 | dsub$condition == "Default") & dsub$method == "CASTER-pair",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = dsub[(dsub$mut == 5e-07 | dsub$condition == "Default") & dsub$method == "RAxML-ng",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = dsub[(dsub$mut == 5e-07 | dsub$condition == "Default") & dsub$method == "SVDQuartets",], paired = TRUE, alternative = "less")
t.test(rf~condition, data = dsub[(dsub$mut == 5e-07 | dsub$condition == "Default") & dsub$method == "wASTRAL",], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = dsub[dsub$mut == 5e-07 & dsub$method %in% c("RAxML-ng","CASTER-pair"),], paired = TRUE, alternative = "less")


t.test(rf~condition, data = d[d$condition %in% c("Default", "Diploid (unphased)") & d$method == "CASTER-pair",], paired = TRUE, alternative = "greater")
t.test(rf~condition, data = d[d$condition %in% c("Default", "Diploid (unphased)") & d$method == "CASTER-site",], paired = TRUE, alternative = "greater")
t.test(rf~condition, data = d[d$condition %in% c("Default", "Diploid (unphased)") & d$method == "RAxML-ng",], paired = TRUE, alternative = "greater")
t.test(rf~condition, data = d[d$condition %in% c("Default", "Diploid (unphased)") & d$method == "SVDQuartets",], paired = TRUE, alternative = "greater")
t.test(rf~condition, data = d[d$condition %in% c("Default", "Diploid (unphased)") & d$method == "wASTRAL",], paired = TRUE, alternative = "greater")
t.test(rf~method, data = d[d$condition == "Diploid (unphased)" & d$method %in% c("RAxML-ng (GTGTR4)", "RAxML-ng"),], paired = TRUE, alternative = "greater")


t.test(rf~ind, data = d[d$data == "M21" & d$pop == 1e6 & d$method == "CASTER-site" & d$ind %in% c(5, 20),], paired = TRUE, alternative = "greater")
t.test(rf~ind, data = d[d$data == "M21" & d$pop == 1e6 & d$method == "CASTER-pair" & d$ind %in% c(5, 20),], paired = TRUE, alternative = "greater")

t.test(rf~ind, data = d[d$data == "M21" & d$method == "CASTER-site" & d$ind %in% c(1, 5),], paired = TRUE, alternative = "greater")
t.test(rf~ind, data = d[d$data == "M21" & d$method == "CASTER-pair" & d$ind %in% c(1, 5),], paired = TRUE, alternative = "greater")
t.test(rf~ind, data = d[d$data == "M21" & d$method == "RAxML-ng" & d$ind %in% c(1, 5),], paired = TRUE, alternative = "greater")
dsub = d[!d$rep %in% c(1,48),]
t.test(rf~ind, data = dsub[dsub$data == "M21" & dsub$method == "SVDQuartets" & dsub$ind %in% c(1, 5),], paired = TRUE, alternative = "greater")
t.test(rf~ind, data = d[d$data == "M21" & d$method == "wASTRAL" & d$ind %in% c(1, 5),], paired = TRUE, alternative = "greater")

spf = ggplot(aes(rf,color=method,shape=method),data=d[d$data == "M201" & (d$method %in% mainMethods | d$method == "ASTRAL") & d$condition == "Default",])+
  stat_ecdf()+theme_classic()+
  scale_color_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "SVDQuartets" = "#E7298A",
                                "wASTRAL" = "#D95F02", "ASTRAL" = "#E6AB02"), name="")+
  scale_x_continuous(name="Species tree error (FN)", breaks=seq(0, 16, 2))+
  scale_y_continuous(name="ECDF") +
  theme(legend.position=c(0.7,0.4),panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="F"); spf

d2 = read.csv("caster-scale.tsv", sep="\t", stringsAsFactors = TRUE)
d2$method = factor(d2$method, levels = c("CASTER-site", "CASTER-pair", "SVDQuartets", "RAxML-ng"))
ph = ggplot(aes(y=as.factor(core),x=method,fill=size,label=paste(size,"Mb")),data=d2)+
  geom_tile()+
  theme_classic()+
  #facet_grid(core~.)+
  #geom_bar(stat="identity")+theme_bw()+
  scale_fill_gradient(trans = "log2", breaks = c(1,4,32,256,2048), 
                      high = "#233B53", low = "#56B1F7",guide="none")+
  scale_x_discrete(name="")+
  scale_y_discrete(name="# cores")+geom_text(color="yellow") +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="H"); ph

d3 = read.csv("support.stat", sep="\t", stringsAsFactors = T, header = F)
d3$V3 = factor(d3$V3, c("CASTER-site","CASTER-pair","RAxML-ng (LRT)","wASTRAL"))
pf = ggplot(aes(V1, y = after_stat(y)*20,color=V3),data=d3[d3$V2 == "Correct" & d3$V3 %in% c("CASTER-site", "CASTER-pair"),])+
  stat_ecdf()+theme_classic()+
  scale_shape_manual(values = myshape, name="")+
  geom_line(aes(V1, y = 1 - after_stat(y)), data=d3[d3$V2 == "Incorrect" & d3$V3 %in% c("CASTER-site", "CASTER-pair"),], stat='ecdf', linetype = "dashed") +
  scale_color_manual(values = mycolor, name="")+
  scale_x_continuous(name="Local Bootstrap support (%)") +
  scale_y_continuous(name="Survival function for\nincorrect branches (1-ECDF)", breaks = c(0,0.2,0.4,0.6,0.8,1),
                     sec.axis = sec_axis(~./20, name="ECDF for correct branches")) + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position="Bottom") +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
          legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="F"); pf
pg = ggplot(aes(m = V1, d = (V2 == "Correct"), color = V3, shape = V3), data = d3) + 
  geom_roc(cutoffs.at = c(0.9, 1, 90, 95, 99, 100), labelsize = 3) +
  scale_color_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "wASTRAL" = "#D95F02", "RAxML-ng (LRT)" = "#7570B3"), name="")+
  scale_shape_manual(values = c("CASTER-site" = 17, "CASTER-pair" = 19, "wASTRAL" = 88, "RAxML-ng (LRT)" = 15), name="")+
  scale_x_continuous(name="False positive rate (1-specificity)", breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_y_continuous(name="True positive rate (sensitivity)", limits = c(0.9,1), breaks = c(0.9,0.92,0.94,0.96,0.98,1))+
  theme_classic()+theme(legend.position=c(0.7,0.5)) +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size=14,face = "bold"), legend.title = element_blank(), 
         legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="G"); pg


d4 = read.csv("S101.tsv", sep = "\t", stringsAsFactors = T)
d4$method = factor(d4$method, c("CASTER-site","CASTER-pair","FastTree","SVDQuartets","wASTRAL","ASTRAL"))
s4 = d4 %>% group_by(len, method) %>% summarise(rf = mean(rf))

pe = ggplot(aes(y=rf,x=as.factor(len),group=method,color=method,shape=method),data=d4)+
  stat_summary(fun.y = mean, geom = "point", size = 1.5) + 
  stat_summary(fun.y = mean, geom = "line") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, alpha = 0.3)+
  #geom_line()+geom_point()+
  theme_classic()+
  scale_color_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "FastTree" = "#7570B3", "SVDQuartets" = "#E7298A",
                                          "wASTRAL" = "#D95F02"), name="")+
  scale_shape_manual(values = c("CASTER-site" = 17, "CASTER-pair" = 19, "FastTree" = 15, "SVDQuartets" = 18,
                                "wASTRAL" = 88), name="")+
  scale_x_discrete(name="Gene length (bp)") +
  scale_y_continuous(name="Species tree error (FN)") +
  theme(legend.position=c(0.67,0.52), ) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"), legend.title = element_blank(),
        legend.key.size = unit(10, 'pt'),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt")) + labs(tag="E"); pe

t.test(rf ~ method, data = d4[d4$len == 1600 & d4$method %in% c("CASTER-site","FastTree"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d4[d4$len == 1600 & d4$method %in% c("CASTER-pair","FastTree"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d4[d4$len == 1600 & d4$method %in% c("FastTree","wASTRAL"),], paired = TRUE, alternative = "greater")

t.test(rf ~ method, data = d4[d4$len == 200 & d4$method %in% c("CASTER-site","FastTree"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d4[d4$len == 200 & d4$method %in% c("CASTER-pair","FastTree"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d4[d4$len == 200 & d4$method %in% c("CASTER-site","wASTRAL"),], paired = TRUE, alternative = "less")
t.test(rf ~ method, data = d4[d4$len == 200 & d4$method %in% c("CASTER-pair","wASTRAL"),], paired = TRUE, alternative = "less")

dds = read.csv("additionalSim/downsample/caster_downsample.tsv", sep = '\t', stringsAsFactors = T)
dds$method = factor(dds$method, c("CASTER-site", "CASTER-pair", "RAxML-ng", "SVDQuartets", "wASTRAL", "ASTRAL"))

pc = ggplot(aes(y=rf, fill=method, x=as.factor(gene), group=method), data=dds[dds$method %in% mainMethods,])+
  stat_summary(geom = "bar", fun = mean,position = position_dodge(width=0.8),width=0.8, color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, alpha = 1,position = position_dodge(width=0.8))+
  theme_classic()+
  scale_fill_manual(values = mycolor, name="")+
  scale_x_discrete(name="Subsampled data", labels = c("2.5%", "5%" ,"10%", "20%", "100%"))+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 21, 2)) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"),
        legend.title = element_blank(), legend.position=c(0.8, 0.8))+labs(tag="C"); pc

spc = ggplot(aes(y=rf, fill=method, x=as.factor(gene), group=method), data=dds)+
  stat_summary(geom = "bar", fun = mean,position = position_dodge(width=0.8),width=0.8, color = "black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, alpha = 1,position = position_dodge(width=0.8))+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair" = "#1B9E77", "RAxML-ng" = "#7570B3", "SVDQuartets" = "#E7298A",
                               "wASTRAL" = "#D95F02", "ASTRAL" = "#E6AB02"), name="")+
  scale_x_discrete(name="Subsampled data", labels = c("2.5%", "5%" ,"10%", "20%", "100%"))+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 21, 2)) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"),
        legend.title = element_blank(), legend.position=c(0.8, 0.8))+labs(tag="C"); spc

d6 = read.csv("additionalSim/version/version.tsv", sep = '\t', stringsAsFactors = T)
d6$method = factor(d6$method, c("CASTER-site","CASTER-pair (ALL)","CASTER-pair (RY)","CASTER-pair (WS)"))

spg = ggplot(aes(y=rf,x="Default",fill=method),data=d6)+
  stat_summary(geom="bar",position = position_dodge(width = 0.8),width=0.8,color="black")+
  stat_summary(geom="errorbar",position = position_dodge(width = 0.8),color="black",width=0.4)+
  theme_classic()+
  scale_fill_manual(values = c("CASTER-site" = "#66A61E", "CASTER-pair (ALL)" = "#1B9E77", "CASTER-pair (RY)" = "#BBBBBB", "CASTER-pair (WS)" = "#666666"), name="")+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Species tree error (FN)", breaks=seq(0, 10, 0.5)) + 
  theme(legend.position="right",panel.border = element_rect(fill=NA))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="G"); spg


g = arrangeGrob(pa, pb, pc, pd, pe, pf, pg, ph, widths = c(3, 1, 2, 1, 3),
                layout_matrix = rbind(c(1, 1, 2, 2, 2), c(3, 3, 4, 4, 5), c(6, 7, 7, 8, 8)))
ggsave(file="Simulation.pdf", g, width = 10, height = 9)

g = arrangeGrob(spa,spb,spc,spd,spe,spf, widths = c(1.8,1), heights = c(1,1,1), layout_matrix = rbind(c(1,4), c(2,5), c(3,6)))
ggsave(file="Simulation_sup.pdf", g, width = 11, height = 9)






dq = read.csv("additionalSim/quartet/quartet.tsv", sep = '\t', stringsAsFactors = T)
dq$method = factor(dq$method, levels = c("BPP", "CASTER-site", "CASTER-pair", "RAxML-ng", "SVDQuartets", "wASTRAL", "ASTRAL"))
dq1 = dq[dq$rf == 1 & dq$case == "hard" & dq$len >= 1 & dq$method %in% c("CASTER-site", "CASTER-pair", "RAxML-ng", "SVDQuartets", "wASTRAL", "ASTRAL"),]
partS1 = ggplot(aes(x=method, fill=method, alpha=top), data=dq1)+theme_classic()+
  facet_grid(mut~as.factor(len),switch = "x")+
  geom_bar(color = "black")+
  scale_x_discrete(name="Sequence length (Mbps)")+
  scale_y_continuous(labels = function(x) percent(x/300),name="Species tree error")+
  scale_fill_manual(values=c("#66A61E","#1B9E77","#7570B3","#E7298A","#D95F02","#E6AB02"), name = "Method")+
  scale_alpha_discrete(name = "Species tree")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        panel.spacing.x = unit(15,"pt"))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt")) + labs(tag="A"); partS1

part1 = ggplot(aes(x=method, fill=method, alpha=top), data=dq1[dq1$len < 20 & dq1$mut == "1X" & dq1$method %in% c("CASTER-site", "CASTER-pair", "RAxML-ng"),])+theme_classic()+
  facet_grid(.~as.factor(len),switch = "x")+
  geom_bar(color = "black")+
  scale_x_discrete(name="Sequence length (Mbps)")+
  scale_y_continuous(labels = function(x) percent(x/300),name="Species tree error")+
  scale_fill_manual(values=c("#66A61E","#1B9E77","#7570B3"), name = "Method")+
  scale_alpha_manual(values=c(0.4,0.7,1), name = "Species tree")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        panel.spacing.x = unit(15,"pt"))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt")) + labs(tag="B"); part1


dq2 = dq[ dq$mut == "1X" & dq$len <= 1 & dq$method %in% c("CASTER-site", "CASTER-pair", "BPP"),]
part2 = ggplot(aes(x=as.factor(1000*len),y=rf, color=method,group=method,linetype=top), data=dq2)+theme_classic()+
  facet_grid(.~case, labeller = labeller(case = c(easy="x = 0.001", hard="x = 0.0001")))+
  scale_x_discrete(name="Sequence length (kbps)")+
  scale_y_continuous(labels = function(x) percent(x),name="Species tree error")+
  scale_color_manual(values=c("#E6AB02","#66A61E","#1B9E77"), name = "")+
  stat_summary(geom="line", fun = mean)+
  stat_summary(geom="errorbar",width=0.2)+
  theme(legend.position=c(0.15,0.75))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="D"); part2

partS3 = ggplot(aes(x=method, fill=method, alpha=top), data=dq2[dq2$rf == 1,])+theme_classic()+
  facet_grid(case~as.factor(len),switch = "x", labeller = labeller(case = c(easy="x = 0.001", hard="x = 0.0001")))+
  geom_bar(color = "black")+
  scale_x_discrete(name="Sequence length (Mbps)")+
  scale_y_continuous(labels = function(x) percent(x/300),name="Species tree error")+
  scale_fill_manual(values=c("#E6AB02","#66A61E","#1B9E77"), name = "Method")+
  scale_alpha_manual(values=c(0.4,0.7,1), name = "Species tree")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        panel.spacing.x = unit(15,"pt"))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt")) + labs(tag="C"); partS3


t.test(rf~method, data = dq[dq$case == "easy" & dq$method %in% c("CASTER-pair", "BPP") & dq$len == 0.05,], paired = TRUE, alternative = "greater")
t.test(rf~method, data = dq[dq$case == "easy" & dq$method %in% c("CASTER-pair", "BPP") & dq$len == 0.1,], paired = TRUE, alternative = "greater")
summary(aov(rf~method*as.factor(len)+Error(top*rep), data=dq[dq$case == "easy" & dq$method %in% c("CASTER-pair", "BPP"),]))
summary(aov(rf~method*as.factor(len)+Error(top*rep), data=dq[dq$case == "easy" & dq$method %in% c("CASTER-site", "BPP"),]))




dgf = read.csv("additionalSim/quartet/quartet_introgression.tsv", sep = '\t', stringsAsFactors = T)
dgf$total = dgf$AB.CD + dgf$AC.BD + dgf$BC.AD
dgf$AB.CD = dgf$AB.CD / dgf$total
dgf$AC.BD = dgf$AC.BD / dgf$total
dgf$BC.AD = dgf$BC.AD / dgf$total

dgfm = melt(dgf, 1:4)
dgfm = dgfm[dgfm$variable != "total",]

ggplot(dgfm, aes(x=top, y=value, color=variable))+theme_bw()+
  facet_grid(mut~method, scale="free_y")+
  geom_boxplot()+
  scale_color_discrete(name = "Gene tree topology", labels = c("AB|CD", "AC|BD", "BC|AD"))+
  geom_hline(yintercept=0.3, linetype="dashed")+
  geom_hline(yintercept=0.4, linetype="dashed")

dgfms = dgfm %>% group_by(mut, top, variable, method) %>% summarise(value = mean(value))
dgfms$diff = ifelse(sub("\\|","",dgfms$top)==sub("\\.","",dgfms$variable),dgfms$value-0.4,dgfms$value-0.3)
partS2 = ggplot(dgfms, aes(x=top, fill= diff, label=round(diff,2), y=variable))+
  facet_grid(mut~method)+
  theme_classic()+
  geom_tile()+
  geom_text(aes(color=(abs(diff)>0.35)), show.legend = FALSE)+
  scale_color_manual(values=c("black","white"))+
  scale_x_discrete(name = "Species tree")+
  scale_y_discrete(name = "Error in support by topology", labels = c("AB|CD", "AC|BD", "BC|AD"))+
  scale_fill_gradient2(name="Error",midpoint = 0)+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="B"); partS2

part3 = ggplot(dgfms[dgfms$mut == "1X",], aes(x=top, fill= diff, label=round(diff,2), y=variable))+
  facet_grid(method~.)+
  theme_classic()+
  geom_tile(show.legend = FALSE)+
  geom_tile(data = dgfms[c(1,2,9,13,14,18,25,20,27),], fill = NA, color = "black", size = 1.5) +
  geom_text(aes(color=(abs(diff)>0.35)), show.legend = FALSE)+
  scale_color_manual(values=c("black","white"))+
  scale_x_discrete(name = "Species tree")+
  scale_y_discrete(name = "Error in support by topology", labels = c("AB|CD", "AC|BD", "BC|AD"))+
  scale_fill_gradient2(name="Error",midpoint = 0)+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="C"); part3

t = data.frame(species=c("AB|CD","AB|CD","AB|CD","AC|BD","AC|BD","AC|BD","BC|AD","BC|AD","BC|AD"),
               gene=c("AB|CD","AC|BD","BC|AD","AB|CD","AC|BD","BC|AD","AB|CD","AC|BD","BC|AD"),
               freq=c(0.4,0.3,0.3,0.3,0.4,0.3,0.3,0.3,0.4))

part0 = ggplot(t, aes(x=species, y=gene, label=freq, fill = freq))+
  theme_classic()+
  geom_tile(show.legend = FALSE)+
  geom_text(aes(color=(freq>0.35)), show.legend = FALSE)+
  scale_color_manual(values=c("black","white"))+
  scale_x_discrete(name = "Species tree")+
  scale_y_discrete(name = "Gene tree frequency   ", labels = c("AB|CD", "AC|BD", "BC|AD"))+
  scale_fill_gradient(low = "white", high = "black")+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt")); part0

g = arrangeGrob(part1, part2, part3, part0, ggplot()+theme_bw()+
                  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
                        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="A"), widths = c(2.5, 2.5, 2.5), heights = c(2, 2.5, 2.5),
                layout_matrix = rbind(c(5, 5, 4), c(1, 1, 3), c(2, 2, 3)))
ggsave(file="Quartet.pdf", g, width = 7.5, height = 7)

g = arrangeGrob(partS1, partS2, partS3, widths = c(1, 1), heights = c(1, 1),
                layout_matrix = rbind(c(1, 1), c(2, 3)))
ggsave(file="Quartet_sup.pdf", g, width = 12, height = 10)
