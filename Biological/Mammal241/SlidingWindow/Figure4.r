setwd("C:/OneDrive - UC San Diego/data/CASTER/mammal")
require(data.table); require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(ggbeeswarm); require(plotROC); require(gridExtra); require(tidyquant); require(MASS); require(purrr)

d5 = read.csv("stat_test_outlier.tsv", sep="\t", stringsAsFactors = TRUE)
d5$iRow = 1
d5$realValue = 1
setorder(d5, Len)
d5[d5$Len == 1,]$realValue = p.adjust(10 ^ -d5[d5$Len == 1,]$Value, "BH", n=53452284*2)
d5[d5$Len == 5,]$realValue = p.adjust(10 ^ -d5[d5$Len == 5,]$Value, "BH", n=11400178*2)
d5[d5$Len == 20,]$realValue = p.adjust(10 ^ -d5[d5$Len == 20,]$Value, "BH", n=2983610*2)
d5[d5$Len == 100,]$realValue = p.adjust(10 ^ -d5[d5$Len == 100,]$Value, "BH", n=605634*2)
d5[d5$Len == 500,]$realValue = p.adjust(10 ^ -d5[d5$Len == 500,]$Value, "BH", n=121189*2)
d5 = d5[d5$realValue < 0.01,]
d5$Color = "Others"
d5[d5$Branch %in% c(196:205,207:237),]$Color = "Other primates"
d5[d5$Branch %in% 217:219,]$Color = "Presbytini monkeys"
d5[d5$Branch == 231,]$Color = "Macaques"
d5[d5$Branch == 21,]$Color = "Southern seals"
d5[d5$Branch %in% c(56,57,59,60),]$Color = "Caprini"
d5[d5$Branch %in% c(65,67),]$Color = "Cattle"
d5[d5$Branch %in% c(32:34,36:37),]$Color = "Cats"
d5[d5$Branch %in% c(28,29),]$Color = "Dogs"
d5[d5$Branch %in% c(97,98),]$Color = "Camelids"
d5[d5$Branch == 170,]$Color = "Cavies"
d5[d5$Branch %in% c(78:88,90:92),]$Color = "Whales"

d5$Color <- factor(d5$Color, levels = c("Macaques", "Presbytini monkeys", "Other primates", 
                                        "Cats", "Dogs", "Cattle", "Camelids",
                                        "Caprini", "Cavies",
                                        "Whales", "Southern seals", "Others"))
mycolors = c("#006060", "#60A0A0", "#C0E0E0", 
             "#FFB0B0", "#F08010", "#C05050", "#D0B050", "#905010", "#8060A0",
             "#5080C0", "#45C0F0", "#808080")

d5$Mpos = floor(d5$Pos / 10000000) * 10
s5 = d5[d5$Len == 1,] %>% group_by(Chr, Mpos, Color) %>% summarise(count = length(realValue))
s5$iRow = 2

temp = d5 %>% group_by(Color, Len) %>% summarise(count = length(realValue), color = Color[1])
ggplot()+
  facet_grid(iRow~Chr,scales="free", space='free_x', labeller = labeller(iRow = c(`1` = "-log10(P-value)", `2` = "Frequency (10kbp)")))+
  geom_hline(data = data.frame(iRow = c(1), h = c(-log10(0.01))), aes(yintercept = h), linetype="dashed") +
  geom_point(aes(y=-log10(realValue),x=Pos/1000000,color=Color,size=as.factor(Len),shape=as.factor(Len)), stroke=0.5, data=d5[d5$Len %in% c(5,20,100),])+
  #geom_segment(aes(y=-log10(realValue),yend=-log10(realValue),x=Pos/1000000,xend=Pos/1000000+10,color=Color,group=Branch), data=d5, linewidth=1.25)+
  geom_bar(aes(y=count,x=Mpos+5,fill=Color), position="stack", stat="identity", data=s5)+
  theme_classic()+
  scale_size_manual(name="", values = c(`5`=1, `20`=1, `100`=2))+
  scale_shape_manual(name="", values = c(`5`=4, `20`=19, `100`=10))+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:6)*100), expand = c(0, 10), limits = c(0, NA))+
  scale_y_continuous(name="") +
  scale_color_manual(name="", values = mycolors)+
  scale_fill_manual(name="", values = mycolors)+
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines")) +
  guides(color = guide_legend(nrow = 1))
ggsave("mammal_main.pdf", width = 15, height = 4.5)

temp1 = d5[d5$Len %in% c(5,20,100),]

d5 = read.csv("stat_test_outlier.tsv", sep="\t", stringsAsFactors = TRUE)
d5$Chr <- factor(d5$Chr, levels = c("1", "2", "3" ,"4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))
d5$iRow = 1
d5$Color = as.factor(ifelse(d5$Branch %in% c(21,170,217,231), d5$Branch, 0))
d5$Color = "Others"
d5[d5$Branch %in% c(224:225,227:231),]$Color = "Cercopithecinae"
d5[d5$Branch %in% c(217:220,222),]$Color = "Colobinae"
d5[d5$Branch %in% 234:237,]$Color = "Apes"
d5[d5$Branch %in% c(207:213,215),]$Color = "New World monkeys"

d5[d5$Branch %in% c(78:87),]$Color = "Toothed whales"
d5[d5$Branch == 21,]$Color = "Southern seals"
d5[d5$Branch == 97,]$Color = "Camels"
#d5[d5$Branch %in% c(56,57,59,60),]$Color = "Caprini"
d5[d5$Branch %in% c(32:34,36:37),]$Color = "Cats"
d5[d5$Branch == 170,]$Color = "Cavies"
#d5[d5$Branch %in% c(28,29),]$Color = "Dogs"
#d5[d5$Branch == 47,]$Color = "Horses"
#d5[d5$Branch %in% c(65,67),]$Color = "Cattle"

#d5$Color <- factor(d5$Color, levels = c("New World monkeys", "Colobinae", "Cercopithecinae", "Apes", 
#                                        "Cats", "Dogs", "Cattle", "Camels",
#                                        "Caprini", "Cavies", "Horses",
#                                        "Toothed whales", "Southern seals", "Others"))
d5$Color <- factor(d5$Color, levels = c("New World monkeys", "Colobinae", "Cercopithecinae", "Apes", 
                                        "Cats", "Camels", "Cavies",
                                        "Toothed whales", "Southern seals", "Others"))
#mycolors = c("#C0F0F0", "#80C0C0", "#409090", "#006060",
#             "#FFB0B0", "#F08010", "#C05050", "#D0B050", "#905010", "#8060A0", "#7080C0",
#             "#5080C0", "#45C0F0", "#808080")
mycolors = c("#C0F0F0", "#80B0B0", "#407070", "#003030",
             "#FFB0B0", "#F08010", "#905010", "#8060A0",
             "#5080C0", "#808080")
d5$Mpos = floor(d5$Pos / 10000000) * 10
d5$Coverage = ifelse(d5$relativeQ >= 1, "high", "low")
d5 = d5[d5$Category != "C" & d5$Coverage == "high",]
d5$realValue = p.adjust(10 ^ -d5$Value, "BH", n=17814773*2)
d5old = d5
d5 = d5[d5$realValue < 0.05,]
d5filter = d5 %>% group_by(Chr, Pos) %>% summarise(count = length(realValue))
#d5filter = d5filter[d5filter$count <= 1,]
d5 = d5[interaction(d5$Chr, d5$Pos) %in% interaction(d5filter$Chr, d5filter$Pos),]
s5 = d5 %>% group_by(Chr, Mpos, Color) %>% summarise(count = length(realValue))
s5$iRow = 2
temp1 = d5 %>% group_by(Branch) %>% summarise(count = length(realValue), color = Color[1])
temp2 = d5[d5$Chr != 'X',] %>% group_by(Branch) %>% summarise(count = length(realValue), color = Color[1])

ggplot()+
  facet_grid(iRow~Chr,scales="free", space='free_x', labeller = labeller(iRow = c(`1` = "-log10(P-value)", `2` = "Frequency (10M bin)")))+
  geom_hline(data = data.frame(iRow = c(1), h = c(-log10(0.05))), aes(yintercept = h), linetype="dashed") +
  geom_point(aes(y=-log10(realValue),x=Pos/1000000,color=Color,group=Branch), data=d5, linewidth=1.25)+
  #geom_segment(aes(y=-log10(realValue),yend=-log10(realValue),x=Pos/1000000,xend=Pos/1000000+10,color=Color,group=Branch), data=d5, linewidth=1.25)+
  #geom_bar(aes(y=count,x=Mpos+5,fill=Color), position="stack", stat="identity", data=s5)+
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:6)*50), expand = c(0, 15), limits = c(NA, NA))+
  scale_y_continuous(name="") +
  scale_color_manual(name="", values = mycolors)+
  scale_fill_manual(name="", values = mycolors)+
  theme(legend.position="bottom", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines")) +
  guides(color = guide_legend(nrow = 1))
ggsave("mammal_main.pdf", width = 15, height = 5)

d2 = read.csv("sliding.tsv", sep="\t", stringsAsFactors = TRUE)
d2$Chr <- factor(d2$Chr, levels = c("1", "2", "3" ,"4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))
d2$Topology <- factor(d2$Topology, levels = c("A", "C", "B"))
d2$Score <- ifelse(d2$Score > 0, d2$Score, 0)
pattern = c(10, 12, 152, 162, 187, 217, 237)
patternChr = c(1, 7, 9, 16:22, "X")
ggplot(aes(y=Score,x=Pos/1000000+1.25,fill=Topology),data=d2[d2$Branch %in% pattern & d2$Chr %in% patternChr & d2$Coverage == "high",])+ # & d$Chr != "X" &
  facet_grid(Branch~Chr,scales="free", space='free_x')+
  geom_hline(yintercept=0.333,linetype=2,linewidth=0.4) +
  geom_hline(yintercept=0.667,linetype=2,linewidth=0.4) +
  geom_bar(position="fill", stat="identity")+ # TODO add color
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(name="", breaks = c(0,0.33,0.67,1)) +
  scale_fill_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#D7E4F0","#FDB96B","#DA382A"),name="") + #TODO
  theme(legend.position="bottom", strip.text.y = element_blank(), panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))
ggsave("mammal_bar_chr.pdf", width = 9, height = 7)

#21: S. seal, 34: Cats, 170: Cavies, 208: Platyrrhini, 218: Presbytini, 231: Macaca
d3=d2[d2$Branch %in% c(21,28,87,218,231) & d2$Coverage == "high",]
d3$Branch <- factor(d3$Branch, levels = c(231,218,28,21,87))
ggplot(aes(y=Score,x=Pos/1000000+1.25,fill=Topology),data=d3[d3$Chr %in% c(1,3,4,6,7,9,10,12,16,19,20,21,22),])+ #
#ggplot(aes(y=Score+abs(Score),x=Pos/1000000,fill=Topology),data=d3)+ # & d$Chr != "X" &
  facet_grid(Branch~Chr,scales="free", space='free_x',
    labeller = labeller(Branch = c(`218` = "Presbytini", `231` = "Macaca", `21` = "Monachinae", `28` = "Caninae", `87` = "Odontoceti")))+
  geom_hline(yintercept=0.333,linetype=2,linewidth=0.4) +
  geom_hline(yintercept=0.667,linetype=2,linewidth=0.4) +
  geom_bar(position="fill", stat="identity")+ # TODO add color
  theme_classic()+
  scale_x_continuous(name="Position (Mbp)", breaks = c((0:5)*50), expand = c(0, 0), limits = c(0, NA))+
  scale_y_continuous(name="", breaks = c(0,0.33,0.67,1)) +
  scale_fill_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"), values=c("#D7E4F0","#FDB96B","#DA382A"),name="") + #TODO
  theme(legend.position="none", panel.border = element_rect(fill=NA), panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0.5, "lines"))
ggsave("mammal_bar_hgt.pdf", width = 15, height = 5)

d8 = read.csv("../sim/stat_test_all.tsv", sep="\t", stringsAsFactors = TRUE)
sumd8Count = sum(d8$Count)
d8$Ncount = d8$Count / sumd8Count
s8 = d8[0.5**d8$Bin>=1/sumd8Count,] %>% arrange(Bin) %>% mutate(sumNcount = map_dbl(Bin, ~ sum(Ncount[Bin > .x])), sumCount = map_dbl(Bin, ~ sum(Count[Bin > .x])))
( ggplot(d8, aes(x=0.5**Bin, y=Count)) + geom_bar(stat="identity") + theme_classic() +
  scale_x_log10(name = "Fitted p-value", lim = c(1e-6, 1e5)) + scale_y_log10(name = "Empirical frequency", lim = c(1, 1e4)) +
  geom_function(fun = function(x) 0.5 * sumd8Count * x/((1-2**-1.01)/(1-2**-0.01)), color = "red", linewidth = 2, linetype = 2) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="a)") ) +
( ggplot(s8, aes(x=0.5**Bin, y=sumNcount)) + theme_classic() + geom_point() + 
  scale_x_log10(name = "Fitted p-value", lim = c(1e-7, 1e2)) + scale_y_log10(name = "Empirical cumulative distribution function", lim = c(1e-7, 1e0)) + 
  geom_function(fun = function(x) x, color = "blue", linewidth = 2, linetype = 2) +
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"))+labs(tag="b)") )
ggsave("pq_plot.pdf", width = 14, height = 7)

d4 = read.csv("example1k.tsv", sep="\t", stringsAsFactors = TRUE)
ggplot(aes(y=Score,x=Pos/1000,color=Topology),data=d4) +
  facet_wrap(Branch~Chr, scale = "free", labeller = "label_both") +
  geom_hline(yintercept=0,linetype=2,linewidth=0.4) +
  #geom_point() +
  geom_segment(aes(y=Score,yend=Score,x=Pos/1000,xend=Pos/1000+10,color=Topology),data=d4, linewidth=1.5) + 
  theme_classic()+
  scale_x_continuous(name="Position (kbp)")+
  scale_y_continuous(name="Normalized CASTER score (10kbp window)", breaks = c(0)) +
  scale_color_manual(breaks = c("C", "A", "B"), labels = c("Main", "Alternative 1", "Alternative 2"),values=c("#AFC9E1","#FDB96B","#DA382A"),name="") + #TODO
  theme(legend.position="bottom", panel.border = element_rect(fill=NA))
ggsave("example1k.pdf", width = 14, height = 6)
