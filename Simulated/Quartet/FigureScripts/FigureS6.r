setwd("C:/OneDrive - UC San Diego/data/CASTER/additionalSim/quartet")
require(reshape2); require(ggplot2); require(scales); require(patchwork); require(dplyr); require(gridExtra)

d0 = read.csv("bppLogL.stat", sep = '\t', stringsAsFactors = T)
d1 = read.csv("bppMCMC.stat", sep = '\t', stringsAsFactors = T)

(ggplot(d0[d0$mix == "well",], aes(x=percent, color=start, y=logL))+
    facet_wrap(.~num, scales="free_y", labeller = labeller(num=c(`1`="(a) 0.2Mbps, Hard, Rep #14",`2`="(b) 0.1Mbps, Easy, Rep #40",`3`="(c) 1Mbps, Easy, Rep #1")))+
    theme_classic()+
    scale_x_continuous(name = "")+
    scale_y_continuous(name = "Log-Likelihood")+
    scale_color_discrete(name = "Starting tree")+
    geom_line() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + labs(tag="A")) / 
(ggplot(d1[d1$mix == "well",], aes(x=chain/10000, fill=start, alpha=as.factor(active), y=interaction(state)))+
  facet_grid(start~num, scales="free_y")+
  theme_classic()+
  geom_tile()+
  scale_alpha_manual(values = c(0, 1), guide = "none")+
  scale_fill_discrete(name = "Starting tree")+
  scale_x_continuous(name = "Million iterations (1 record per 100 MCMC iterations)")+
  scale_y_discrete(name = "MCMC topology", labels = c("AB|CD", "AC|BD", "BC|AD", "AB|CD", "AC|BD", "BC|AD", "AB|CD", "AC|BD", "BC|AD"))+
  theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"), panel.border = element_rect(fill="NA"),
        panel.spacing.y = unit(2,"pt"),panel.spacing.x = unit(50,"pt"),
        legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"), strip.text = element_blank())) / 
(ggplot(d0[d0$mix == "poor",], aes(x=percent, color=start, y=logL))+
  facet_wrap(.~num, scales="free_y", labeller = labeller(num=c(`4`="(d) 1Mbps, Hard, Rep #100", `5`="(e) 1Mbps, Hard, Rep #57", `6`="(f) 1Mbps, Hard, Rep #89")))+
  theme_classic()+
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "Log-Likelihood")+
  scale_color_discrete(name = "Starting tree")+
  geom_line() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(tag="B")) / 
(ggplot(d1[d1$mix == "poor",], aes(x=chain/2000, fill=start, alpha=as.factor(active), y=interaction(state)))+
   facet_grid(start~num, scales="free_y")+
   theme_classic()+
   geom_tile()+
   scale_alpha_manual(values = c(0, 1), guide = "none")+
   scale_fill_discrete(name = "Starting tree")+
   scale_x_continuous(name = "Million iterations (1 record per 500 MCMC iterations)")+
   scale_y_discrete(name = "MCMC topology", labels = c("AB|CD", "AC|BD", "BC|AD", "AB|CD", "AC|BD", "BC|AD", "AB|CD", "AC|BD", "BC|AD"))+
   geom_vline(xintercept=0.4, linetype="dashed")+
   theme(plot.tag.position = c(0, 1),plot.tag = element_text(size=14,face = "bold"), panel.border = element_rect(fill="NA"),
         panel.spacing.y = unit(2,"pt"),panel.spacing.x = unit(50,"pt"),
         legend.margin = margin(0,0,0,0,"pt"),legend.box.margin = margin(0,0,0,0,"pt"), strip.text = element_blank()))
ggsave(file="../../bpp_mcmc_trace.pdf", width = 12, height = 8)
