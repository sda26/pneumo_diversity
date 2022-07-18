#Haplotype Analysis

require(MASS)
require(tidyverse)
require(lmerTest)
require(gridExtra)
require(robustlmm)
require(latex2exp)
require(ggthemes)
require(ggridges)
require(reshape2)
require(fitdistrplus)

theme_set(theme_grey())

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

Diversity = function(Counts,Type) {
  p_i = Counts / sum(Counts)
  if (Type == "Shannon") {
    
    return(-sum(p_i * log(p_i)))
    
  }
  
  if (Type == "Hill_N1") {
    return(exp(-sum(p_i * log(p_i))))
    
  }
  
  if (Type == "Simpson") {
    
    return(sum(p_i * p_i))
    
  }
  
  if (Type == "Hill_N2") {
    
    return(1 / sum(p_i * p_i))
    
  }
  
  if (Type == "Shannon") {
    
    return(-sum(p_i * log(p_i)))
    
  }
  
  if (Type == "Unique") {
    
    return(length(Counts))
    
  }
  
}

Proj.Home = "/Users/aggars05/Desktop/Analysis_Tools"
setwd(Proj.Home)

Base.Plot = ggplot() +
  theme(text=element_text(size=14,face="bold"),
        plot.title = element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=0,vjust=0,size=rel(1.5),color = "black"),
        axis.text.y=element_text(size=rel(1.5),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5, color = "black"),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.title.y = element_text(size=rel(1.2),color = "black"),
        axis.title.x = element_text(size=rel(1.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Read data ----
Haplotypes = read.csv(file = file.path(Proj.Home,"output","SDA281","Spn_Haplotypes.csv"), header = TRUE)

# Summarize Haplotype data ----

Hapl.Sum = Haplotypes %>%
  summarise(Simpson = sum((Count / sum(Count)) * (Count / sum(Count))),
            Shannon = -sum((Count / sum(Count)) * log(Count / sum(Count))),
            Menhinick = length(Count) / sqrt(sum(Count)),
            Hill.N1 = exp(-sum((Count / sum(Count)) * log(Count / sum(Count)))),
            Hills.N2 = 1 / sum((Count / sum(Count)) * (Count / sum(Count))),
            Unique.Genotypes = Count %>% length,
            Singlets = which(Count == 1) %>% length,
            Dominant.Number = max(Count),
            Dominant.Proportion = max(Count)/sum(Count),
            Dominant.Percent = Dominant.Proportion * 100)
Hapl.Sum

write.csv(Hapl.Sum,
          file = "output/SDA281/Haplotypes_Summary.csv",
          row.names = FALSE)

# Fit Poisson (expected for even barcode frequencies) and negative binomial (overdispersed Poisson, deviation from evenness) distributions

mP = fitdist(Haplotypes$Count,distr = "pois")
summary(mP)

fitP = cbind(reads = 1:max(Haplotypes$Count),
             expect = dpois(1:max(Haplotypes$Count), lambda = mP$estimate[1])) %>%
  data.frame

fitP$expect = fitP$expect / max(fitP$expect) * max(as.numeric(table(Haplotypes$Count)))

mNB = fitdist(Haplotypes$Count,distr = "nbinom")
summary(mNB)

fitNB = cbind(reads = 1:max(Haplotypes$Count),
              expect = dnbinom(1:max(Haplotypes$Count), size = mNB$estimate[1], mu = mNB$estimate[2])) %>%
  data.frame

fitNB$expect = fitNB$expect / max(fitNB$expect) * max(as.numeric(table(Haplotypes$Count)))

Base.Plot +
  geom_histogram(data = Haplotypes,
                 aes(x = Count),
                 color = "black",
                 fill = "darkgray",
                 binwidth = 1) +
  geom_line(data = fitP,
            aes(x = reads,
                y = expect),
            color = "black",
            lty = 2,
            lwd = 0.8) +
  geom_line(data = fitNB,
            aes(x = reads,
                y = expect),
            color = "red",
            lwd = 1) +
  labs(x = "Reads/barcode",
       y = "Occurrence count")

ggsave(file = "output/SDA281/Barcode_Histogram_Linear.pdf",
       width = 6,
       height = 6,
       unit = "in")
