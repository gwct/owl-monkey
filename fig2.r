this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
####################
in_data = read.csv("spectra.csv", header=TRUE)
om_spec = data.frame(in_data$Mutation, in_data$Owl.monkey.CpG.proportion, in_data$Owl.monkey.proportion)
names(om_spec) = c("Type", "CpG", "Non-CpG")
om_spec_m = melt(om_spec)

om_p = ggplot(om_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Non-CpG"="#006ddb","CpG"="#db6d00")) +
  labs(x="",y="Proportion of mutations") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.title=element_blank(),
        legend.position=c(0.1,0.87)
  )

human_spec = data.frame(in_data$Mutation, in_data$Human.average.proportion)
names(human_spec) = c("Type", "Human average")
human_spec_m = melt(human_spec)

h_p = ggplot(human_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Human average" = "#920000")) +
  labs(x="",y="Proportion of mutations") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="none"
  )

ch_spec = data.frame(in_data$Mutation, in_data$Chimp.proportion)
names(ch_spec) = c("Type", "Proportion")
ch_spec_m = melt(ch_spec)

c_p = ggplot(ch_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Proportion" = "#920000")) +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="none"
  )

g1 = grid.arrange(h_p,c_p,nrow=1)
g = plot_grid(om_p,g1, nrow=2, rel_heights=c(2, 1))
print(g)
#ggsave(filename="fig2c.pdf", g, width=8, height=6, units="in")

print("--OWL MONKEY Chi-squared test using human proportions as expectations--")
print("Using all mutation categories:")
all_chi = chisq.test(in_data$Owl.monkey+in_data$Owl.monkey.CpG, p=in_data$Human.average.proportion)
print(all_chi)
print("Leaving out A>T:")
human_noAT = in_data$Human.average.proportion[in_data$Mutation!="A>T"]
human_noAT = human_noAT / sum(human_noAT)
om_noAT = in_data$Owl.monkey[in_data$Mutation!="A>T"] + in_data$Owl.monkey.CpG[in_data$Mutation!="A>T"]
noAT_chi = chisq.test(om_noAT, p=human_noAT)
print(noAT_chi)

print("--CHIMP Chi-squared test using human proportions as expectations--")
print("Using all mutation categories:")
all_chi = chisq.test(in_data$Chimp, p=in_data$Human.average.proportion)
print(all_chi)
print("Leaving out A>T:")
ch_noAT = in_data$Chimp[in_data$Mutation!="A>T"]
noAT_chi = chisq.test(ch_noAT, p=human_noAT)
print(noAT_chi)

stop()
context = read.csv("context.csv", header=T)

p = ggplot(context, aes(Base2, Base1)) +
  geom_tile(aes(fill = Proportion), color = "white") +
  facet_grid(. ~ Mutation) +
  scale_fill_continuous(limits=c(0,0.2), breaks=seq(0,0.2,by=0.1), low="#490092", high="#24ff24") +
  labs(x="",y="") +
  #scale_x_discrete(position = "top") +
  theme_classic() +
  theme(
    plot.margin=unit(c(0,0.2,1.2,0.425),"in"),
    axis.text.x=element_text(size=16),
    axis.text.y=element_text(size=16),
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    legend.position="bottom",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.title=element_blank()
  )

g = plot_grid(om_p,p,h_p, nrow=3)#, rel_heights=c(2, 1, 2))
#ggsave(file="fig2d.pdf", g, width=8, height=6, units="in")
