############################################################
# For owl monkey project, last revised 04.18
# Plots and statistical tests for Figure 2 in owl monkey
# paper.
#
# Gregg Thomas
############################################################

library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

####################
## Read the input data
in_data = read.csv("data/mutational-spectra.csv", header=TRUE)
####################

####################
## Generate the owl monkey plot
om_spec = data.frame(in_data$Mutation, in_data$Owl.monkey.CpG.proportion, in_data$Owl.monkey.proportion)
names(om_spec) = c("Type", "CpG", "Non-CpG")
om_spec_m = melt(om_spec)

om_p = ggplot(om_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Non-CpG"="#006ddb","CpG"="#db6d00")) +
  ggtitle("Owl monkey") + 
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
####################

####################
## Generate the human plot
human_spec = data.frame(in_data$Mutation, in_data$Human.average.proportion)
names(human_spec) = c("Type", "Human average")
human_spec_m = melt(human_spec)

h_p = ggplot(human_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Human average" = "#920000")) +
  ggtitle("Human") + 
  labs(x="",y="Proportion of mutations") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="none"
  )
####################

####################
## Generate the chimp plot
ch_spec = data.frame(in_data$Mutation, in_data$Chimp.proportion)
names(ch_spec) = c("Type", "Proportion")
ch_spec_m = melt(ch_spec)

c_p = ggplot(ch_spec_m, aes(Type, value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual("legend", values = c("Proportion" = "#920000")) +
  ggtitle("Chimpanzee") + 
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="none"
  )
####################

####################
## Arrange the three plots and display/save
g1 = arrangeGrob(h_p,c_p,nrow=1)
fig2 = plot_grid(om_p,g1, nrow=2, rel_heights=c(2, 1))
print(fig2)
ggsave(filename="fig2.pdf", fig2, width=8, height=6, units="in")
####################

####################
## Chi-square tests to check for significant differences between human and owl monkey spectra
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
####################

