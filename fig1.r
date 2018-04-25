############################################################
# For owl monkey project, last revised 04.18
# Plots and statistical tests for Figure 1 and Figure S2
# in owl monkey paper.
#
# Gregg Thomas
############################################################

library(ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
cat("----------\n")

figs2_opt = FALSE
# CHANGE TO TRUE TO GENERATE FIG. S2

####################
## Estimating model parameters using Kong (2012) data
human_dm0 = 34; human_df = 31; human_num_f = 14.2; human_haploid = 2630000000; human_diploid = human_haploid * 2;
# Parameters observed from Drost & Lee (1995) and Kong (2012)

human_mu_c = (human_num_f / human_df)
cat("mu_c estimated from Kong params: ", human_mu_c, "\n")
human_mu_c = human_mu_c / human_haploid
cat("mu_c per site estimated from Kong params: ", human_mu_c, "\n")
# Estimating mu_c

# Calculate spermatogenic proportion and expected number of spermatogenic divisions
# per year using human (Kong) observations.
human_sc = 16
human_dy1 = ((2.01 / human_haploid) / human_mu_c)
human_sp = human_dy1 / (365/16)
cat("Human expected # of spermatogenic cycles/year: ", human_dy1, "\n")
cat("Human spermatogenic cell proportion given a cycle length of 16 days: ", human_sp, "\n")
cat("----------\n")
# Estimating the expected number of spermatogenic cycles per year

# 2.01 is the slope of the Kong line -- mutations per year in males after puberty.
# We can calculate the expected number of spermatogenic cycles per year (human_dy1) based on mu_c.
# Then, we can estimate the proportion of cells undergoing spermatogenesis at any given time (human_sp)
####################

####################
## Owl monkey prediction (purple dashed) using Kong parameters
om_puberty = 1; om_sc = 10.2;
om_dy1 = (365/om_sc) * human_sp
cat("Owl monkey expected # of spermatogenic cycles/year: ", om_dy1, "\n")
# Estimate the expected number of spermatogenic cycles per year in owl monkeys using the human proportion
# of cells actively dividing (human_sp).

age_range = seq(om_puberty, 15, by=0.1)
rep_long = age_range - om_puberty
# The owl monkey age of puberty set at 1 year and a range of paternal ages

om_mu_g = ((human_df * human_mu_c) + ((human_dm0 * human_mu_c) + (om_dy1 * rep_long * human_mu_c))) / 2
# Calculating the mutation rate (Equation 9 in paper)

om_pred = data.frame(mu_g, age_range)
names(om_pred) = c("Mutation.rate", "Parental.age")
om_pred_reg = lm(mu_g ~ age_range)
cat("----------\n")
####################

####################
## Owl monkey observed data (purple solid)
in_data = read.csv("owl-monkey-filters.csv", header=TRUE)

if(figs2_opt) {
  cur_filter = subset(in_data, Min.DP==20 & Max.DP==60 & Min.AB==0.3 & Max.AB==0.7)
}else{
  cur_filter = subset(in_data, Min.DP==20 & Max.DP==60 & Min.AB==0.4 & Max.AB==0.6)
}
om_obs = cur_filter$MVs.corrected.for.FN
om_gt = cur_filter$Paternal.GT
om_callable = (cur_filter$Sites - cur_filter$Uncallable.SNP.sites) * 2
om_mu = om_obs / om_callable
cur_filter$CpG.rate = cur_filter$CpG.MVs / om_callable
om_obs_reg = lm(om_mu ~ om_gt)
# Read data and define current filter

cat("PREDICTED SLOPE (d_ym1)\n")
cat(coef(om_pred_reg)[2], " mutations per site per year\n")
cat(coef(om_pred_reg)[2] * human_haploid * 2, " mutations per year\n")
cat("OBSERVED SLOPE (dy_m1)\n")
cat(coef(om_obs_reg)[2], " mutations per site per year\n")
cat(coef(om_obs_reg)[2] * (mean(om_callable)), " mutations per year\n")
cat("----------\n")
# Predicted vs observed slopes (d_ym1)
####################

####################
## The plot for Fig. 1
fig1 = ggplot(cur_filter, aes(x=Paternal.GT, y=Mutation.rate.corrected.for.FN)) + 
  geom_smooth(method='glm', color='#b66dff', fill="#f2e6ff", fullrange=T) +
  geom_point(color='#b66dff', size=3) +
  geom_smooth(aes(x=Paternal.GT, y=CpG.rate), method='glm', color="#6db6ff", fill=NA, fullrange=T) +
  geom_point(aes(x=Paternal.GT, y=CpG.rate), color="#6db6ff", size=3) +
  scale_x_continuous(breaks = seq(1, 60, by=2)) + 
  scale_y_continuous(breaks = seq(0, 1.8e-8, by=0.4e-8)) + 
  geom_line(data=pred_df, aes(x=Parental.age, y=Mutation.rate), linetype=2, size=0.75, color="purple") +
  geom_vline(xintercept=1, linetype=3, color="grey", size=0.75) + 
  labs(x="Age of father at conception (y)",y="Mutation rate per site\nper generation") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=20), 
        axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm")
        )
print(fig1)
if(figs2_opt) {
  ggsave(filename="figS2.pdf", fig1, width=6, height=5, units="in")
}else{
  ggsave(filename="fig1.pdf", fig1, width=6, height=5, units="in")
}
cat("--Regression summary for observed points--\n")
print(summary(om_obs_reg))
cat("----------\n")
####################

####################
## The F test for Fig. 1
pred_actual = coef(om_pred_reg)[1] + om_gt * coef(om_pred_reg)[2]
pred_resids = om_mu - pred_actual
#plot(age_range,om_mu_g, type='l',ylim=c(min(om_mu),max(om_mu)))
#points(om_gt, om_mu)
#segments(x0=om_gt,y0=om_mu,x1=om_gt,y1=om_mu-pred_resids)
#abline(om_obs_reg, lty=2)
#segments(x0=om_gt,y0=om_mu,x1=om_gt,y1=om_mu-resid(om_obs_reg),lty=2)
cat("--F test for differences between residuals of the predicted and observed lines--")
print(var.test(resid(om_obs_reg), pred_resids, alternative="two.sided"))
cat("----------\n")
####################
stop()
####################
## Paternal (orange) vs. Maternal (teal) mutations
pat_df = data.frame(cur_filter$Paternal.MVs[cur_filter$Paternal.MVs!=0], cur_filter$Paternal.GT[cur_filter$Paternal.MVs!=0])
names(pat_df) = c("Paternal.MVs", "Paternal.GT")
mat_df = data.frame(cur_filter$Maternal.MVs[cur_filter$Maternal.MVs!=0], cur_filter$Maternal.GT[cur_filter$Maternal.MVs!=0])
names(mat_df) = c("Maternal.MVs", "Maternal.GT")

p = ggplot(pat_df, aes(x=Paternal.GT, y=Paternal.MVs)) +
  geom_smooth(method='glm', color='#db6d00', fill='#ffe6cc', fullrange=T) +
  geom_point(color='#db6d00', size=3) +
  geom_vline(xintercept=1, linetype=3, color="grey", size=0.75) + 
  geom_smooth(data=mat_df, aes(x=Maternal.GT, y=Maternal.MVs), method='glm', color='#009292', fill='#ccffff', fullrange=T) +
  geom_point(data=mat_df, aes(x=Maternal.GT, y=Maternal.MVs), color='#009292', size=3) +
  #scale_y_continuous(limits=c(0,16)) + 
  labs(x="Age of parent at conception (y)",y="Mutation rate per site") +
  theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=20), 
        axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm")
  )
#print(p)
pat_reg = lm(pat_df$Paternal.MVs ~ pat_df$Paternal.GT)
cat("--Regression summary for paternal mutations--\n")
print(summary(pat_reg))
cat("--Regression summary for maternal mutations--\n")
mat_reg = lm(mat_df$Maternal.MVs ~ mat_df$Maternal.GT)
print(summary(mat_reg))
####################
