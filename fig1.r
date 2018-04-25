library(ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
##########
# This function calculates mutation rates given a set of parameters
mutCalc <- function(df,d0,dy1,P,mu_c)
{
  P_max = P + 100
  mu_g = c(); gtm = c();
  
  i = P
  while(i <= P_max)
  {
    gm = i
    gf = i
    mu_gF = (df * mu_c)
    mu_g0 = d0 * (mu_c)
    mu_g1 = dy1 * (gm - P) * (mu_c)
    mu_gM = mu_g0 + mu_g1
    cur_mu_g = (mu_gF + mu_gM) / 2
    mu_g = c(mu_g, cur_mu_g)
    gtm = c(gtm, i)
    
    i = i + 0.1
  }
  return(list(mu_g, gtm))
}

####################
## Calculate mu_c using human (Kong) observations.
human_d0 = 34; human_df = 31; human_num_f = 14.2; human_haploid = 2630000000; human_diploid = human_haploid * 2;
human_mu_c = (human_num_f / human_df) / human_haploid
print(paste("mu_c estimated from Kong params: ", human_mu_c))
####################
## Calculate spermatogenic proportion and expected number of spermatogenic divisions
## per year using human (Kong) observations.
human_sc = 16
human_dy1 = ((2.01 / human_haploid) / human_mu_c)
human_sp = human_dy1 / (365/16)
print(paste("Human expected # of spermatogenic cell divisions/year: ", human_dy1))
print(paste("Human spermatogenic cell proportion given a cycle length of 16 days: ", human_sp))
print("----------")
# 2.01 is the slope of the Kong line -- mutations per year in males after puberty.
# We can calculate the expected number of spermatogenic cycles per year (human_dy1) based on mu_c.
# Then, we can estimate the proportion of cells undergoing spermatogenesis at any given time (human_sp)

####################
## Owl monkey prediction (purple dotted) using Kong parameters
om_p = 1; om_sc = 10.2;
om_dy1 = (365/om_sc) * human_sp
print(paste("Owl monkey expected # of spermatogenic cell divisions/year: ", om_dy1))



output = mutCalc(human_df,human_d0,om_dy1,om_p,human_mu_c)
mu_g_alt = unlist(output[1]); gtm_alt = unlist(output[2]);
print("----------")
####################
## Owl monkey data (purple solid)
in_data = read.csv("owl-monkey-sites.csv", header=TRUE)
cur_filter = subset(in_data, Min.DP==20 & Max.DP==60 & Min.AB==0.3 & Max.AB==0.7)
# Read data and define current filter

om_obs = cur_filter$MVs.corrected.for.FN
om_l = (cur_filter$Sites - cur_filter$Uncallable.SNP.sites) * 2
om_gt = cur_filter$Paternal.GT
om_mu = om_obs / om_l
#cur_filter$CpG.rate = cur_filter$CpG.MVs / ((1 - cur_filter$Min.FN) * om_l)
cur_filter$CpG.rate = cur_filter$CpG.MVs / om_l
# Parse the data for a given filter.

om_pred_reg = lm(mu_g_alt ~ gtm_alt)
om_obs_reg = lm(om_mu ~ om_gt)
print("PREDICTED SLOPE")
print(paste("", (coef(om_pred_reg)[2])))
print(paste("", coef(om_pred_reg)[2] * human_haploid * 2))
print("OBSERVED SLOPE")
print(paste("", coef(om_obs_reg)[2]))
print(paste("", coef(om_obs_reg)[2] * (mean(om_l))))
# Some predictions

####################
## The plot
gtm_alt = gtm_alt[which(gtm_alt<=15)]
mu_g_alt = mu_g_alt[which(gtm_alt<=15)]

pred_df = data.frame(gtm_alt,mu_g_alt)
names(pred_df) = c("Parental.age", "Mutation.rate")

cur_mu = cur_filter$Mutation.rate.corrected.for.FN
cur_gt = cur_filter$Paternal.GT

obs_reg = lm(cur_filter$Mutation.rate.corrected.for.FN ~ cur_filter$Paternal.GT)
pred_reg = lm(pred_df$Mutation.rate ~ pred_df$Parental.age)

p = ggplot(cur_filter, aes(x=Paternal.GT, y=Mutation.rate.corrected.for.FN)) + 
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
print(p)
print(summary(obs_reg))
#"#cce6ff"
# Call the plot fjuijuiiiizunction
## This is the main figure 1
stop()
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
print(p)
pat_reg = lm(pat_df$Paternal.MVs ~ pat_df$Paternal.GT)
print(summary(pat_reg))
mat_reg = lm(mat_df$Maternal.MVs ~ mat_df$Maternal.GT)
print(summary(mat_reg))
# This shows Paternal and Maternal MVs
#stop()
pred_actual = coef(pred_reg)[1] + cur_gt * coef(pred_reg)[2]
pred_resids = cur_mu - pred_actual
plot(gtm_alt,mu_g_alt, type='l',ylim=c(min(cur_mu),max(cur_mu)))
points(cur_gt, cur_mu)
segments(x0=cur_gt,y0=cur_mu,x1=cur_gt,y1=cur_mu-pred_resids)
abline(obs_reg, lty=2)
segments(x0=cur_gt,y0=cur_mu,x1=cur_gt,y1=cur_mu-resid(obs_reg),lty=2)
print(var.test(resid(obs_reg), pred_resids, alternative="two.sided"))
# This is the F-test for differences between the two lines in Fig.1 