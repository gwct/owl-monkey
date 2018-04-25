library(plotrix)
library(ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
##########
# This function calculates mutation rates given a set of parameters
mutCalc <- function(df,d0,dy1,P,mu_c,P_max)
{
  mu_g = c();mu_y = c();alpha = c();gtm = c();rp = c();rt = c();nrt = c();nrp = c()
  i = P
  while(i <= P_max)
  #for(i in P:P_max)
  {
    gm = i
    gf = i
    mu_gF = df * (mu_c)
    mu_g0 = d0 * (mu_c)
    mu_g1 = dy1 * (gm - P) * mu_c
    mu_gM = mu_g0 + mu_g1
    cur_mu_g = (mu_gF + mu_gM) / 2
    cur_mu_y = (mu_gF + mu_gM) / (gm + gf)
    cur_alpha = mu_gM / mu_gF
    mu_g = c(mu_g, cur_mu_g)
    mu_y = c(mu_y, cur_mu_y)
    alpha = c(alpha, cur_alpha)
    gtm = c(gtm, i)
    
    cur_rt = gm - P
    cur_nrt = gf + P
    cur_rp = cur_rt / (cur_nrt)
    cur_nrp = cur_nrt / (gm + gf)
    rt = c(rt, cur_rt)
    nrt = c(nrt, cur_nrt)
    rp = c(rp, cur_rp)
    nrp = c(nrp, cur_nrp)
    i = i + 0.1
  }
  return(list(mu_g, mu_y, alpha, gtm, rt, nrt, rp, nrp))
}
##########

plot_var = "mu_g"
# One of: mu_g, mu_y, alpha, rt, nrt, rp, nrp

####################
## Calculate mu_c using human (Kong) observations.
human_d0 = 34; human_df = 31; human_num_f = 14.2; human_haploid = 2630000000; human_diploid = human_haploid * 2;
human_mu_c = (human_num_f / human_df) / human_haploid
print(paste("mu_c estimated from Kong params: ", human_mu_c))

####################
## Calculate spermatogenic proportion and expected number of spermatogenic divisions
## per year using human (Kong) observations.
human_sc = 16
human_dy1 = (2.01 / human_haploid) / human_mu_c
human_sp = human_dy1 / (365/16)
print(paste("Human expected # of spermatogenic cell divisions/year: ", human_dy1))
print(paste("Human spermatogenic cell proportion given a cycle length of 16 days: ", human_sp))
print("----------")
# 2.01 is the slope of the Kong line -- mutations per year in males after puberty.
# We can calculate the expected number of spermatogenic cycles per year (human_dy1) based on mu_c.
# Then, we can estimate the proportion of cells undergoing spermatogenesis at any given time (human_sp)

# NOTE: average slope of 3 studies: 1.74
# NOTE: average haploid genome size from 3 studies: 2587630000

####################
## Human parameters (orange)
human_p = 13.4; human_max_age = 50.4;
output = mutCalc(human_df,human_d0,human_dy1,human_p,human_mu_c,human_max_age); mu_g_null = unlist(output[1]); mu_y_null = unlist(output[2]); alpha_null = unlist(output[3]); gtm_null = unlist(output[4]); rt_null = unlist(output[5]); nrt_null = unlist(output[6]); rp_null = unlist(output[7]);
####################
## Owl monkey prediction (purple) using Kong parameters
om_p = 1; om_sc = 10.2; om_max_age = 16;
om_dy1 = (365/om_sc) * human_sp
print(paste("Owl monkey expected # of spermatogenic cell divisions/year: ", om_dy1))
t = "Owl monkey mutation rate per generation"
output = mutCalc(human_df,human_d0,om_dy1,om_p,human_mu_c,om_max_age); mu_g_om = unlist(output[1]); mu_y_om = unlist(output[2]); alpha_om = unlist(output[3]); gtm_om = unlist(output[4]); rt_om = unlist(output[5]); nrt_om = unlist(output[6]); rp_om = unlist(output[7])
print("----------")
####################
## Chimp parameters (red solid)
ch_p = 7.5; ch_sc = 14; ch_max_age = 35.5;
ch_dy1 = (365/ch_sc) * human_sp
output = mutCalc(human_df,human_d0,ch_dy1,ch_p,human_mu_c,ch_max_age); mu_g_ch = unlist(output[1]); mu_y_ch = unlist(output[2]); alpha_ch = unlist(output[3]); gtm_ch = unlist(output[4]); rt_ch = unlist(output[5]); nrt_ch = unlist(output[6]); rp_ch = unlist(output[7]); nrp_ch = unlist(output[8]);
####################
## Vervet prediction (aquamarine3 solid)
vv_p = 5; vv_sc = 10.5; vv_max_age = 13;
vv_dy1 = (365/vv_sc) * human_sp
output = mutCalc(human_df,human_d0,vv_dy1,vv_p,human_mu_c,vv_max_age); mu_g_vv = unlist(output[1]); mu_y_vv = unlist(output[2]); alpha_vv = unlist(output[3]); gtm_vv = unlist(output[4]); rt_vv = unlist(output[5]); nrt_vv = unlist(output[6]); rp_vv = unlist(output[7]); nrp_vv = unlist(output[8]);
####################
## Gibbon prediction (blue)
gib_p = 4; gib_sc = 16; gib_max_age = 30;
gib_dy1 = (365/gib_sc) * human_sp
output = mutCalc(human_df,human_d0,gib_dy1,gib_p,human_mu_c,gib_max_age); mu_g_gib = unlist(output[1]); mu_y_gib = unlist(output[2]); alpha_gib = unlist(output[3]); gtm_gib = unlist(output[4]); rt_gib = unlist(output[5]); nrt_gib = unlist(output[6]); rp_gib = unlist(output[7]); nrp_gib = unlist(output[8]);
####################
## Mouse lemur prediction (grey)
ml_p = 0.51; ml_sc = 16; ml_max_age = 18;
ml_dy1 = (365/ml_sc) * human_sp
output = mutCalc(human_df,human_d0,ml_dy1,ml_p,human_mu_c,ml_max_age); mu_g_ml = unlist(output[1]); mu_y_ml = unlist(output[2]); alpha_ml = unlist(output[3]); gtm_ml = unlist(output[4]); rt_ml = unlist(output[5]); nrt_ml = unlist(output[6]); rp_ml = unlist(output[7]); nrp_ml = unlist(output[8]);
####################
## Gorilla prediction (green)
go_p = 7; go_sc = 16; go_max_age = 40;
go_dy1 = (365/go_sc) * human_sp
output = mutCalc(human_df,human_d0,go_dy1,go_p,human_mu_c,go_max_age); mu_g_go = unlist(output[1]); mu_y_go = unlist(output[2]); alpha_go = unlist(output[3]); gtm_go = unlist(output[4]); rt_go = unlist(output[5]); nrt_go = unlist(output[6]); rp_go = unlist(output[7]); nrp_go = unlist(output[8]);
####################
## Orang prediction (yellow)
or_p = 7; or_sc = 16; or_max_age = 50;
or_dy1 = (365/or_sc) * human_sp
output = mutCalc(human_df,human_d0,or_dy1,or_p,human_mu_c,or_max_age); mu_g_or = unlist(output[1]); mu_y_or = unlist(output[2]); alpha_or = unlist(output[3]); gtm_or = unlist(output[4]); rt_or = unlist(output[5]); nrt_or = unlist(output[6]); rp_or = unlist(output[7]); nrp_or = unlist(output[8]);
####################
## plot
if(plot_var=="mu_g"){
  #mu_g_null = mu_g_null/human_diploid; mu_g_om = mu_g_om/human_diploid; mu_g_ch = mu_g_ch/human_diploid; mu_g_vv = mu_g_vv/human_diploid;
  # Convert from rate to rate per site
  
  human_pred = data.frame(gtm_null,mu_g_null)
  names(human_pred) = c("Parental.age", "Mutation.rate")
  ph = lm(mu_g_null ~ gtm_null)
  print(paste("Human slope: ", coef(ph)[2], " -> ", coef(ph)[2]*human_diploid))
  # Human prediction
  
  in_data = read.csv("../owl-monkey-sites.csv", header=TRUE)
  cur_filter = subset(in_data, Min.DP==20 & Max.DP==60 & Min.AB==0.4 & Max.AB==0.6)
  # Read data and define current filter
  
  om_obs = cur_filter$MVs.corrected.for.FN
  om_l = (cur_filter$Sites - cur_filter$Uncallable.SNP.sites) * 2
  om_gt = cur_filter$Paternal.GT
  om_mu = om_obs / om_l
  # Parse the data for a given filter.
  ## OWL MONKEY DATA
  
  om_pred = data.frame(gtm_om,mu_g_om)
  names(om_pred) = c("Parental.age", "Mutation.rate")
  pom = lm(mu_g_om ~ gtm_om)
  print(paste("Owl monkey slope: ", coef(pom)[2], " -> ", coef(pom)[2]*human_diploid))
  # Owl monkey prediction
  
  ch_pred = data.frame(gtm_ch,mu_g_ch)
  names(ch_pred) = c("Parental.age", "Mutation.rate")
  pch = lm(mu_g_ch ~ gtm_ch)
  print(paste("Chimp slope: ", coef(pch)[2], " -> ", coef(pch)[2]*human_diploid))
  # Chimp prediction
  
  vv_pred = data.frame(gtm_vv,mu_g_vv)
  names(vv_pred) = c("Parental.age", "Mutation.rate")
  pvv = lm(mu_g_vv ~ gtm_vv)
  print(paste("Vervet slope: ", coef(pvv)[2], " -> ", coef(pvv)[2]*human_diploid))
  # Vervet prediction
  
  gib_pred = data.frame(gtm_gib,mu_g_gib)
  names(gib_pred) = c("Parental.age", "Mutation.rate")
  pgib = lm(mu_g_gib ~ gtm_gib)
  print(paste("Gibbon slope: ", coef(pgib)[2], " -> ", coef(pgib)[2]*human_diploid))
  # Gibbon prediction
  
  ml_pred = data.frame(gtm_ml,mu_g_ml)
  names(ml_pred) = c("Parental.age", "Mutation.rate")
  pml = lm(mu_g_ml ~ gtm_ml)
  print(paste("Mouse lemur slope: ", coef(pml)[2], " -> ", coef(pml)[2]*human_diploid))
  # Mouse lemur prediction
  
  go_pred = data.frame(gtm_go,mu_g_go)
  names(go_pred) = c("Parental.age", "Mutation.rate")
  pgo = lm(mu_g_go ~ gtm_go)
  print(paste("Gorilla slope: ", coef(pgo)[2], " -> ", coef(pgo)[2]*human_diploid))
  # Gorilla prediction
  
  or_pred = data.frame(gtm_or,mu_g_or)
  names(or_pred) = c("Parental.age", "Mutation.rate")
  por = lm(mu_g_or ~ gtm_or)
  print(paste("Orang slope: ", coef(por)[2], " -> ", coef(por)[2]*human_diploid))
  # Orang prediction
  
  study_data = read.csv("mutation-studies.csv", header=TRUE)
  study_data = subset(study_data, Species !="Gorilla" & Species != "Orang" & Species != "Chimp" & Species!="Mouse lemur" & Species != "Gibbon" & Species != "Vervet")
  #study_data$Species = factor(study_data$Species, levels = c("Owl monkey", "Chimp", "Human", "Vervet", "Mouse lemur", "Gibbon", "Chimp2", "Gorilla", "Orang"))
  #study_data$Species = factor(study_data$Species, levels = c("Owl monkey", "Chimp", "Human"))
  # Read the data from various studies
  
  study_data_no_om = subset(study_data, Species!="Owl monkey")
  cur_filter$Species = "Owl monkey"
  cur_filter$Avg.parental.age = (cur_filter$Paternal.GT + cur_filter$Maternal.GT) / 2
  
  human_data = data.frame(study_data$Avg.parental.age[study_data$Species=="Human"], study_data$Mutation.rate[study_data$Species=="Human"])
  names(human_data) = c("Age","Mu")
  human_data = na.omit(human_data)
  human_study_reg = lm(human_data$Mu ~ human_data$Age)
  human_pred_reg = lm(human_pred$Mutation.rate ~ human_pred$Parental.age)
  
  p = ggplot(study_data, aes(x=Paternal.age, y=Mutation.rate, color=Species)) +
    geom_line(data=human_pred, aes(x=Parental.age, y=Mutation.rate), color='#db6d00', size=1) +
    geom_segment(aes(x=human_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=human_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="#db6d00", size=0.5) +
    geom_line(data=om_pred, aes(x=Parental.age, y=Mutation.rate), color='#b66dff', size=1) +
    geom_segment(aes(x=om_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=om_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="#b66dff", size=0.5) +
    geom_line(data=ch_pred, aes(x=Parental.age, y=Mutation.rate), color='#920000', size=1) +
    geom_segment(aes(x=ch_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=ch_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="#920000", size=0.5) +
    #geom_line(data=vv_pred, aes(x=Parental.age, y=Mutation.rate), color='#009292', size=1) +
    #geom_segment(aes(x=vv_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=vv_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="#009292", size=0.5) +
    #geom_line(data=gib_pred, aes(x=Parental.age, y=Mutation.rate), color='blue', size=1) +
    #geom_segment(aes(x=gib_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=gib_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="grey", size=0.5) +
    #geom_line(data=ml_pred, aes(x=Parental.age, y=Mutation.rate), color='grey', size=1) +
    #geom_segment(aes(x=ml_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=ml_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="grey", size=0.5) +
    #geom_line(data=go_pred, aes(x=Parental.age, y=Mutation.rate), color='green', size=1) +
    #geom_segment(aes(x=go_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=go_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="grey", size=0.5) +
    #geom_line(data=or_pred, aes(x=Parental.age, y=Mutation.rate), color='black', size=1) +
    #geom_segment(aes(x=or_p,y=mu_g_null[1]-mu_g_null[1]/10,xend=or_p,yend=mu_g_null[1]+mu_g_null[1]/10), linetype=1, color="grey", size=0.5) +
    #geom_hline(yintercept=1.022e-8) +
    #geom_vline(xintercept=13) + 
    #geom_vline(xintercept=25.32) + 
    #geom_hline(yintercept=5.656814e-09, color="#333333", linetype=3) + 
    geom_point(size=3) +
    #geom_smooth(data=human_data, aes(x=Age, y=Mu), color="#db6d00", se=F, method='glm', linetype=2) + 
    #geom_point(data=cur_filter, aes(x=Avg.parental.age, y=Mutation.rate.corrected.for.FN, color=Species)) +
    #geom_point(aes(x=29.7, y=1.2e-8), colour="blue") +
    scale_x_continuous(limits=c(0,50), breaks = seq(1, 50, by=5)) + 
    #scale_y_continuous(breaks = seq(0, 1.8e-8, by=0.4e-8)) + 
    #geom_vline(xintercept=1, linetype=3, color="grey", size=0.75) + 
    labs(x="Paternal age (y)",y="Mutation rate per site") +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=20), 
          axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),color="black"), 
          axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0),color="black"),
          axis.line=element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959",size = 1),
          axis.ticks.length=unit(0.2,"cm"),
          legend.title=element_blank(),
          legend.position=c(0.1,0.87),
          legend.text=element_text(size=18)
    ) +
    scale_color_manual(breaks=c("Owl monkey", "Vervet", "Chimpanzee", "Human"), values=c("Human"='#db6d00',"Owl monkey"='#b66dff',"Chimpanzee"='#920000',"Vervet"='#009292'))#,"Gibbon"='blue',"Mouse lemur"='grey',"Gorilla"='green',"Orang"='black'))
  print(p)
  # Call the plot fjuijuiiiizunction
  stop()
  pred_actual = coef(human_pred_reg)[1] + human_data$Age * coef(human_pred_reg)[2]
  pred_resids = human_data$Mu - pred_actual
  plot(human_pred$Parental.age,human_pred$Mutation.rate, type='l',ylim=c(0.8e-8,1.8e-8), xlim=c(13.4,50))
  points(human_data$Age, human_data$Mu)
  segments(x0=human_data$Age,y0=human_data$Mu,x1=human_data$Age,y1=human_data$Mu-pred_resids)
  abline(human_study_reg, lty=2)
  segments(x0=human_data$Age,y0=human_data$Mu,x1=human_data$Age,y1=human_data$Mu-resid(human_study_reg),lty=2)
  print(var.test(resid(human_study_reg), pred_resids, alternative="two.sided"))
  # This is the F-test for differences between the predicted and observed human lines in Fig.3
  
}
if(plot_var=="mu_y"){
  mu_y_null = mu_y_null/human_diploid; mu_y_om = mu_y_om/human_diploid; mu_y_ch = mu_y_ch/human_diploid; mu_y_vv = mu_y_vv/human_diploid;
  # Convert from rate to rate per site
  
  human_pred = data.frame(gtm_null,mu_y_null)
  names(human_pred) = c("Generation.time", "Substitution.rate")
  om_pred = data.frame(gtm_om,mu_y_om)
  names(om_pred) = c("Generation.time", "Substitution.rate")
  ch_pred = data.frame(gtm_ch,mu_y_ch)
  names(ch_pred) = c("Generation.time", "Substitution.rate")
  vv_pred = data.frame(gtm_vv,mu_y_vv)
  names(vv_pred) = c("Generation.time", "Substitution.rate")
  p = ggplot() +
    geom_line(data=human_pred, aes(x=Generation.time, y=Substitution.rate), color='#db6d00', size=1) +
    geom_line(data=om_pred, aes(x=Generation.time, y=Substitution.rate), color='#b66dff', size=1) +
    geom_line(data=ch_pred, aes(x=Generation.time, y=Substitution.rate), color='#920000', size=1) +
    geom_line(data=vv_pred, aes(x=Generation.time, y=Substitution.rate), color='#009292', size=1) +
    scale_x_continuous(limits=c(1,60), breaks = seq(1, 60, by=5)) + 
    #scale_y_continuous(breaks = seq(0, 1.8e-8, by=0.4e-8)) + 
    #geom_vline(xintercept=1, linetype=3, color="grey", size=0.75) + 
    labs(x="Average generation time (y)",y="Substitution rate per site") +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=20), 
          axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),color="black"), 
          axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0),color="black"),
          axis.line=element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959",size = 1),
          axis.ticks.length=unit(0.2,"cm"),
          legend.title=element_blank(),
          legend.position=c(0.1,0.87),
          legend.text=element_text(size=14)
    )
  print(p)
}
if(plot_var=="alpha"){
  human_pred = data.frame(gtm_null,alpha_null)
  names(human_pred) = c("Generation.time", "Alpha")
  om_pred = data.frame(gtm_om,alpha_om)
  names(om_pred) = c("Generation.time", "Alpha")
  ch_pred = data.frame(gtm_ch,alpha_ch)
  names(ch_pred) = c("Generation.time", "Alpha")
  vv_pred = data.frame(gtm_vv,alpha_vv)
  names(vv_pred) = c("Generation.time", "Alpha")
  
  p = ggplot() +
    geom_line(data=human_pred, aes(x=Generation.time, y=Alpha), color='#db6d00', size=1) +
    geom_line(data=om_pred, aes(x=Generation.time, y=Alpha), color='#b66dff', size=1) +
    geom_line(data=ch_pred, aes(x=Generation.time, y=Alpha), color='#920000', size=1) +
    geom_line(data=vv_pred, aes(x=Generation.time, y=Alpha), color='#009292', size=1) +
    scale_x_continuous(limits=c(1,60), breaks = seq(1, 60, by=5)) + 
    #scale_y_continuous(breaks = seq(0, 1.8e-8, by=0.4e-8)) + 
    #geom_vline(xintercept=1, linetype=3, color="grey", size=0.75) + 
    labs(x="Average generation time (y)",y="M:F mutation ratio") +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=20), 
          axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),color="black"), 
          axis.title.x=element_text(margin=margin(t=20,r=0,b=0,l=0),color="black"),
          axis.line=element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959",size = 1),
          axis.ticks.length=unit(0.2,"cm"),
          legend.title=element_blank(),
          legend.position=c(0.1,0.87),
          legend.text=element_text(size=14)
    )
  print(p)
}
if(plot_var=="rt"){
  t = "Reproductive period predictions"
  yl = "Reproductive period"
  plotMu(gtm_null,rt_null,gtm_ch,mu_g_ch,gtm_om,rt_om,gtm_mac,rt_mac,gtm_vv,rt_vv,gtm_ml,rt_ml,gtm_mo,rt_mo,yl,t)
}
if(plot_var=="nrt"){
  t = "Non-reproductive period predictions"
  yl = "Non-reproductive period"
  plotMu(gtm_null,nrt_null,gtm_ch,mu_g_ch,gtm_om,nrt_om,gtm_mac,nrt_mac,gtm_vv,nrt_vv,gtm_ml,nrt_ml,gtm_mo,nrt_mo,yl,t)
}
if(plot_var=="rp"){
  t = "Reproductive proportion predictions"
  yl = "Reproductive proportion"
  plotMu(gtm_null,rp_null,gtm_om,rp_om,gtm_mac,rp_mac,gtm_vv,rp_vv,gtm_ml,rp_ml,gtm_mo,rp_mo,yl,t)
}
if(plot_var=="nrp"){
  t = "Non-reproductive proportion predictions"
  yl = "Non-reproductive proportion"
  plotMu(gtm_null,nrp_null,gtm_ch,mu_g_ch,gtm_om,nrp_om,gtm_mac,nrp_mac,gtm_vv,nrp_vv,gtm_ml,nrp_ml,gtm_mo,nrp_mo,yl,t)
}






#abline(v=P, lty=2, col="grey")
# Plot the alt line.
#rect(2,mu_g_om[which(gtm_om==2)],14,mu_g_om[which(gtm_om==14)],col=rgb(0.76,0.76,0.84,alpha=0.2), border=NA)











