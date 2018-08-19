############################################################
# For owl monkey project, last revised 04.18
# Plots and statistical tests for Figure 3, Figure S3, and
# Figure S4 in owl monkey paper.
#
# Gregg Thomas
############################################################

library(ggplot2)
library(akima)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
cat("----------\n")

figs2_opt = F
# CHANGE TO TRUE TO GENERATE FIG. S2

figs3_opt = T
# CHANGE TO TRUE TO GENERATE FIG. S3

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
## Human prediction (orange)
human_puberty = 13.4
age_range = seq(human_puberty, 50.4, by=0.1)
rep_long = age_range - human_puberty
# The human age of puberty set at 13.4 years and a range of paternal ages

human_mu_g = ((human_df * human_mu_c) + ((human_dm0 * human_mu_c) + (human_dy1 * rep_long * human_mu_c))) / 2
# Calculating the mutation rate (Equation 9 in paper)

human_pred = data.frame(human_mu_g, age_range)
names(human_pred) = c("Mutation.rate", "Parental.age")
human_pred_reg = lm(human_mu_g ~ age_range)
cat("Human slope (d_ym1): ", coef(human_pred_reg)[2], " -> ", coef(human_pred_reg)[2]*human_diploid, "\n")
cat("----------\n")
####################

####################
## Owl monkey prediction (purple)
if(figs2_opt){
  om_sc = 10.2
}else{
  om_sc = 16
}
om_dy1 = (365/om_sc) * human_sp
cat("Owl monkey expected # of spermatogenic cycles/year: ", om_dy1, "\n")
# Estimate the expected number of spermatogenic cycles per year in owl monkeys using the human proportion
# of cells actively dividing (human_sp).

om_puberty = 1
age_range = seq(om_puberty, 16, by=0.1)
rep_long = age_range - om_puberty
# The owl monkey age of puberty set at 1 year and a range of paternal ages

om_mu_g = ((human_df * human_mu_c) + ((human_dm0 * human_mu_c) + (om_dy1 * rep_long * human_mu_c))) / 2
# Calculating the mutation rate (Equation 9 in paper)

om_pred = data.frame(om_mu_g, age_range)
names(om_pred) = c("Mutation.rate", "Parental.age")
om_pred_reg = lm(om_mu_g ~ age_range)
cat("Owl monkey slope (d_ym1): ", coef(om_pred_reg)[2], " -> ", coef(om_pred_reg)[2]*human_diploid, "\n")
cat("----------\n")
####################

####################
## Chimp prediction (red)
if(figs2_opt){
  chimp_sc = 14
}else{
  chimp_sc = 16
}
chimp_dy1 = (365/chimp_sc) * human_sp
cat("Chimpanzee expected # of spermatogenic cycles/year: ", chimp_dy1, "\n")
# Estimate the expected number of spermatogenic cycles per year in chimpss using the human proportion
# of cells actively dividing (human_sp).

chimp_puberty = 7.5
age_range = seq(chimp_puberty, 35.5, by=0.1)
rep_long = age_range - chimp_puberty
# The owl monkey age of puberty set at 1 year and a range of paternal ages

chimp_mu_g = ((human_df * human_mu_c) + ((human_dm0 * human_mu_c) + (chimp_dy1 * rep_long * human_mu_c))) / 2
# Calculating the mutation rate (Equation 9 in paper)

chimp_pred = data.frame(chimp_mu_g, age_range)
names(chimp_pred) = c("Mutation.rate", "Parental.age")
chimp_pred_reg = lm(chimp_mu_g ~ age_range)
cat("Chimp slope (d_ym1): ", coef(chimp_pred_reg)[2], " -> ", coef(chimp_pred_reg)[2]*human_diploid, "\n")
cat("----------\n")
####################

####################
# Read the study data
study_data = read.csv("data/mutation-studies.csv", header=TRUE)
####################

####################
## The plot for Fig. 3
  if(!figs3_opt){
  fig3 = ggplot(study_data, aes(x=Paternal.age, y=Mutation.rate, color=Species)) +
    geom_line(data=human_pred, aes(x=Parental.age, y=Mutation.rate), color='#db6d00', size=1) +
    geom_segment(aes(x=human_puberty,y=human_mu_g[1]-human_mu_g[1]/10,xend=human_puberty,yend=human_mu_g[1]+human_mu_g[1]/10), linetype=1, color="#db6d00", size=0.5) +
    geom_line(data=om_pred, aes(x=Parental.age, y=Mutation.rate), color='#b66dff', size=1) +
    geom_segment(aes(x=om_puberty,y=om_mu_g[1]-om_mu_g[1]/10,xend=om_puberty,yend=om_mu_g[1]+om_mu_g[1]/10), linetype=1, color="#b66dff", size=0.5) +
    geom_line(data=chimp_pred, aes(x=Parental.age, y=Mutation.rate), color='#920000', size=1) +
    geom_segment(aes(x=chimp_puberty,y=chimp_mu_g[1]-chimp_mu_g[1]/10,xend=chimp_puberty,yend=chimp_mu_g[1]+chimp_mu_g[1]/10), linetype=1, color="#920000", size=0.5) +
    geom_point(size=3) +
    scale_x_continuous(limits=c(0,50), breaks = seq(1, 50, by=5)) + 
    labs(x="Paternal age (y)",y="Mutation rate per site\nper generation") +
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
          legend.position=c(0.13,0.87),
          legend.text=element_text(size=18)
    ) +
    scale_color_manual(breaks=c("Owl monkey", "Chimpanzee", "Human"), values=c("Human"='#db6d00',"Owl monkey"='#b66dff',"Chimpanzee"='#920000'))
  print(fig3)
  if(figs2_opt) {
    ggsave(filename="figS2.pdf", fig3, width=10, height=6, units="in")
  }else{
    ggsave(filename="fig3.pdf", fig3, width=10, height=6, units="in")
  }
  ####################
  
  ####################
  ## The F test for Fig. 3
  human_data = data.frame(study_data$Paternal.age[study_data$Species=="Human"], study_data$Mutation.rate[study_data$Species=="Human"])
  names(human_data) = c("Age","Mu")
  human_data = na.omit(human_data)
  
  human_study_reg = lm(human_data$Mu ~ human_data$Age)
  cat("--Regression summary for human study data--\n")
  print(summary(human_study_reg))
  
  human_pred_reg = lm(human_pred$Mutation.rate ~ human_pred$Parental.age)
  
  pred_actual = coef(human_pred_reg)[1] + human_data$Age * coef(human_pred_reg)[2]
  pred_resids = human_data$Mu - pred_actual
  #plot(human_pred$Parental.age,human_pred$Mutation.rate, type='l',ylim=c(0.8e-8,1.8e-8), xlim=c(13.4,50))
  #points(human_data$Age, human_data$Mu)
  #segments(x0=human_data$Age,y0=human_data$Mu,x1=human_data$Age,y1=human_data$Mu-pred_resids)
  #abline(human_study_reg, lty=2)
  #segments(x0=human_data$Age,y0=human_data$Mu,x1=human_data$Age,y1=human_data$Mu-resid(human_study_reg),lty=2)
  cat("--F test for differences between residuals of the predicted and best fit lines for the human studies--")
  print(var.test(resid(human_study_reg), pred_resids, alternative="two.sided"))
}
####################

####################
## Fig. S3
if(figs3_opt){
  colfunc <- colorRampPalette(c("#009292", "#920000"))
  #colfunc <- colorRampPalette(c("#009292", "#db6d00", "#920000"))
  #colfunc <- colorRampPalette(c("#009292", "#006ddb", "#920000"))
  #colfunc <- colorRampPalette(c("#018571", "#80cdc1", "#a6611a"))
  #colfunc <- colorRampPalette(c("#eecc82", "#926544", "#000000"))
  # A nice(?) color palette
  
  rt = c();nrt = c();
  for(i in seq(5, 51, 1))
  {
    for(j in seq(5, 51, 1))
    {
      rt = c(rt, i)
      nrt = c(nrt, j)
    }
  }
  mu_y = (((human_df * human_mu_c) + (human_dm0 * human_mu_c) + (human_dy1 * rt * human_mu_c)) / (rt + nrt));
  # Calculate mu_y for a range of life history values
  
  am = 20; af = 6; p = 3; rt1 = am-p; nrt1 = af+p
  mu_y_m1 = (((human_df * human_mu_c) + (human_dm0 * human_mu_c) + (human_dy1 * (am-p) * human_mu_c)) / (am + af));
  am = 35; af = 6; p = 3; rt2 = am-p; nrt2 = af+p
  mu_y_m2 = (((human_df * human_mu_c) + (human_dm0 * human_mu_c) + (human_dy1 * (am-p) * human_mu_c)) / (am + af));
  
  am = 20; af = 27; p = 3; rt3 = am-p; nrt3 = af+p
  mu_y_m3 = (((human_df * human_mu_c) + (human_dm0 * human_mu_c) + (human_dy1 * (am-p) * human_mu_c)) / (am + af));
  am = 35; af = 27; p = 3; rt4 = am-p; nrt4 = af+p
  mu_y_m4 = (((human_df * human_mu_c) + (human_dm0 * human_mu_c) + (human_dy1 * (am-p) * human_mu_c)) / (am + af));
  
  rts = c(rt1, rt2, rt3, rt4)
  nrts = c(nrt1, nrt2, nrt3, nrt4)
  muys = c(signif(mu_y_m1, 3), signif(mu_y_m2,3), signif(mu_y_m3,3), signif(mu_y_m4,3))
  # Generate a couple points with varying male age at conception
  
  pdf(NULL)
  dev.control(displaylist="enable")
  # pdf seems to save with some weird grid marks.
  par(mar=c(5,5,4,2))
  pred_data = interp(rt,nrt,mu_y)
  filled.contour(pred_data, xlab='RL (years)', ylab=expression('F + M'[0]*' (years)'), cex.lab=2, col=colfunc(26),
                 key.title = {par(cex.main=2);title(main=expression(mu[y]))},
                 #key.axes = axis(4, seq(2.5e-10, 1.415e-9, by = 0.5e-10)),
                 plot.axes = { contour(pred_data, drawlabels=T, axes=F, frame.plot=F, add=T, col="#999999", lwd=1.5, lty=3, labcex=1);
                   axis(1);
                   axis(2);
                   #segments(rt1, nrt1, rt2, nrt2)
                   #segments(rt3, nrt3, rt4, nrt4)
                   #points(rts, nrts)
                   #points(rts, nrts)
                   #text(rts, nrts, labels=muys, pos=1)
                   
                   points(16.3,43.1,pch=21,bg='#db6d00',col="white",cex=2)
                   # Kong human point
                   points(5.64,7.53,pch=21,bg='#b66dff',col="white",cex=2)
                   # Our owl monkey point
                   points(16.8,33.8,pch=21,bg='#920000',col="white",cex=2)
                   # Venn chimp point
                   #points(5.15,4.19,pch=3,col='#009292',cex=1)
                   # Pfeiffer vervet point
                   
                   #rect(xleft=min(rt_h), ybottom=min(nrt_h), xright=max(rt_h), ytop=max(nrt_h), border="orange")
                   # Human range
                   #rect(xleft=min(rt_ch), ybottom=min(nrt_ch), xright=max(rt_ch), ytop=max(nrt_ch), border="red")
                   # Chimp range
                   #rect(xleft=min(rt_om), ybottom=min(nrt_om), xright=max(rt_om), ytop=max(nrt_om), border="purple")
                   # Owl monkey range
                   #rect(xleft=min(rt_vv), ybottom=min(nrt_vv), xright=max(rt_vv), ytop=max(nrt_vv), border="#009292")
                   # Vervet range
                   
                   #points(rt_p,nrt_p,type='l',col="red")
                   # Varying P line
                   #points(rt_gm,nrt_gm,type='l',col="green")
                   # Varying Gm line
                   #points(rt_gf,nrt_gf,type='l',col="blue")
                   # Varying Gf line
                 })
  plot_obj = recordPlot()
  dev.off()
  # Generate the plot and save it as an object. This allows simultaneous saving and displaying.
  
  png(file="figS3.png", width=8, height=7, units="in", res=400)
  print(plot_obj)
  dev.off()
  # Save plot to png file. pdf saves with weird grid lines.
  
  print(plot_obj)
  # Display plot
}










