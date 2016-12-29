
# ----------------------------------------------------------------------------------------------------------------------------------
#                                                     confounder estimation
#                                         by: JC Barnes, University of Cincinnati
#                                                   version date: 7/2015
# ----------------------------------------------------------------------------------------------------------------------------------

# this simulation draws on the beta distribution, see the accompanying betaDistribution.R file for more detail
# the user is encouraged to explore the beta distribution properties before proceeding
# this will ensure that the estimates gleaned below are meaningful

# clear the workspace
remove(list=ls())

# ----------------------------------------------------------------------------------------------
# write the function to estimate the level of confounding that exists
# ----------------------------------------------------------------------------------------------

confound<-function(a_x,b_x,a_y,b_y,a_p,b_p,rg,num_calcs) {
  
  # need to adapt the h2_trait1 to adjust a_trait1 and b_trait1
    
  # set the seed for reproducibility
  set.seed(2217)
  
  # generate "population" values for the distribution of h2s and rps
  popN<-100000
  h2_x_pop<-   rbeta(popN,shape1=a_x,shape2=b_x)
  h2_y_pop<-   rbeta(popN,shape1=a_y,shape2=b_y)
  r_p_pop<-    rbeta(popN,shape1=a_p,shape2=b_p)
  
  # number of times to run the calculation below
  n<-num_calcs
  
  # generate object to hold the estimates
  estimates<-NA
  estimates<-rep(NA,length(n))
  
  # run a loop to calculate % of rp due to share genetic overlap 
  for (i in 1:n) {
    h2_x<-          sample(h2_x_pop,1,replace=F)
    h2_y<-          sample(h2_y_pop,1,replace=F)
    rp<-            sample(r_p_pop,1,replace=F)
    # bound the estimates at 0
    estimates[i]<-(sqrt(h2_x)*rg*sqrt(h2_y))/rp
    ifelse(estimates[i]<0,estimates[i]<-0,estimates[i]<-estimates[i])
    # bound the estimates at 1
    ifelse(estimates[i]>1,estimates[i]<-1,estimates[i]<-estimates[i])
  } # end of the loop
  
  # -----------------
  # output 
  # -----------------
  
  # function to estimate the mode 
  # (from Rasmus Baath: http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode),
  # to be used below.
  estimate_mode <- function(x) {
    d <- density(x)
    d$x[which.max(d$y)]
  }
  
  # ----
  # beta distributions used for the estimates
  # ----
  
  #par(mfrow=c(3,1))
  
  # h2 for x
  #hist(h2_x_pop,main=expression(italic(h^2)~"of menarche timing"),ylab="",xlab="",xlim=c(0,1),axes=F);axis(1)
  # text for the mean, mode, and standard deviation
  #text(ifelse(min(h2_x_pop)<=.2,.9,.1),popN*.25,"empirical estimates:")
  #text(ifelse(min(h2_x_pop)<=.2,.9,.1),popN*.20,paste("mean=",round(mean(h2_x_pop),3)))
  #text(ifelse(min(h2_x_pop)<=.2,.9,.1),popN*.15,paste("mode=",round(estimate_mode(x=h2_x_pop),3))) 
  #text(ifelse(min(h2_x_pop)<=.2,.9,.1),popN*.10,paste("sd=",round(sd(h2_x_pop),3))) 
  
  # h2 for y
  #hist(h2_y_pop,main=expression(italic(h^2)~"of father absence"),ylab="",xlab="",xlim=c(0,1),axes=F);axis(1)
  # text for the mean, mode, and standard deviation
  #text(ifelse(min(h2_y_pop)<=.2,.9,.1),popN*.25,"empirical estimates:")
  #text(ifelse(min(h2_y_pop)<=.2,.9,.1),popN*.20,paste("mean=",round(mean(h2_y_pop),3)))
  #text(ifelse(min(h2_y_pop)<=.2,.9,.1),popN*.15,paste("mode=",round(estimate_mode(x=h2_y_pop),3))) 
  #text(ifelse(min(h2_y_pop)<=.2,.9,.1),popN*.10,paste("sd=",round(sd(h2_y_pop),3))) 
  
  # rp
  #hist(r_p_pop,main=expression("phenotypic correlation"~(italic(r[p]))),ylab="",xlab="",xlim=c(0,1),axes=F);axis(1)
  # text for the mean, mode, and standard deviation
  #text(ifelse(min(r_p_pop)<=.2,.9,.1),popN*.25,"empirical estimates:")
  #text(ifelse(min(r_p_pop)<=.2,.9,.1),popN*.20,paste("mean=",round(mean(r_p_pop),3)))
  #text(ifelse(min(r_p_pop)<=.2,.9,.1),popN*.15,paste("mode=",round(estimate_mode(x=r_p_pop),3))) 
  #text(ifelse(min(r_p_pop)<=.2,.9,.1),popN*.10,paste("sd=",round(sd(r_p_pop),3))) 
  
  # numerical output of the estimates
  print(paste("rg =",rg))
  print(paste("mean =",round(mean(estimates,na.rm=T),4)))
  print(paste("mode =",round(estimate_mode(estimates),4)))
  print(paste("variance =",round(var(estimates,na.rm=T),4)))
  print(paste("95% credible interval =",paste(round(quantile(estimates,probs=.025),4),round(quantile(estimates,probs=.975),4),sep=",")))
  print("----------------")
  
  # grapchical output of the estimates for h2cov
  #layout(matrix(c(1,2),2,1))
  #par(mfrow=c(1,1))
  rgi<-round(rg[1],2)
  hist(estimates,main=substitute(paste(italic(r[g])~"= ",rgi)),,ylab="",xlab="",axes=F,xlim=c(-.05,1.05));axis(1,at=c(0,.5,1));abline(v=mean(estimates),lwd=3,col="blue");abline(v=(quantile(estimates,probs=.025)),lwd=1,col="green");abline(v=(quantile(estimates,probs=.975)),lwd=1,col="green")
    # text for the mean, mode, and standard deviation
    #text(ifelse(min(estimates)<=.2,.9,.1),num_calcs*.25,"empirical estimates:")
    #text(ifelse(min(estimates)<=.2,.9,.1),num_calcs*.20,paste("mean=",round(mean(estimates),3)))
    #text(ifelse(min(estimates)<=.2,.9,.1),num_calcs*.15,paste("mode=",round(estimate_mode(x=estimates),3))) 
    #text(ifelse(min(estimates)<=.2,.9,.1),num_calcs*.10,paste("sd=",round(sd(estimates),3))) 
 
  #plot(ts(estimates),ylab="",xlab="iteration",main="estimates of % confounded",axes=F,ylim=c(0,1.5));axis(1);axis(2)
  #  abline(h=mean(estimates),lwd=3,col="blue")
  #  abline(h=(mean(estimates)+1.96*(sd(estimates))),lwd=1,col="green")
  #  abline(h=(mean(estimates)-1.96*(sd(estimates))),lwd=1,col="green")
}
# ----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------
# run the function that was written above 
# with user-specified values for the beta distributions and a specific value of rg
# ----------------------------------------------------------------------------------------------

setwd("/Users/jcbarnes/Box Sync/manuscripts/_underReview/lifeHistory&EarlyMenarche_--_ehb/data")

# enter rg values to loop over
rg<-seq(.0,.20,length.out=16)

# test a for loop
pdf("figure2.pdf")
par(mfrow=c(4,4))
for(i in rg) { 
confound(a_x=10,b_x=10,       # this will produce distribution centered on 0.50 for h2 of menarche
         a_y=10,b_y=10,       # this will produce distribution centered on 0.50 for h2 of father-absence
         a_p=35,b_p=300,      # this will produce distribution centered on 0.10 for h2 of menarche
         rg=i,
         num_calcs=10000)
}
dev.off()
#















# ----------------------------------------------------------------------------------------------
# test function by hand
# ----------------------------------------------------------------------------------------------
# rg=0.01
#confound(a_x=10,b_x=10,       # this will produce distribution centered on 0.50 for h2 of menarche
#         a_y=10,b_y=10,       # this will produce distribution centered on 0.50 for h2 of father-absence
#         a_p=35,b_p=300,      # this will produce distribution centered on 0.10 for h2 of menarche
#         rg=0.01,
#         num_calcs=10000)
#
# rg=0.2
#confound(a_x=10,b_x=10,       # this will produce distribution centered on 0.50 for h2 of menarche
#         a_y=10,b_y=10,       # this will produce distribution centered on 0.50 for h2 of father-absence
#         a_p=35,b_p=300,      # this will produce distribution centered on 0.10 for h2 of menarche
#         rg=0.2,
#         num_calcs=10000)
#