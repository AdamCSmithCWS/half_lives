

### using earlier version of this before ad-hoc fix for zero-value estimates of Q 
orig = read.csv("output/original trajectories extinction risks_pre_adhoc.csv")
orig = orig[which(orig$Q > 0.00001 & orig$R > 0.00001),]
q = orig$Q
r = orig$R/2 ## half of the observation variance see rationale below
u = abs(orig$U)

lq = log(q)
lr = log(r)
plot(lr,lq)
abline(0,1)
# approximately linear and very close to 1:1 
# probably very close to 1, since we know there is error in both estimates
## therefore the mean of the ratio of the logs is a reasonable way to approximate one from the other
fix = mean(lq/lr)
## so on average across species, on the log scale, the MARSS function is partitioning half of the variance to 
## process and half to observation
## this relatively equal split is not surprising given that the model has no data that support
## a clear partition

### Given that this relationship is for trajectories where some of the error has been
### allocated to each process and observation
### to fix the zero-estimated process variance, our best estimate is approximately half of the variance estimate


### function to fix values of q, based on half of r
fix_fx <- function(r){
  q <- exp(fix*log(r/2))
  return(q)
}

# directory to place results
res.dir.gen <- 'results'

#############  Load Data  ################################################################ 
spp.indall <- read.csv('C:/Estimating_Change_in_NorthAmerican_Birds/Rosenberg et al annual indices of abundance.csv',stringsAsFactors = F) 


# lookup file relating AOU# to species name
spdf=read.fwf(file='species 2011.txt',
              width = c(11,30,31),
              col.names = c("AOU", "English_Common_Name", "Species"),
              colClasses = c("integer","character","character"),
              strip.white = T)
load('outnames.R')

library(MARSS)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
##########################################################################################
###################################################
###   LOAD FUNCTIONS
an=function(x) as.numeric(x)  #as a number
ran=function(x,dec=2) round(as.numeric(x),dec)  #round a number

#######       Q.EX Function
#Prob. of decline to given threshold (quasi-extinction)
# qe.level = ln(threshold abundance value)
# qeyrs = #of years (timesteps) to project (quasiextinction with qeyrs years)
# u = trend
# q = process error
# x0 = ln(current abundance index)

QE=function(qe.level,qeyrs,u,q, x0){  #Quasi-Extinction risk estimate based on Dennis et al. '91
  u=as.numeric(u)
  q=as.numeric(q)
  q=ifelse(q==0,0.00001,q)
  if(length(qeyrs)==1){tyrs=1:qeyrs}else tyrs=qeyrs
  xd=x0-qe.level 
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/q)) #q=sigma2
  pis=NULL
  for(i in seq_along(tyrs)){
    if(is.finite(exp(2*xd*abs(u)/q))){
      sec.part=exp(2*xd*abs(u)/q)*
        pnorm((-xd-abs(u)*
                 tyrs[i])/sqrt(q*tyrs[i]))}else sec.part=0
    
    pis[i]=round(min(c(1,p.ever * pnorm((-xd+abs(u)*tyrs[i])/sqrt(q*tyrs[i]))+
                   sec.part)),3)
  }##for i
  return(pis)
}## End QE


#Prob. of decline by PD percent from current abundance after qeyrs years
# PD = percent decline from current abundance

HL=function(qeyrs,u,q,PD){  #Quasi-Extinction risk estimate based on Dennis et al. '91
  u=as.numeric(u)
  q=as.numeric(q)
  q=ifelse(q==0,0.00001,q)
  if(length(qeyrs)==1){tyrs=1:qeyrs}else tyrs=qeyrs
  xd=log(100/(100-PD)) 
  p.ever = ifelse(u<=0,1,exp(-2*u*xd/q)) #q=sigma2
  pis=NULL
  for(i in seq_along(tyrs)){
    if(is.finite(exp(2*xd*abs(u)/q))){
      sec.part=exp(2*xd*abs(u)/q)*
        pnorm((-xd-abs(u)*
                 tyrs[i])/sqrt(q*tyrs[i]))}else sec.part=0
    
    pis[i]= round(min(c(1,p.ever * pnorm((-xd+abs(u)*tyrs[i])/sqrt(q*tyrs[i]))+
                   sec.part)),3)
  }##for i
  return(pis)
}## End HL
##########################################################################################
out.names[3] <- 'loca'
#Set up Datafile to export results at each loop
write.table(t(out.names),file=paste(res.dir.gen, '/sp.state.csv', sep=''), quote=FALSE, 
            sep=',', row.names=FALSE, col.names=FALSE)

spp.indall$English_Common_Name <- as.character(spp.indall$species)
spp.indall$yr <- spp.indall$year-1969
sp1 = unique(spp.indall$English_Common_Name)
#sp2 = unique(spdf$English_Common_Name)




newpal <- viridis::viridis_pal(option = "D")

cl.orig = newpal(6)[6]
cl.pred = newpal(6)[5]
cl.qe50 = newpal(6)[4]
cl.hl30 = newpal(6)[3]
cl.hl50 = newpal(6)[2]
cl.hl70 = newpal(6)[1]

spp.ind12 <- spp.indall[which(spp.indall$year > 2017-13),] #includes years 2005 - 2017, so minimum trajectory includes 11-years and 10-year trend

sp.states.i <- unique(spp.indall[,c("species","English_Common_Name","firstyear","lastyear")])
# 


 
spp.ind <- spp.indall

plotsout <- list()
length(plotsout) <- length(sp1)
names(plotsout) <- sp1

# Long-term loop ----------------------------------------------------------

projj_out <- NULL

for (i in 1:length(sp1)){ #loop through each species
  
  projj_sp <- NULL
  
  
  print("NEXT SPECIES")  
  print(Sys.time())
  ####  i<-round(runif(1,1,nrow(aou)),0) ### For Testing
  aou.i <- sp1[i]

  worig <- which(spp.ind$English_Common_Name == aou.i & !is.na(spp.ind$index.raw))
  d <- spp.ind[worig,]
  
  d$RCI <- (d$uci.raw-d$lci.raw)/d$index.raw #Relative CV-ish ~ SE*4/mean 
  maxyear <- max(d$yr)
  maxyearn <- unique(d$lastyear)
  yrs <- d$yr
  sp<- as.character(aou.i)
  sp<- str_trim(sp,side='right')
  reg = "Continental"
  df <- d$index.raw
  ci <- d$RCI

  ############## defined as mean relative confidence index greater than 2
  w.out = which(sp.states.i$sp == aou.i)
  sp.states.i[w.out,"ave.ind"] <- mean(df,na.rm = T)
  sp.states.i[w.out,"max.ind"] <- max(df,na.rm = T)
  sp.states.i[w.out,"RCI"] <- mean(ci,na.rm = T)
  sp.states.i[w.out,"Uncertain"] <- sp.states.i[w.out,"RCI"]>2|!is.finite(sp.states.i[w.out,"RCI"])
  

  df<-log(df)
  ###### MARSS #########
  mod.struc <- list(R='diagonal and equal', Q='diagonal and unequal')
  print(sp)
  mod.out <- MARSS(df, model=mod.struc,
                   control=list(conv.test.slope.tol=0.1,maxit=5000, trace=1, allow.degen=TRUE))
  

# ad hoc fix for process variance estimates = 0 ---------------------------

#### solution here is to use the geometric mean ratio of process variance to trend = 0.06151012
  ### this was calculated before this loop was added, using only the 395 species for which
  ### non-zero values of Q were estimated in the long-term process
  
  if(mod.out$par$Q < 0.0001){
  newq = fix_fx(mod.out$par$R)
  
  mq <- matrix(newq,nrow = 1,ncol = 1)
  
  mod.struc <- list(R='diagonal and equal', Q=mq)
  
  mod.out <- MARSS(df, model=mod.struc,
                   control=list(conv.test.slope.tol=0.1,maxit=5000, trace=1, allow.degen=TRUE))
  
  mod.out.CI=MARSSparamCIs(mod.out)
  
  mod.out$par$Q <- newq
  mod.out.CI$par.lowCI$Q <- NA
  mod.out.CI$par.upCI$Q<- NA
  
  }else{
  
  mod.out.CI=MARSSparamCIs(mod.out)
  }

  ### Add MARSS results to output
  sp.states.i[w.out,"conv"] <- mod.out$convergence
  sp.states.i[w.out,"u"] <- mod.out$par$U
  sp.states.i[w.out,"u.lowCI"] <- mod.out.CI$par.lowCI$U
  sp.states.i[w.out,"u.upCI"] <- mod.out.CI$par.upCI$U
  sp.states.i[w.out,"R"] <- mod.out$par$R
  sp.states.i[w.out,"R.lowCI"] <-mod.out.CI$par.lowCI$R
  sp.states.i[w.out,"R.upCI"] <-mod.out.CI$par.upCI$R
  sp.states.i[w.out,"Q"] <- mod.out$par$Q
  sp.states.i[w.out,"Q.lowCI"] <-mod.out.CI$par.lowCI$Q
  sp.states.i[w.out,"Q.upCI"] <-mod.out.CI$par.upCI$Q
  st.indx <- exp(mod.out$states)
  st.indx.uci <- exp(mod.out$states + 1.96*mod.out$states.se)
  st.indx.lci <- exp(mod.out$states - 1.96*mod.out$states.se)
  spp.ind[worig,"st.indx"] <- as.numeric(st.indx)
  spp.ind[worig,"st.indx.uci"] <- as.numeric(st.indx.uci)
  spp.ind[worig,"st.indx.lci"] <- as.numeric(st.indx.lci)
  

  #####  Calculate QE and HL projections  ###########
  qex.thresh <- 0.01 # Quasi-extinction threshold
    qe.pred.j <- QE(qe.level=log(qex.thresh),
                    qeyrs=50, 
                    u=mod.out$par$U, 
                    q=mod.out$par$Q, 
                    x0=log(st.indx[length(st.indx)]))
    spp.ind[worig,"qe.pred"] <- qe.pred.j[rev(51-(yrs))]   
    
    sp.states.i[w.out,"predicted_year"] <- maxyearn+50
    sp.states.i[w.out,paste0("qe.pred_50")] <-qe.pred.j[length(qe.pred.j)]
    
   
    for(hlp in c(30,50,70)){
    hl.pred.j <- HL(qeyrs=50, 
                    u=mod.out$par$U, 
                    q=mod.out$par$Q,
                    PD = hlp)
    projj <- tibble(species = aou.i,year = c((maxyearn+1):(maxyearn+50)),
                    p.dec = as.numeric(hl.pred.j),
                    Decline = paste0(hlp,"_%"))
    projj_sp <- bind_rows(projj_sp,projj)
    
    spp.ind[worig,paste0("p.hl.",hlp)] <- hl.pred.j[rev(51-(yrs))]
  sp.states.i[w.out,paste0("p.",hlp,"_percent_decline_50_years")] <-hl.pred.j[length(hl.pred.j)]
  if(any(hl.pred.j > 0.5)){
  sp.states.i[w.out,paste0("y.hl.",hlp,"_percent_decline")] <- min(which(hl.pred.j > 0.5))
  }else{
    sp.states.i[w.out,paste0("y.hl.",hlp,"_percent_decline")] <- NA 
  }
    }
  ####### Make and Export PLots
  # n.states <- 17+ncol(df)
  # res.fit <- sp.states.i[,c(3, 18:n.states)]
  # res.fit <- melt(res.fit, id.vars='loca',variable.name='year',value.name='fit.indx' )
  # res.fit$year <- an(gsub("[^0-9]", "",res.fit$year)) + 1969
  # res.indx <- data.frame(loca = sp.states.i[,3], df)
  # res.indx <- melt(res.indx)
  # names(res.indx)<-c('loca','year', 'bbs.indx')
  # res.indx$bbs.indx <- exp(res.indx$bbs.indx)
  # res.indx$year <-an(res.indx$year)+1969
  # res <- merge(res.indx, res.fit)
 
  res = spp.ind[worig,]
  upy <- 0.95*(max(max(res$uci.raw),max(res$st.indx.uci)))
  
  res$QE <- res$qe.pred*(upy)
  
  for(hlp in c(30,50,70)){
  res[,paste0("HL_",hlp)] <- res[,paste0("p.hl.",hlp)]*(upy)
  }
  
  lbl = res[nrow(res),]
  lbl$year = 2020
  
  plot.name <- Sp.output <-paste(res.dir.gen,'/',sp,'.png', sep="")
  ppi <- 400
  png(file=plot.name, width=8*ppi, height=8*ppi, res=ppi)
  print(
    p <- ggplot(res, aes(x=year, y=index.raw))+ #geom_point(shape=21, size=2.5) +
      #facet_wrap(~ loca, nrow=4, scales='free_y') + 
      theme_bw() + 
      geom_ribbon(aes(x = year,ymax = uci.raw,ymin = lci.raw),alpha = 0.2,fill = cl.orig)+
      geom_ribbon(aes(x = year,ymax = st.indx.uci,ymin = st.indx.lci),alpha = 0.2,fill = cl.pred)+
      geom_line(aes(y=index.raw),colour = cl.orig) +
      geom_line(aes(y=st.indx),colour = cl.pred) +
      # geom_line(aes(y=QE),colour = cl.qe50) +
      # geom_line(aes(y=HL_30),colour = cl.hl30) +
      # geom_line(aes(y=HL_50),colour = cl.hl50) +
      # geom_line(aes(y=HL_70),colour = cl.hl70) +
      annotate(geom = "text",x = lbl$year,y = lbl$QE,label = "QE",colour = cl.qe50)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_30,label = "HL_30",colour = cl.hl30)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_50,label = "HL_50",colour = cl.hl50)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_70,label = "HL_70",colour = cl.hl70)+
      annotate(geom = "text",x = 2010,y = upy*1.02,label = "Probability = 1.0",colour = cl.hl70)+
      annotate(geom = "line",x = c(2000:2020),y = rep(upy,21),colour = cl.hl70)+
      scale_y_continuous(limits = c(0,upy*1.06))+
      
      ylab('Original index') + ggtitle(sp)
  )
  dev.off()
  
  pp <- ggplot(data = projj_sp, aes(x=year))+ 
    theme_bw() + 
     geom_line(aes(y=p.dec,colour = Decline)) +
      scale_y_continuous(limits = c(0,1))+
    scale_color_viridis_d(end = 0.8)+
    
    ylab('probability of decline') + ggtitle(sp)
  
  
  po <- p + pp
  
  plotsout[[i]] <- po
  
  projj_out <- bind_rows(projj_out,projj_sp)
  
  
  write.csv(sp.states.i,"output/original trajectories extinction risks.csv")
  
  
}

# end of species loop -----------------------------------------------------


# sink()
spp.ind$Prediction_year = spp.ind$year + 50
write.csv(spp.ind,"output/original data w annual predictions.csv")
pdf("output/Original trajectories and half-life projections.pdf",
    width = 11,height = 7)
for(i in 1:length(plotsout)){
  print(plotsout[[i]])
}
dev.off()








# short-term projections --------------------------------------------------

spp.ind <- spp.ind12
spp.ind$yr <- spp.ind$year-2004

sp.states.i <- unique(spp.indall[,c("species","English_Common_Name","firstyear","lastyear")])


plotsout <- list()
length(plotsout) <- length(sp1)
names(plotsout) <- sp1


projj_out <- NULL

for (i in 1:length(sp1)){ #loop through each species
  
  projj_sp <- NULL
  
  
  print("NEXT SPECIES")  
  print(Sys.time())
  ####  i<-round(runif(1,1,nrow(aou)),0) ### For Testing
  aou.i <- sp1[i]
  
  worig <- which(spp.ind$English_Common_Name == aou.i & !is.na(spp.ind$index.raw))
  d <- spp.ind[worig,]
  
  d$RCI <- (d$uci.raw-d$lci.raw)/d$index.raw #Relative CV-ish ~ SE*4/mean 
  maxyear <- max(d$yr)
  maxyearn <- unique(d$lastyear)
  yrs <- d$yr
  sp<- as.character(aou.i)
  sp<- str_trim(sp,side='right')
  reg = "Continental"
  df <- d$index.raw
  ci <- d$RCI
  
   w.out = which(sp.states.i$sp == aou.i)
  
  sp.states.i[w.out,"ave.ind"] <- mean(df,na.rm = T)
  sp.states.i[w.out,"max.ind"] <- max(df,na.rm = T)
  sp.states.i[w.out,"RCI"] <- mean(ci,na.rm = T)
  sp.states.i[w.out,"Uncertain"] <- sp.states.i[w.out,"RCI"]>2|!is.finite(sp.states.i[w.out,"RCI"])
  
    df<-log(df)
  ###### MARSS #########
  mod.struc <- list(R='diagonal and equal', Q='diagonal and unequal')
  print(sp)
  mod.out <- MARSS(df, model=mod.struc,
                   control=list(conv.test.slope.tol=0.1,maxit=5000, trace=1, allow.degen=TRUE))
  
  # ad hoc fix for process variance estimates = 0 ---------------------------
  
  #### solution here is to use the geometric mean ratio of process variance to trend = 0.06151012
  ### this was calculated before this loop was added, using only the 395 species for which
  ### non-zero values of Q were estimated in the long-term process above
  
  
  if(mod.out$par$Q < 0.0001){
    newq = fix_fx(mod.out$par$R)
    
    mq <- matrix(newq,nrow = 1,ncol = 1)
    
    mod.struc <- list(R='diagonal and equal', Q=mq)
    
    mod.out <- MARSS(df, model=mod.struc,
                     control=list(conv.test.slope.tol=0.1,maxit=5000, trace=1, allow.degen=TRUE))
    
    mod.out.CI=MARSSparamCIs(mod.out)
    
    mod.out$par$Q <- newq
    mod.out.CI$par.lowCI$Q <- NA
    mod.out.CI$par.upCI$Q<- NA
    
  }else{
    
    mod.out.CI=MARSSparamCIs(mod.out)
  }
  
  ### Add MARSS results to output
  sp.states.i[w.out,"conv"] <- mod.out$convergence
  sp.states.i[w.out,"u"] <- mod.out$par$U
  sp.states.i[w.out,"u.lowCI"] <- mod.out.CI$par.lowCI$U
  sp.states.i[w.out,"u.upCI"] <- mod.out.CI$par.upCI$U
  sp.states.i[w.out,"R"] <- mod.out$par$R
  sp.states.i[w.out,"R.lowCI"] <-mod.out.CI$par.lowCI$R
  sp.states.i[w.out,"R.upCI"] <-mod.out.CI$par.upCI$R
  sp.states.i[w.out,"Q"] <- mod.out$par$Q
  sp.states.i[w.out,"Q.lowCI"] <-mod.out.CI$par.lowCI$Q
  sp.states.i[w.out,"Q.upCI"] <-mod.out.CI$par.upCI$Q
  st.indx <- exp(mod.out$states)
  st.indx.uci <- exp(mod.out$states + 1.96*mod.out$states.se)
  st.indx.lci <- exp(mod.out$states - 1.96*mod.out$states.se)
  spp.ind[worig,"st.indx"] <- as.numeric(st.indx)
  spp.ind[worig,"st.indx.uci"] <- as.numeric(st.indx.uci)
  spp.ind[worig,"st.indx.lci"] <- as.numeric(st.indx.lci)
  
  
  #####  Calculate QE and HL projections  ###########
  qex.thresh <- 0.01 # Quasi-extinction threshold
  qe.pred.j <- QE(qe.level=log(qex.thresh),
                  qeyrs=50, 
                  u=mod.out$par$U, 
                  q=mod.out$par$Q, 
                  x0=log(st.indx[length(st.indx)]))
  spp.ind[worig,"qe.pred"] <- qe.pred.j[rev((50+1)-(yrs))]   
  
  sp.states.i[w.out,"predicted_year"] <- maxyearn+50
  sp.states.i[w.out,paste0("qe.pred_50")] <-qe.pred.j[length(qe.pred.j)]
  
  
  for(hlp in c(30,50,70)){
    hl.pred.j <- HL(qeyrs=50, 
                    u=mod.out$par$U, 
                    q=mod.out$par$Q,
                    PD = hlp)
    
    projj <- tibble(species = aou.i,year = c((maxyearn+1):(maxyearn+50)),
                    p.dec = as.numeric(hl.pred.j),
                    Decline = paste0(hlp,"_%"))
    projj_sp <- bind_rows(projj_sp,projj)
    
    spp.ind[worig,paste0("p.hl.",hlp)] <- hl.pred.j[rev((50+1)-(yrs))]
    sp.states.i[w.out,paste0("p.",hlp,"_percent_decline_50_years")] <-hl.pred.j[length(hl.pred.j)]
    if(any(hl.pred.j > 0.5)){
      sp.states.i[w.out,paste0("y.hl.",hlp,"_percent_decline")] <- min(which(hl.pred.j > 0.5))
    }else{
      sp.states.i[w.out,paste0("y.hl.",hlp,"_percent_decline")] <- NA 
    }
  }
  ####### Make and Export PLots
  # n.states <- 17+ncol(df)
  # res.fit <- sp.states.i[,c(3, 18:n.states)]
  # res.fit <- melt(res.fit, id.vars='loca',variable.name='year',value.name='fit.indx' )
  # res.fit$year <- an(gsub("[^0-9]", "",res.fit$year)) + 1969
  # res.indx <- data.frame(loca = sp.states.i[,3], df)
  # res.indx <- melt(res.indx)
  # names(res.indx)<-c('loca','year', 'bbs.indx')
  # res.indx$bbs.indx <- exp(res.indx$bbs.indx)
  # res.indx$year <-an(res.indx$year)+1969
  # res <- merge(res.indx, res.fit)
  
  res = spp.ind[worig,]
  upy <- 0.95*(max(max(res$uci.raw),max(res$st.indx.uci)))
  
  res$QE <- res$qe.pred*(upy)
  
  for(hlp in c(30,50,70)){
    res[,paste0("HL_",hlp)] <- res[,paste0("p.hl.",hlp)]*(upy)
  }
  
  lbl = res[nrow(res),]
  lbl$year <- 2020
  
  plot.name <- Sp.output <-paste(res.dir.gen,'/',sp,'short.png', sep="")
  ppi <- 400
  png(file=plot.name, width=8*ppi, height=8*ppi, res=ppi)
  print(
    p <- ggplot(res, aes(x=year, y=index.raw))+ #geom_point(shape=21, size=2.5) +
      #facet_wrap(~ loca, nrow=4, scales='free_y') + 
      theme_bw() + 
      geom_ribbon(aes(x = year,ymax = uci.raw,ymin = lci.raw),alpha = 0.2,fill = cl.orig)+
      geom_ribbon(aes(x = year,ymax = st.indx.uci,ymin = st.indx.lci),alpha = 0.2,fill = cl.pred)+
      geom_line(aes(y=index.raw),colour = cl.orig) +
      geom_line(aes(y=st.indx),colour = cl.pred) +
      # geom_line(aes(y=QE),colour = cl.qe50) +
      # geom_line(aes(y=HL_30),colour = cl.hl30) +
      # geom_line(aes(y=HL_50),colour = cl.hl50) +
      # geom_line(aes(y=HL_70),colour = cl.hl70) +
      annotate(geom = "text",x = lbl$year,y = lbl$QE,label = "QE",colour = cl.qe50)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_30,label = "HL_30",colour = cl.hl30)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_50,label = "HL_50",colour = cl.hl50)+
      annotate(geom = "text",x = lbl$year,y = lbl$HL_70,label = "HL_70",colour = cl.hl70)+
      annotate(geom = "text",x = 2010,y = upy*1.02,label = "Probability = 1.0",colour = cl.hl70)+
      annotate(geom = "line",x = c(2005:2020),y = rep(upy,16),colour = cl.hl70)+
      scale_y_continuous(limits = c(0,upy*1.06))+
      
      ylab('Original index') + ggtitle(sp)
  )
  dev.off()
  
  
  pp <- ggplot(data = projj_sp, aes(x=year))+ 
    theme_bw() + 
    geom_line(aes(y=p.dec,colour = Decline)) +
    scale_y_continuous(limits = c(0,1))+
    scale_color_viridis_d(end = 0.8)+
    
    ylab('probability of decline') + ggtitle(sp)
  
  
  po <- p + pp
  
  plotsout[[i]] <- po
  
  projj_out <- bind_rows(projj_out,projj_sp)
  
  write.csv(sp.states.i,"output/original trajectories extinction risks short-term.csv")
}
# sink()
spp.ind$Prediction_year = spp.ind$year + 50
write.csv(spp.ind,"output/original data w annual predictions short-term.csv")
pdf("output/Original trajectories and half-life projections short-term.pdf")
for(i in 1:length(plotsout)){
  print(plotsout[[i]])
}
dev.off()




