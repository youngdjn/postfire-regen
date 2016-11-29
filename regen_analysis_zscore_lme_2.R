setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

source("R/regen_utils_2_zscore.R")


library(arm)
library(plyr)
library(MuMIn)
library(snow)
library(lme4)
library(ggplot2)
library(gridExtra)

d <- read.csv("Data/regen_clim_full.csv",header=TRUE)
d$FORB <- d$FORBE


d$firesev <- d$FIRE_SEV # use field-assessed fire sev
d <- d[d$firesev %in% c(4,5),] # only medium-high and high severity
d$firesev <- as.factor(d$firesev)

d <- d[!(d$Fire %in% c("ELK","CHINABACK","ZACAM","SHOWERS","DEEP","STRAYLOR","HARDING")),] # rich has some uncertainty in fire and sampling year could be added later; just removed because of uncertainty in fire and sampling year
# kevin says deep, straylor, harding, showers were anomalous (showers was already removed, not sure why)
d$Fire <- as.factor(d$Fire)

# ## remove plots with over 2000 mm precip (most BTU plots)
#d <- d[d$ppt.normal < 2000,]



d$seedtr.comb <- ifelse((is.na(d$seedtr.dist) | (d$seedtr.dist >= 150)),d$dist.to.low,d$seedtr.dist)
d$seedtr.any.comb <- ifelse(is.na(d$seed.tree.any) | (d$seed.tree.any >= 150),d$dist.to.low,d$seed.tree.any)

d$regen.presab <- ifelse(d$regen.count > 0,TRUE,FALSE)

## calculate northness and eastness
d$northness <- cos(deg2rad(d$aspect))
d$eastness <- sin(deg2rad(d$aspect))

d$ppt.normal.sq <- d$ppt.normal^2
d$tmean.normal.sq <- d$tmean.normal^2


#(optional) thin to Sierra fires
sierra.fires <- c("STRAYLOR","CUB","RICH","DEEP","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER")
d$in.sierra <- ifelse(d$Fire %in% sierra.fires,TRUE,FALSE)
d <- d[d$in.sierra,]

# d.a <- d[d$species == "PSME",]
# write.csv(d.a,"welch_plots.csv")

#### Summary plots of climate distribution ####

#early vs. late average precip: plots that were dry during the first three years were wet during the last two years
ggplot(d,aes(x=diff.norm.ppt.z,y=diff.norm.ppt.z.late,color=Fire)) +
  geom_point(size=5) +
  theme_bw(20) +
  labs(x="Average precip anomaly first 3 years",y="Average precip anomaly last 2 years")

#early vs. late minimum precip: not very much correlation here
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.min.z.late,color=Fire)) +
  geom_point(size=5) +
  theme_bw(20) +
  labs(x="Minimum precipitation anomaly first 3 years",y="Minimum precip anomaly last 2 years")

#early minimum precip vs. normal climate: #!may want to leave out Rich
ggplot(d,aes(x=diff.norm.ppt.min.z,y=ppt.normal,color=Fire)) +
  geom_point(size=5)  +
  theme_bw(20) +
  labs(x="Minimum precipitation anomaly first 3 years (SD from mean)",y="Normal precipitation (mm)")

summary(d$ppt.normal)
summary(d$diff.norm.ppt.min.z)

#early minimum precip vs. early average precip
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.z,color=Fire)) +
  geom_point(size=5)  +
  theme_bw(20) +
  labs(x="Minimum precip anomaly first 3 years",y="Average precip anomaly first 3 years")

#early minimum precip vs. late average precip: no correlatio here but #! may want to leave out Rich
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.z.late,color=Fire)) +
  geom_point(size=5) +
  theme_bw(20) +
  labs(x="Minimum precipitation anomaly first 3 years (SD from mean)",y="Average precipitation anomaly in the last 2 years (SD from mean)")




### Calc percentage of plots with hardwoods by fire
d.hw <- d[d$species=="HARDWOOD",]
prop.hw <- aggregate(d.hw$regen.presab,by=list(d.hw$Fire),FUN=mean)
names(prop.hw) <- c("Fire","Prop.hw")

d.conif <- d[d$species=="CONIFER",]
prop.conif <- aggregate(d.conif$regen.presab,by=list(d.conif$Fire),FUN=mean)
names(prop.conif) <- c("Fire","Prop.conif")

prop <- cbind(prop.hw,prop.conif$Prop.conif)
names(prop) <- c("Fire","prop.hw","prop.conif")




#### fire-level summaries ####
library(reshape2)

focal.vars <- c("GRASS","SHRUB","diff.norm.ppt.min.z","diff.norm.ppt.z","northness","seedtr.any.comb","ppt.normal","regen.count","regen.presab")

fire.mean <- aggregate(d[,focal.vars],by=list(d$Fire,d$species),FUN=mean,na.rm=TRUE)
fire.high <- aggregate(d[,focal.vars],by=list(d$Fire,d$species),quantile,probs=.75,na.rm=TRUE)
fire.low <- aggregate(d[,focal.vars],by=list(d$Fire,d$species),quantile,probs=.25,na.rm=TRUE)

names(fire.mean)[1:2] <- c("Fire","species")
names(fire.high)[1:2] <- c("Fire","species")
names(fire.low)[1:2] <- c("Fire","species")

fire.mean.conif <- fire.mean[fire.mean$species=="CONIFER",]
fire.high.conif <- fire.high[fire.high$species=="CONIFER",]
fire.low.conif <- fire.low[fire.low$species=="CONIFER",]

fire.mean.hw <- fire.mean[fire.mean$species=="HARDWOOD",c("regen.count","regen.presab")]
fire.high.hw <- fire.high[fire.high$species=="HARDWOOD",c("regen.count","regen.presab")]
fire.low.hw <- fire.low[fire.low$species=="HARDWOOD",c("regen.count","regen.presab")]

names(fire.mean.hw) <- paste(names(fire.mean.hw),"hw",sep=".")
names(fire.high.hw) <- paste(names(fire.high.hw),"hw",sep=".")
names(fire.low.hw) <- paste(names(fire.low.hw),"hw",sep=".")


fire.nplots <- aggregate(d[,focal.vars[1]],by=list(d$Fire,d$species),FUN=length)
fire.nplots <- fire.nplots[fire.nplots$Group.2=="ABCO","x"]


fire.mean.full <- cbind(fire.mean.conif,fire.mean.hw)
fire.low.full <- cbind(fire.low.conif,fire.low.hw)
fire.high.full <- cbind(fire.high.conif,fire.high.hw)

firenames <- fire.mean.full$Fire

mean.q <- fire.mean.full[,3:ncol(fire.mean.full)]
low.q <- fire.low.full[,3:ncol(fire.low.full)]
high.q <- fire.high.full[,3:ncol(fire.high.full)]

# mean.q <- signif(mean.q,4)
# low.q <- signif(low.q,4)
# high.q <- signif(high.q,4)

mean.q <- as.matrix(round(mean.q,2))
low.q <- as.matrix(round(low.q,2))
high.q <- as.matrix(round(high.q,2))

mean.q <- (signif(mean.q,4))
low.q <- (signif(low.q,4))
high.q <- (signif(high.q,4))


bounds <- paste("(",low.q,",",high.q,")",sep="")
full <- paste(mean.q,bounds)
full <- matrix(full,nrow=nrow(mean.q))
rownames(full) <- firenames
colnames(full) <- colnames(mean.q)
full <- as.data.frame(full)
full$nplots <- fire.nplots


write.csv(full,"fire_summary.csv")


# ## optional: remove BTU, AmRiv, Cub, Rich as potential confounders
# remove.fires <- c("BTU LIGHTENING","AMERICAN RIVER","CUB","RICH")
# d$keep.fires <- ifelse(!(d$Fire %in% remove.fires),TRUE,FALSE)
# d <- d[d$keep.fires,]







#### Prepare to fit models and make counterfactual plots ####


#rescale parameters
d$regen.presab <- ifelse(d$regen.count > 0,TRUE,FALSE)
d$ppt.post.sq <- d$ppt.post^2
d$ppt.post.min.sq <- d$ppt.post.min^2

## center and standardize
leave.cols <- c("regen_presab","regen_count","Regen_Plot","Fire","Year","survey.years.post","SHRUB","FORB","GRASS","diff.norm.ppt.min.z","diff.norm.ppt.z","diff.norm.ppt.z.late")
d.c <- center.df(d,leave.cols=leave.cols)
stan_d1_mean <- mean.df(d,leave.cols=leave.cols)
stan_d1_sd <- sd.df(d,leave.cols=leave.cols)

d.c.focal <- d.c #alternate

# # aggregate by fire and look for correlation among predictors (potential for confounding)
# d.agg <- aggregate(d.c.focal,by=list(d.c.focal$Fire),FUN=mean)
# d.year <- data.frame(Fire=d$Fire,Year=d$Year)
# d.year <- unique(d.year)
# d.agg <- merge(d.agg,d.year,by.x="Group.1",by.y="Fire")
# d.agg$in.sierra <- ifelse(d.agg$Group.1 %in% sierra.fires,1,2)
# plot(diff.norm.ppt.z_c~diff.norm.tmean.z_c,col=in.sierra,d=d.agg)
# plot(diff.norm.ppt.min.z_c~Year,col=in.sierra,d=d.agg)
# d.agg$fire.year <- paste(d.agg$Group.1,d.agg$Year.x-5)
# library(ggplot2)
# library(directlabels)
# p <- ggplot(d.agg, aes(x=diff.norm.ppt.z_c,y=diff.norm.ppt.min.z_c,label=fire.year))+
#   geom_point(size=5) +
#   geom_text(check_overlap=TRUE) +
#   scale_x_continuous(limits=c(-2.5,1.5))
# p + direct.label(p, 'last.points')
# m <- lm(diff.norm.ppt.z_c~diff.norm.ppt.min.z_c,d=d.agg)






### figure out precip counterfactual bins
d.c.btu <- d.c[d.c$Fire == "BTU LIGHTENING",]
b1.upper <- max(d.c.btu$ppt.normal_c)

b1.b2.boundary <- max(d.c.btu$ppt.normal_c)

d.c.rich <- d.c[d.c$Fire == "RICH",]
gap.upper <- min(d.c.rich$ppt.normal_c)

d.c.moon <- d.c[d.c$Fire == "MOONLIGHT",]
gap.lower <- max(d.c.moon$ppt.normal_c)

b2.b3.boundary <- mean(c(gap.upper,gap.lower))

d.c.str <- d.c[d.c$Fire == "STRAYLOR",]
b3.bottom <- min(d.c.str$ppt.normal_c)

## bounds are -1.59, -0.5, 1.252, 3.06

ppt.level1 <- mean(c(b1.upper,b1.b2.boundary))
ppt.level2 <- mean(c(b2.b3.boundary,b1.b2.boundary))
ppt.level3 <- mean(c(b2.b3.boundary,b3.bottom))




#Plot the two climatic (weather) predictors on the transformed scale
ggplot(d.c,aes(x=diff.norm.ppt.min.z,y=ppt.normal_c,color=Fire)) +
  geom_point(size=5)  +
  theme_bw(20) +
  labs(x="Minimum precip anomaly first 3 years",y="Average precip anomaly first 3 years")











#### Fit model and make counterfactual plots ####


p <- list() #for holding counterfactual plots for each vegetation type
coef.df <- data.frame() # for holding coefficient summaries

types <- c("CONIFER","HARDWOOD","SHRUB","GRASS")
for(type in types) {

  if(type=="HARDWOOD" | type=="CONIFER") {
    pres.ab <- TRUE
    d.sp <- d[d$species==type,]
  } else {
    pres.ab <- FALSE
    d.sp <- d[d$species=="ALL",] #species is irrelevant here; just need to reduce to one row per plot
  }
  
  #center and standardize
  leave.cols <- c("regen_presab","regen_count","Regen_Plot","Fire","Year","survey.years.post","SHRUB","FORB","GRASS","diff.norm.ppt.z.late","diff.norm.ppt.z","diff.norm.ppt.min.z")
  d.c.sp <- center.df(d.sp,leave.cols=leave.cols)
  
  
  d.c.sp.focal <- with(d.c.sp,data.frame(regen.presab,northness_c,ppt.normal_c,ppt.normal.sq_c,ppt.post_c,diff.norm.ppt_c,ppt.post.min_c,ppt.post.sq_c,ppt.post.min.sq_c,tmean.normal_c,tmean.normal.sq_c,diff.norm.ppt.z,diff.norm.ppt.min.z,Fire,seedtr.any.comb_c,seedtr.comb_c,eastness_c,diff.norm.tmean.JJA.mean.z_c,diff.norm.tmean.DJF.mean.z_c,diff.norm.tmean.z_c,SHRUB,FORB,GRASS,Year,diff.norm.ppt.z.late))
  d.c.sp.incomplete <- d.c.sp.focal[!complete.cases(d.c.sp.focal),]
  d.c.sp.focal <- d.c.sp.focal[complete.cases(d.c.sp.focal),]
  
  
  
  ### Prepare percent cover data ###
  
  hist(d.c.sp.focal$SHRUB,nclass=30)
  hist(d.c.sp.focal$FORB)
  hist(d.c.sp.focal$GRASS)
  
  d.c.sp.focal$shrub.prop <- d.c.sp.focal$SHRUB/100
  d.c.sp.focal$grass.prop <- d.c.sp.focal$GRASS/100 # could do same for FORB
  
  d.c.sp.focal$shrub.prop <- ifelse(d.c.sp.focal$shrub.prop==0,0.005,d.c.sp.focal$shrub.prop)
  d.c.sp.focal$shrub.prop <- ifelse(d.c.sp.focal$shrub.prop>=1,0.995,d.c.sp.focal$shrub.prop)
  
  d.c.sp.focal$grass.prop <- ifelse(d.c.sp.focal$grass.prop==0,0.005,d.c.sp.focal$grass.prop)
  d.c.sp.focal$grass.prop <- ifelse(d.c.sp.focal$grass.prop==1,0.995,d.c.sp.focal$grass.prop)
  
  d.c.sp.focal$SHRUB <- ifelse(d.c.sp.focal$SHRUB == 0,0.5,d.c.sp.focal$SHRUB)
  d.c.sp.focal$GRASS <- ifelse(d.c.sp.focal$GRASS == 0,0.5,d.c.sp.focal$GRASS)
  
  d.c.sp.focal$shrub.l <- logit(d.c.sp.focal$shrub.prop)
  d.c.sp.focal$grass.l <- logit(d.c.sp.focal$grass.prop)
  
  
  
  # calculate what proportion of plots have the focal species
  #sum(d.c.sp.focal$regen.presab)/nrow(d.c.sp.focal)
  
  m.full <- NULL
  
  if(pres.ab) {
    
    if(type=="CONIFER") {
  
      m.full <- glmer(regen.presab ~ 
                    seedtr.any.comb_c +
                    seedtr.comb_c +
                    northness_c +
                    ppt.normal_c * diff.norm.ppt.z +
                    ppt.normal_c * diff.norm.ppt.min.z +
                    ppt.normal.sq_c +
                    (1|Fire),
                  data=d.c.sp.focal, family=binomial, na.action="na.fail")
      

    } else if(type=="HARDWOOD") {
      
      m.full <- glmer(regen.presab ~ 
                        #seedtr.any.comb_c +
                        #seedtr.comb_c +
                        northness_c +
                        ppt.normal_c * diff.norm.ppt.z +
                        ppt.normal_c * diff.norm.ppt.min.z +
                        ppt.normal.sq_c +
                        (1|Fire),
                      data=d.c.sp.focal, family=binomial, na.action="na.fail")
      
    }

  } else {  #this is a % cover model
    
    if(type=="SHRUB") {
  
      m.full <-lmer(shrub.l ~ 
                      northness_c +
                      ppt.normal_c * diff.norm.ppt.z +
                      ppt.normal_c * diff.norm.ppt.min.z +
                      ppt.normal.sq_c +
                      #tmean.normal_c +
                      (1|Fire),
                  data=d.c.sp.focal, na.action="na.fail", REML=FALSE)
      
    } else if(type=="GRASS") {
      m.full <-lmer(grass.l ~ 
                      northness_c +
                      ppt.normal_c * diff.norm.ppt.z +
                      ppt.normal_c * diff.norm.ppt.min.z +
                      ppt.normal.sq_c +
                      #tmean.normal_c +
                      (1|Fire),
                    data=d.c.sp.focal, na.action="na.fail", REML=FALSE)
    }
   
  
     
  }
  
  
  cl <- makeCluster(4, type = "SOCK")
  clusterEvalQ(cl, library(lme4))
  clusterExport(cl,"d.c.sp.focal")
  a <- pdredge(m.full,cluster=cl, fixed=c("ppt.normal_c"),
               subset=
                 #(!(`seedtr.any.comb_c`&`seedtr.comb_c`)) &
                 #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.JJA.mean.z_c`)) &
                 #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.DJF.mean.z_c`)) &
                 (!(`ppt.normal.sq_c` & !`ppt.normal_c`))
               #!(`tmean.normal.sq_c` & !`tmean.normal_c`))
  )
  stopCluster(cl)
  
  imp <- data.frame(var=names(importance(a)),imp=importance(a))
  best.model.num <- as.character(row.names(a[1,]))
  best.model <- get.models(a,best.model.num)[[1]]
  #best.model <- m.full
  
  print(summary(best.model))
  #best.model.tree <- best.model
  best.coef <- fixef(best.model)
  best.coef <- data.frame(var=names(best.coef),coef=best.coef)
  #imp.coef <- merge(imp,best.coef,by="var",all.x=TRUE)
  #write.csv(imp.coef,"mumin_out.csv",row.names=FALSE)

  
#   mod.names <-row.names(a)
#   mod.weights <- a[,"weight"]
  
  N <- 1000 #number of samples to take
  nsims <- N
  coef.sim <- data.frame()
  
  # # get all fit models  (this section for multi-model inference only)
  # models <- list()
  # 
  # for(i in 1:length(mod.names)) {
  #   
  #   mod.name <- mod.names[i]
  #   models[[mod.name]] <- get.models(a,mod.name)[[1]]
  #   
  # }
  #  this section for multi-model averaging
  # for(i in 1:N) {
  #   
  #   mod.name <- sample(mod.names,prob=mod.weights,replace=TRUE,size=1)
  #   mod.name <- mod.names[1] # if using only the best model (otherwise comment out for multi-model averaging)
  #   model <- models[[mod.name]]
  #   coef.sim.row <- fixef(sim(model,n.sims=1)) # change to get just one simulation, store in matrix
  #   coef.sim.row <- as.data.frame(coef.sim.row)
  #   coef.sim <- rbind.fill(coef.sim,coef.sim.row)
  #   
  # }
  
  # simulate coef for best model. this section for using best-model only
  coef.sim <- fixef(sim(best.model,n.sims=nsims))
  coef.sim[is.na(coef.sim)] <- 0
  
  
  coef.mean <- apply(coef.sim,2,mean)
  coef.low <- apply(coef.sim,2,quantile,prob=0.025)
  coef.high <- apply(coef.sim,2,quantile,prob=0.975)
  
  coef.mean <- signif(coef.mean,3)
  coef.low <- signif(coef.low,3)
  coef.high <- signif(coef.high,3)
  
  coef.summary <- rbind(coef.mean,coef.low,coef.high)
  
  coef.bounds <- paste("(",coef.low,",",coef.high,")",sep="")
  coef.full <- paste(coef.mean,coef.bounds)
  names(coef.full) <- names(coef.mean)
  coef.full <- as.data.frame(t(coef.full))
  coef.full$type <- type
  coef.df <- rbind.fill(coef.df,coef.full)
  
  # could adapt for plotting coef plots (would require standardized coefs to interpret)
#   bin.full.num <- as.data.frame(rbind(bin.median,bin.lwr,bin.upr))
#   bin.full.num$val <- c("median","lower","upper")
#   bin.full.num$year <- year
#   bin.coef.num.df <- rbind(bin.coef.num.df,bin.full.num)
  
  
  
  
  ################################
  #### Make predictions ##########
  ################################
  
  #link ppt.normal to ppt.sq
  
  ppt.link <- data.frame(ppt.normal_c=d.c.sp.focal$ppt.normal_c,ppt.normal.sq_c=d.c.sp.focal$ppt.normal.sq_c)
  m <- lm(ppt.normal.sq_c ~ ppt.normal_c + I(ppt.normal_c^2),ppt.link)
  
  ppt.link.a <- coef(m)[1]
  ppt.link.b <- coef(m)[2]
  ppt.link.b2 <- coef(m)[3]
  
  
  # -0.1721 + 0.8811x + 0.1722x^2
  
  
  # prediction matrix (newdata)
  
  #determine appropriate "counterfactual" values for predictors
  summary(d.c.sp.focal$diff.norm.ppt.min.z) # -1.475 to 2.015
  summary(d.c.sp.focal$ppt.normal_c) # -1.5 to 3
  summary(d.c.sp.focal$seedtr.any.comb_c) #-1 to 4.9
  hist(d.c.sp.focal$seedtr.any.comb_c) #-1 to 4.9
  hist(d.c.sp.focal$northness_c) #-1 to 4.9
  
  # set "counterfactual" values for predictors #-1.045, 0.376, 2.156
  ppt.normal.c.low <- ppt.level3
  ppt.normal.c.mid <- ppt.level2
  ppt.normal.c.high <- ppt.level1
  
  plot(ppt.normal_c~diff.norm.ppt.min.z,col=Fire,data=d.c.sp.focal)
  
  #get range of min precip anomaly for each ppt value
  d.c.sp.focal.lowprecip <- d.c.sp.focal[(d.c.sp.focal$ppt.normal_c > b3.bottom) & (d.c.sp.focal$ppt.normal_c < b2.b3.boundary),]
  low.min <- min(d.c.sp.focal.lowprecip$diff.norm.ppt.min.z)
  low.max <- max(d.c.sp.focal.lowprecip$diff.norm.ppt.min.z)
  
  d.c.sp.focal.midprecip <- d.c.sp.focal[(d.c.sp.focal$ppt.normal_c < b1.b2.boundary) & (d.c.sp.focal$ppt.normal_c > b2.b3.boundary),]
  mid.min <- min(d.c.sp.focal.midprecip$diff.norm.ppt.min.z)
  mid.max <- max(d.c.sp.focal.midprecip$diff.norm.ppt.min.z)
  
#   d.c.sp.focal.highprecip <- d.c.sp.focal[(d.c.sp.focal$ppt.normal_c > b1.b2.boundary) & (d.c.sp.focal$ppt.normal_c < b1.upper),]
#   high.min <- min(d.c.sp.focal.highprecip$diff.norm.ppt.min.z)
#   high.max <- max(d.c.sp.focal.highprecip$diff.norm.ppt.min.z)
#   
  
  
  diff.norm.min.seq <- seq(from=low.min,to=low.max,length.out=100)
  
  pred.df.lowprecip <- data.frame(
    
    diff.norm.ppt.min.z = diff.norm.min.seq,
    ppt.normal_c =ppt.normal.c.low,
    ppt.normal.sq_c = ppt.link.a + ppt.link.b*ppt.normal.c.low + ppt.link.b2*(ppt.normal.c.low^2),
    northness_c= 0,
    diff.norm.ppt.z = 0,
    seedtr.comb_c = 0,
    seedtr.any.comb_c = 0,
    tmean.normal_c = 0
  
  )
  
  diff.norm.min.seq <- seq(from=mid.min,to=mid.max,length.out=100)
  
  pred.df.midprecip <- data.frame(
    
    diff.norm.ppt.min.z = diff.norm.min.seq,
    ppt.normal_c =ppt.normal.c.mid,
    ppt.normal.sq_c = ppt.link.a + ppt.link.b*ppt.normal.c.mid + ppt.link.b2*(ppt.normal.c.mid^2),
    northness_c= 0,
    diff.norm.ppt.z = 0,
    seedtr.comb_c = 0,
    seedtr.any.comb_c = 0,
    tmean.normal_c = 0
    
  )
  
#   diff.norm.min.seq <- seq(from=high.min,to=high.max,length.out=100)
#   
#   pred.df.highprecip <- data.frame(
#     
#     diff.norm.ppt.min.z = diff.norm.min.seq,
#     ppt.normal_c =ppt.normal.c.high,
#     ppt.normal.sq_c = -0.1721 + 0.8811*ppt.normal.c.high + 0.1722*ppt.normal.c.high^2,
#     northness_c= 0,
#     diff.norm.ppt.z = 0,
#     seedtr.comb_c = 0,
#     seedtr.any.comb_c = 0,
#     tmean.normal_c = 0
#     
#   )
  
  pred.df <- rbind(pred.df.lowprecip,pred.df.midprecip)#,pred.df.highprecip)
  
  # 
  # #### alternative pred df for 3-d plotting: all potential combinations of min and average (need to restrict z.seq ranges)
  # pairwise <- expand.grid(z.seq,z.seq)
  # pred.df <- data.frame(
  #   diff.norm.ppt.min.z = pairwise[,1],
  #   ppt.normal_c = 0,
  #   ppt.normal.sq_c = -0.186 + 0.846*0 + 0.186*0^2,
  #   northness_c= 0,
  #   diff.norm.ppt.z = pairwise[,2]
  # )
  # ### end alternative pred df ####
  
  
  # interact.df <- matrix(data=NA,nrow=nsims,ncol=ncol(pred.df)^2)
  # names.df <- vector(length=ncol(pred.df)^2)
  
  #### add interaction prediction terms to the df (all pairsiwe interactions)
  
  interact.df <- data.frame("Fake"=rep(NA,nrow(pred.df)))
  
  for(i in 1:ncol(pred.df)) { # for each col
    
    for(j in 1:ncol(pred.df)) {
      
      name.i <- names(pred.df)[i]
      name.j <- names(pred.df)[j]
      name.inter <- paste(name.i,name.j,sep=":")
      val.inter <- pred.df[,i] * pred.df[,j]
      
      interact.df <- cbind(interact.df,val.inter)
      names(interact.df)[ncol(interact.df)] <- name.inter
      
      
    }
  
  }
  
  pred.df <- cbind(pred.df,interact.df)
  pred.df$"(Intercept)" <- 1 # this is the "predictor" value to multiple the intercept by
  
  
  # ### predict for all tree seedlings
  # library(merTools)
  # 
  # predFun <- function(x)predict(x, pred.df, re.form=NA)
  # 
  # p <- predict(best.model.tree,newdata=pred.df,type="response")
  # p2 <- bootMer(best.model.tree,nsim=100,FUN=predFun)
  # #predMat <- p2$t
  # 
  
  coef.names <- c("(Intercept)",colnames(coef.sim)[-1])
  
  pred.order <- pred.df[,coef.names]
  
  #pred.terms <- coef.sim * pred.order
  
  pred.fit <- data.frame(mean=rep(NA,nrow(pred.df)),lwr=NA,upr=NA)
  
  for (i in 1:nrow(pred.df)) {
    
    # for each row in the prediction DF, calculate the distribution of fitted response
    pred.mat <- matrix(rep(as.matrix(pred.order)[i,],each=nrow(coef.sim)),nrow=nrow(coef.sim))
    a <- as.matrix(coef.sim)
    b <- as.matrix(pred.mat)
    
    pred.terms <- a * b
    pred.y <- rowSums(pred.terms)
    pred.y.resp <- invlogit(pred.y) #both cover and probability were logit-transformed, so transform back to response scale
    pred.y.resp.mean <- mean(pred.y.resp)
    pred.y.resp.lwr <- quantile(pred.y.resp,probs=0.025)
    pred.y.resp.upr <- quantile(pred.y.resp,probs=0.975)
    
    ret <- data.frame(mean=pred.y.resp.mean,lwr=pred.y.resp.lwr,upr=pred.y.resp.upr)
    
    pred.fit[i,] <- ret
  
  }
  
  pred.fit <- pred.fit*100 # transform to percent
  pred.fit <- cbind(pred.df,pred.fit) # merge predictors and simulated responses
  
  ylab <- NULL
  if(pres.ab) {
    ylab <- "Probability of regeneration"
    ylim <- 60
  } else {
    ylab <- "Percent cover"
    ylim <- 60
  }
  
  
  #### plot counterfactual predictions ####
  
  title <- paste(type,"- Low mean precip")
  
  
  # low normal precip
  d.plot <- pred.fit[pred.fit$ppt.normal_c==ppt.normal.c.low,]
  p.low <- ggplot(d.plot,aes(x=diff.norm.ppt.min.z,y=mean)) +
           geom_line(color="darkgreen") +
           geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
            coord_cartesian(ylim=c(0, ylim)) +
           scale_x_continuous(limits = c(-1.25,0)) +
           xlab("Postfire precipitation (3-yr minimum)") +
           ylab(ylab) +
           theme_bw(13) +
           ggtitle(title)
  
  title <- paste(type,"- High mean precip")
  
  
  # mid normal precip
  d.plot <- pred.fit[pred.fit$ppt.normal_c==ppt.normal.c.mid,]
  p.mid <- ggplot(d.plot,aes(x=diff.norm.ppt.min.z,y=mean)) +
    geom_line(color="darkgreen") +
    geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
    coord_cartesian(ylim=c(0, ylim)) +
    scale_x_continuous(limits = c(-1.25,0)) +
    xlab("Postfire precipitation (3-yr minimum)") +
    ylab(ylab) +
    theme_bw(13) +
    ggtitle(title)
  
#   # high normal precip
#   d.plot <- pred.fit[pred.fit$ppt.normal_c==ppt.normal.c.high,]
#   p.high <- ggplot(d.plot,aes(x=diff.norm.ppt.min.z,y=mean)) +
#     geom_line(color="darkgreen") +
#     geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
#     scale_x_continuous(limits = c(-1.25,0)) +
#     coord_cartesian(ylim=c(0, 50)) +
#     xlab("Postfire precipitation (3-yr minimum)") +
#     ylab(ylab) +
#     theme_bw(13)
  
  
  p[[type]] <- grid.arrange(p.low,p.mid,nrow=1)

}


p.full <- grid.arrange(p[["CONIFER"]],p[["HARDWOOD"]],p[["SHRUB"]],p[["GRASS"]],nrow=4)

coef.df
write.csv(coef.df,"model_coefs.csv",row.names=FALSE)


library(Cairo)
Cairo(file="Fig2_moderate_sev.png",typ="png",width=600,height=1000)
p.full <- grid.arrange(p[["CONIFER"]],p[["HARDWOOD"]],p[["SHRUB"]],p[["GRASS"]],nrow=4)
dev.off()








#### 3-d plot no longer using ####

# 
# #plot 3-d
# library(ggplot2)
# library(viridis)
# 
# p <- ggplot(pred.fit,aes(x=diff.norm.ppt.z,y=diff.norm.ppt.z_c)) +
#   geom_tile(aes(fill=mean)) +
#   #geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
#   #scale_y_continuous(limits = c(-6,100)) +
#   #xlab("Postfire precipitation (3-yr average)") +
#   #ylab("Percent cover") +
#   theme_bw(13) +
#   scale_fill_viridis(limits=c(0,100))
# 
# 
# library(Cairo)
# Cairo(file="test2.png",typ="png",width=400,height=300)
# p
# dev.off()





