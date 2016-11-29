setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

source("R/regen_utils_2_zscore.R")

d <- read.csv("Data/regen_clim_full.csv",header=TRUE)
d$FORB <- d$FORBE


d$firesev <- d$FIRE_SEV # use field-assessed fire sev
d <- d[d$firesev %in% c(4,5),] # only medium-high and high severity
d$firesev <- as.factor(d$firesev)

d <- d[!(d$Fire %in% c("ELK","CHINABACK","ZACAM","SHOWERS")),] # rich has some uncertainty in fire and sampling year could be added later; just removed because of uncertainty in fire and sampling year
d$Fire <- as.factor(d$Fire)

d$seedtr.comb <- ifelse((is.na(d$seedtr.dist) | (d$seedtr.dist >= 150)),d$dist.to.low,d$seedtr.dist)
d$seedtr.any.comb <- ifelse(is.na(d$seed.tree.any) | (d$seed.tree.any >= 150),d$dist.to.low,d$seed.tree.any)

d$regen.presab <- ifelse(d$regen.count > 0,TRUE,FALSE)

## calculate northness and eastness
d$northness <- cos(deg2rad(d$aspect))
d$eastness <- sin(deg2rad(d$aspect))

d$ppt.normal.sq <- d$ppt.normal^2
d$tmean.normal.sq <- d$tmean.normal^2


# (optional) thin to Sierra fires

sierra.fires <- c("STRAYLOR","CUB","RICH","DEEP","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER")
d$in.sierra <- ifelse(d$Fire %in% sierra.fires,TRUE,FALSE)
d <- d[d$in.sierra,]




#### Summary plots of climate distribution ####

#early vs. late normal precip: plots that were dry during the first three years were wet during the last two years
ggplot(d,aes(x=diff.norm.ppt.z,y=diff.norm.ppt.z.late,color=Fire)) +
  geom_point()

#early vs. late minimum precip: not very much correlation here
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.min.z.late,color=Fire)) +
  geom_point()

#early minimum precip vs. normal climate: #!may want to leave out Rich
ggplot(d,aes(x=diff.norm.ppt.min.z,y=ppt.normal,color=Fire)) +
  geom_point()  

#early minimum precip vs. early average precip
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.z,color=Fire)) +
  geom_point() 

#early minimum precip vs. late average precip: no correlation here but #! may want to leave out Rich
ggplot(d,aes(x=diff.norm.ppt.min.z,y=diff.norm.ppt.z.late,color=Fire)) +
  geom_point() 




#### run through for each tree species (and all tree species) and shrub, grass, forb
library(lme4)

## first with all species
d.sp <- d[d$species=="HARDWOOD",]



#rescale parameters
d.sp$regen.presab <- ifelse(d.sp$regen.count > 0,TRUE,FALSE)
d.sp$ppt.post.sq <- d.sp$ppt.post^2
d.sp$ppt.post.min.sq <- d.sp$ppt.post.min^2

leave.cols <- c("regen_presab","regen_count","Regen_Plot","Fire","Year","survey.years.post","SHRUB","FORB","GRASS")

d.sp.c <- center.df(d.sp,leave.cols=leave.cols)
stan_d1_mean <- mean.df(d.sp,leave.cols=leave.cols)
stan_d1_sd <- sd.df(d.sp,leave.cols=leave.cols)


d.sp.c.focal <- with(d.sp.c,data.frame(regen.presab,northness_c,ppt.normal_c,ppt.normal.sq_c,ppt.post_c,diff.norm.ppt_c,ppt.post.min_c,ppt.post.sq_c,ppt.post.min.sq_c,tmean.normal_c,tmean.normal.sq_c,diff.norm.ppt.z_c,diff.norm.ppt.min.z_c,Fire,seedtr.any.comb_c,eastness_c,diff.norm.tmean.JJA.mean.z_c,diff.norm.tmean.DJF.mean.z_c,diff.norm.tmean.z_c,seedtr.comb_c,SHRUB,FORB,GRASS,Year))
d.sp.c.focal <- d.sp.c.focal[complete.cases(d.sp.c.focal),]


# calculate what proportion of plots have the focal species
sum(d.sp.c.focal$regen.presab)/nrow(d.sp.c.focal)


# aggregate by fire and look for correlation
d.agg <- aggregate(d.sp.c.focal,by=list(d.sp.c.focal$Fire),FUN=mean)

d.year <- data.frame(Fire=d.sp$Fire,Year=d.sp$Year)
d.year <- unique(d.year)

d.agg <- merge(d.agg,d.year,by.x="Group.1",by.y="Fire")
d.agg$in.sierra <- ifelse(d.agg$Group.1 %in% sierra.fires,1,2)


plot(diff.norm.ppt.z_c~diff.norm.tmean.z_c,col=in.sierra,d=d.agg)
plot(diff.norm.ppt.min.z_c~Year,col=in.sierra,d=d.agg)



d.agg$fire.year <- paste(d.agg$Group.1,d.agg$Year.x-5)

library(ggplot2)
library(directlabels)
p <- ggplot(d.agg, aes(x=diff.norm.ppt.z_c,y=diff.norm.ppt.min.z_c,label=fire.year))+
  geom_point(size=5) +
  geom_text(check_overlap=TRUE) +
  scale_x_continuous(limits=c(-2.5,1.5))

p + direct.label(p, 'last.points')

m <- lm(diff.norm.ppt.z_c~diff.norm.ppt.min.z_c,d=d.agg)


# 
# ##### Attempt at using ppt anom--normal interaction #####
# ## 3-year average
# 
# m.full <- glmer(regen.presab ~ 
#                   #seedtr.any.comb_c +
#                   #seedtr.comb_c +
#                   northness_c +
#                   ppt.post_c * diff.norm.ppt_c +
#                   #ppt.post.sq_c * diff.norm.ppt.z_c +
#                   #ppt.normal_c * diff.norm.ppt.min.z_c +
#                   #ppt.normal.sq_c +
#                   #tmean.normal_c * diff.norm.tmean.z_c +
#                   #tmean.normal_c * diff.norm.tmean.JJA.mean.z_c +
#                   #tmean.normal_c * diff.norm.tmean.DJF.mean.z_c +
#                   #tmean.normal_c +
#                   #tmean.normal.sq_c +
#                   (1|Fire),
#                 data=d.sp.c.focal, family=binomial, na.action="na.fail")
# 
# ## 3-year minimum
# 
# m.full <- glmer(regen.presab ~ 
#                   #seedtr.any.comb_c +
#                   #seedtr.comb_c +
#                   northness_c +
#                   ppt.post.min_c * diff.norm.ppt.min.z_c +
#                   ppt.post.min.sq_c * diff.norm.ppt.min.z_c +
#                   #ppt.normal_c * diff.norm.ppt.min.z_c +
#                   #ppt.normal.sq_c +
#                   #tmean.normal_c * diff.norm.tmean.z_c +
#                   #tmean.normal_c * diff.norm.tmean.JJA.mean.z_c +
#                   #tmean.normal_c * diff.norm.tmean.DJF.mean.z_c +
#                   #tmean.normal_c +
#                   #tmean.normal.sq_c +
#                   (1|Fire),
#                 data=d.sp.c.focal, family=binomial, na.action="na.fail")
# 
# 
# 
# 





#******* average and extreme post-fire, WITH interactions (ppt only)


m.full <- glmer(regen.presab ~ 
              seedtr.any.comb_c +
              #seedtr.comb_c +
              #northness_c +
              ppt.normal_c * diff.norm.ppt.z_c +
              ppt.normal_c * diff.norm.ppt.min.z_c +
              ppt.normal.sq_c +
              #tmean.normal_c * diff.norm.tmean.z_c +
              #tmean.normal_c * diff.norm.tmean.JJA.mean.z_c +
              #tmean.normal_c * diff.norm.tmean.DJF.mean.z_c +
              #tmean.normal_c +
              #tmean.normal.sq_c +
              (1|Fire),
            data=d.sp.c.focal, family=binomial, na.action="na.fail")
summary(m.full)

library(MuMIn)
library(snow)

cl <- makeCluster(4, type = "SOCK")
clusterEvalQ(cl, library(lme4))
clusterExport(cl,"d.sp.c.focal")
a <- pdredge(m.full,cluster=cl, fixed=c("northness_c","ppt.normal_c"),
             subset=
               #(!(`seedtr.any.comb_c`&`seedtr.comb_c`)) &
               #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.JJA.mean.z_c`)) &
               #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.DJF.mean.z_c`)) &
               (!(`ppt.normal.sq_c` & !`ppt.normal_c`))
               #(!(`tmean.normal.sq_c` & !`tmean.normal_c`))
               )
stopCluster(cl)

imp <- data.frame(var=names(importance(a)),imp=importance(a))
best.model.num <- as.character(row.names(a[1,]))
best.model <- get.models(a,best.model.num)[[1]]
best.model.tree <- best.model
best.coef <- fixef(best.model)
best.coef <- data.frame(var=names(best.coef),coef=best.coef)
imp.coef <- merge(imp,best.coef,by="var",all.x=TRUE)

write.csv(imp.coef,"mumin_out.csv",row.names=FALSE)


#### "average" the model and make predictions ####

## iterate N times
## for each iteration, draw a model weighted based on its AICc weight
##    then draw one random set of parameter estimates from it using sim
##    compile into a prediction matrix



##################### Percent cover prediction ##################

library(arm)

hist(d.sp.c.focal$SHRUB,nclass=30)
hist(d.sp.c.focal$FORB)
hist(d.sp.c.focal$GRASS)

d.sp.c.focal$shrub.prop <- d.sp.c.focal$SHRUB/100
d.sp.c.focal$forb.prop <- d.sp.c.focal$FORB/100
d.sp.c.focal$grass.prop <- d.sp.c.focal$GRASS/100

d.sp.c.focal$shrub.prop <- ifelse(d.sp.c.focal$shrub.prop==0,0.005,d.sp.c.focal$shrub.prop)
d.sp.c.focal$shrub.prop <- ifelse(d.sp.c.focal$shrub.prop==1,0.995,d.sp.c.focal$shrub.prop)

d.sp.c.focal$grass.prop <- ifelse(d.sp.c.focal$grass.prop==0,0.005,d.sp.c.focal$grass.prop)
d.sp.c.focal$grass.prop <- ifelse(d.sp.c.focal$grass.prop==1,0.995,d.sp.c.focal$grass.prop)

d.sp.c.focal$SHRUB <- ifelse(d.sp.c.focal$SHRUB == 0,0.5,d.sp.c.focal$SHRUB)

d.sp.c.focal$GRASS <- ifelse(d.sp.c.focal$GRASS == 0,0.5,d.sp.c.focal$GRASS)

d.sp.c.focal$shrub.l <- logit(d.sp.c.focal$shrub.prop)
hist(d.sp.c.focal$shrub.l)
summary(d.sp.c.focal$shrub.l)

d.sp.c.focal$grass.l <- logit(d.sp.c.focal$grass.prop)



#******* average and extreme post-fire
library(glmmADMB)

m.full <-lmer(grass.l ~ 
                 northness_c +
                 ppt.normal_c * diff.norm.ppt.z_c +
                 ppt.normal_c * diff.norm.ppt.min.z_c +
                 ppt.normal.sq_c +
                 #tmean.normal_c +
                 #tmean.normal.sq_c +
                 diff.norm.ppt.z_c +
                 diff.norm.ppt.min.z_c +
                 #diff.norm.tmean.z_c +
                 #diff.norm.tmean.JJA.mean.z_c +
                 #diff.norm.tmean.DJF.mean.z_c +
                 #tmean.normal_c * diff.norm.tmean.z_c +
                 #tmean.normal_c * diff.norm.tmean.JJA.mean.z_c +
                 #tmean.normal_c * diff.norm.tmean.DJF.mean.z_c +

                 (1|Fire),
                data=d.sp.c.focal, na.action="na.fail", REML=FALSE)#,family=gaussian(link="logit"))
summary(m.full)
hist(resid(m.full))





#### Multi-model coefficient simulation ####

library(MuMIn)
library(snow)

cl <- makeCluster(4, type = "SOCK")
clusterEvalQ(cl, library(lme4))
clusterExport(cl,"d.sp.c.focal")
a <- pdredge(m.full,cluster=cl,
             subset=
               #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.JJA.mean.z_c`)) &
               #(!(`diff.norm.tmean.z_c` & `diff.norm.tmean.DJF.mean.z_c`)) &
               (!(`ppt.normal.sq_c` & !`ppt.normal_c`))
               #(!(`tmean.normal.sq_c` & !`tmean.normal_c`))
)
stopCluster(cl)


mod.names <-row.names(a)
mod.weights <- a[,"weight"]


# get all fit models
models <- list()

for(i in 1:length(mod.names)) {
  
  mod.name <- mod.names[i]
  models[[mod.name]] <- get.models(a,mod.name)[[1]]
  
}

library(arm)
library(plyr)
N <- 1000
nsims <- N
coef.sim <- data.frame()


for(i in 1:N) {
  
  mod.name <- sample(mod.names,prob=mod.weights,replace=TRUE,size=1)
  mod.name <- mod.names[1] # if using only the best model (otherwise comment out for multi-model averaging)
  model <- models[[mod.name]]
  coef.sim.row <- fixef(sim(model,n.sims=1)) # change to get just one simulation, store in matrix
  coef.sim.row <- as.data.frame(coef.sim.row)
  coef.sim <- rbind.fill(coef.sim,coef.sim.row)
  
}

coef.sim[is.na(coef.sim)] <- 0




################################
#### Make predictions ##########
################################

#link ppt.normal to ppt.sq

ppt.link <- data.frame(ppt.normal_c=d.sp.c.focal$ppt.normal_c,ppt.normal.sq_c=d.sp.c.focal$ppt.normal.sq_c)

m <- lm(ppt.normal.sq_c ~ ppt.normal_c + I(ppt.normal_c^2),ppt.link)
# -0.186 + 0.846x + 0.1863x^2


# prediction matrix (newdata)

z.seq <- seq(from=-2,to=2,length.out=100)

pred.df <- data.frame(
  
  diff.norm.ppt.min.z_c = z.seq,
  ppt.normal_c =0,
  ppt.normal.sq_c = -0.186 + 0.846*0 + 0.186*0^2,
  northness_c= 0,
  diff.norm.ppt.z_c = 0,
  seedtr.comb_c = 0,
  seedtr.any.comb_c = 0

)

# 
# #### alternative pred df for 3-d plotting: all potential combinations of min and average
# pairwise <- expand.grid(z.seq,z.seq)
# pred.df <- data.frame(
#   diff.norm.ppt.min.z_c = pairwise[,1],
#   ppt.normal_c = 0,
#   ppt.normal.sq_c = -0.186 + 0.846*0 + 0.186*0^2,
#   northness_c= 0,
#   diff.norm.ppt.z_c = pairwise[,2]
# )
# ### end alternative pred df ####


# interact.df <- matrix(data=NA,nrow=nsims,ncol=ncol(pred.df)^2)
# names.df <- vector(length=ncol(pred.df)^2)

interact.df <- data.frame("Fake"=rep(NA,nsims))


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
pred.df$"(Intercept)" <- 1


# ### predict for all tree seedlings
# library(merTools)
# 
# predFun <- function(x)predict(x, pred.df, re.form=NA)
# 
# p <- predict(best.model.tree,newdata=pred.df,type="response")
# p2 <- bootMer(best.model.tree,nsim=100,FUN=predFun)
# #predMat <- p2$t
# 
library(arm)

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
  pred.y.resp <- pred.y
  pred.y.resp <- invlogit(pred.y)
  #pred.y.resp <- exp(pred.y)
  pred.y.resp.mean <- mean(pred.y.resp)
  pred.y.resp.lwr <- quantile(pred.y.resp,probs=0.05)
  pred.y.resp.upr <- quantile(pred.y.resp,probs=0.95)
  
  ret <- data.frame(mean=pred.y.resp.mean,lwr=pred.y.resp.lwr,upr=pred.y.resp.upr)
  
  pred.fit[i,] <- ret

}

pred.fit <- pred.fit*100

pred.fit <- cbind(pred.df,pred.fit)


#plot 2-d
library(ggplot2)
p <- ggplot(pred.fit,aes(x=diff.norm.ppt.min.z_c,y=mean)) +
         geom_line(color="darkgreen") +
         geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
         scale_y_continuous(limits = c(-6,100)) +
         xlab("Postfire precipitation (3-yr minimum)") +
         ylab("Percent cover") +
         theme_bw(13)


library(Cairo)
Cairo(file="test2.png",typ="png",width=400,height=300)
p
dev.off()


#plot 3-d
library(ggplot2)
library(viridis)

p <- ggplot(pred.fit,aes(x=diff.norm.ppt.z_c,y=diff.norm.ppt.z_c)) +
  geom_tile(aes(fill=mean)) +
  #geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,fill="darkgreen") +
  #scale_y_continuous(limits = c(-6,100)) +
  #xlab("Postfire precipitation (3-yr average)") +
  #ylab("Percent cover") +
  theme_bw(13) +
  scale_fill_viridis(limits=c(0,100))


library(Cairo)
Cairo(file="test2.png",typ="png",width=400,height=300)
p
dev.off()





