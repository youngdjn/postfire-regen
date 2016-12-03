library(reshape)
library(plyr)


#### summarize post-fire climate ####
summarize.clim <- function(plot.df,plot.climate.df,years.clim) {

  plot.ids <- unique(plot.df$Regen_Plot)
  
  clim.plot.ret <- data.frame()
  for (plot.id in plot.ids) {
    
    plot.plot <- plot.df[plot.df$Regen_Plot == plot.id,]
    
    ## extract post-fire climate data
    clim.plot <- plot.climate.df[plot.climate.df$Regen_Plot == plot.id,]
    
    if(nrow(clim.plot) == 0) {next()} #no climate data for this plot (bad coords?); skip
    
    ##compute mean and SD of each variable
    
    tmean.DJF.cols <- grep("tmean.DJF.[0-9][0-9][0-9][0-9]",names(clim.plot))
    tmean.JJA.cols <- grep("tmean.JJA.[0-9][0-9][0-9][0-9]",names(clim.plot))
    tmean.cols <- grep("tmean.[0-9][0-9][0-9][0-9]",names(clim.plot))
    ppt.cols <- grep("ppt.[0-9][0-9][0-9][0-9]",names(clim.plot))
    
    tmean.DJF.plot <- unlist(clim.plot[,tmean.DJF.cols])
    tmean.JJA.plot <- unlist(clim.plot[,tmean.JJA.cols])
    tmean.plot <- unlist(clim.plot[,tmean.cols])
    ppt.plot <- unlist(clim.plot[,ppt.cols])
    
    tmean.DJF.mean <- mean(tmean.DJF.plot)
    tmean.JJA.mean <- mean(tmean.JJA.plot)
    tmean.mean <- mean(tmean.plot)
    ppt.mean <- mean(ppt.plot)
    
    tmean.DJF.sd <- sd(tmean.DJF.plot)
    tmean.JJA.sd <- sd(tmean.JJA.plot)
    tmean.sd <- sd(tmean.plot)
    ppt.sd <- sd(ppt.plot)

    
    ## get the post-fire values of each variable
    
    clim.years <- plot.plot$fire.year + years.clim
    
    tmean.cols <- paste("tmean.",clim.years,sep="")
    tmean.JJA.cols <- paste("tmean.JJA.",clim.years,sep="")
    tmean.DJF.cols <- paste("tmean.DJF.",clim.years,sep="")
    ppt.cols <- paste("ppt.",clim.years,sep="")
    
    snow.cols <- paste("snow.",clim.years,sep="")
    rain.cols <- paste("rain.",clim.years,sep="")
    
    tmean.plot <- unlist(clim.plot[,tmean.cols])
    tmean.DJF.plot <- unlist(clim.plot[,tmean.DJF.cols])
    tmean.JJA.plot <- unlist(clim.plot[,tmean.JJA.cols])
    ppt.plot <- unlist(clim.plot[,ppt.cols])
    snow.plot <- unlist(clim.plot[,snow.cols])
    rain.plot <- unlist(clim.plot[,rain.cols])    
    
    tmean.post.mean <- mean(tmean.plot)
    tmean.DJF.post.mean <- mean(tmean.DJF.plot)
    tmean.JJA.post.mean <- mean(tmean.JJA.plot)
    tmean.post.max <- max(tmean.plot)
    tmean.DJF.post.min <- min(tmean.DJF.plot)
    tmean.JJA.post.max <- max(tmean.JJA.plot)
    ppt.post.mean <- mean(ppt.plot)
    ppt.post.min <- min(ppt.plot)
    rain.post.min <- min(rain.plot)
    rain.post.mean <- mean(rain.plot)
    snow.post.min <- min(snow.plot)
    snow.post.max <- max(snow.plot)
    snow.post.mean <- mean(snow.plot)
    
    tmean.normal <- clim.plot$tmean.normal.ann
    ppt.normal <- clim.plot$ppt.normal
    snow.normal <- clim.plot$snow.normal
    rain.normal <- clim.plot$rain.normal
    
    perc.norm.ppt <- ppt.post.mean/ppt.normal
    perc.norm.ppt.min <- ppt.post.min/ppt.normal
    diff.norm.ppt <- ppt.post.mean-ppt.normal
    diff.norm.ppt.min <- ppt.post.min-ppt.normal
    
    diff.norm.ppt.z <- (ppt.post.mean-ppt.mean)/ppt.sd
    diff.norm.ppt.min.z <- (ppt.post.min-ppt.mean)/ppt.sd
    
    diff.norm.tmean.z <- (tmean.post.mean-tmean.mean)/ppt.sd
    diff.norm.tmean.max.z <- (tmean.post.max-tmean.mean)/tmean.sd
    diff.norm.tmean.JJA.mean.z <- (tmean.JJA.post.mean-tmean.JJA.mean)/tmean.JJA.sd
    diff.norm.tmean.DJF.mean.z <- (tmean.DJF.post.mean-tmean.DJF.mean)/tmean.DJF.sd
    diff.norm.tmean.JJA.max.z <- (tmean.JJA.post.max-tmean.JJA.mean)/tmean.JJA.sd
    diff.norm.tmean.DJF.min.z <- (tmean.DJF.post.min-tmean.DJF.mean)/tmean.DJF.sd
    
    
    ### compute total negative departure from normal
    # each years' departure
    ppt.anom <- ppt.plot - ppt.normal
    ppt.anom.neg <- ppt.anom[ppt.anom < 0]
    tot.neg.ppt.anom <- sum(ppt.anom.neg)/length(years.clim)
    
    
    
    diff.norm.tmean <- tmean.post.mean - tmean.normal
    
    clim.plot.out <- data.frame(Regen_Plot=plot.id,tmean.post=tmean.post.mean,ppt.post=ppt.post.mean,ppt.post.min=ppt.post.min,tmean.normal=tmean.normal,ppt.normal=ppt.normal,perc.norm.ppt=perc.norm.ppt,perc.norm.ppt.min=perc.norm.ppt.min,diff.norm.ppt,diff.norm.ppt.min,diff.norm.tmean=diff.norm.tmean,
                                diff.norm.ppt.z,diff.norm.ppt.min.z,diff.norm.tmean.z,diff.norm.tmean.max.z,
                                diff.norm.tmean.JJA.mean.z, diff.norm.tmean.DJF.mean.z,
                                diff.norm.tmean.JJA.max.z, diff.norm.tmean.DJF.min.z,
                                ppt.plot[1],ppt.plot[2],ppt.plot[3],ppt.plot[4],
                                rain.post.min,rain.post.mean,snow.post.min,snow.post.max,snow.post.mean,
                                snow.plot[1],snow.plot[2],snow.plot[3],snow.plot[4],
                                rain.plot[1],rain.plot[2],rain.plot[3],rain.plot[4],
                                snow.normal,rain.normal, tot.neg.ppt.anom)
    # removed these from above: tmean.plot,ppt.plot
  
    clim.plot.ret <- rbind.fill(clim.plot.ret,clim.plot.out)
  }
  return(clim.plot.ret)
}




#### summarize regen counts for given species and ages ####
summarize.regen <- function(plot.df,regen.df,sp,years.regen) {
  
  plot.ids <- unique(plot.df$Regen_Plot)
  
  regen.sp <- regen.df[regen.df$species %in% sp,]
  
  regen.ret <- data.frame()
  for (plot.id in plot.ids) {
    
    ## extract post-fire regen data
    regen.cols <- paste("count.",years.regen,"yr",sep="")
    regen.sp.plot <- regen.sp[regen.sp$Regen_Plot == plot.id,regen.cols]
    regen.sp.plot.tot <- sum(unlist(regen.sp.plot),na.rm=TRUE)
    
    seedtr.sp.plot <- regen.sp[regen.sp$Regen_Plot == plot.id,"seedtr.dist"]
    
    if(all(is.na(seedtr.sp.plot) | is.null(seedtr.sp.plot))) { #if they are all NULL and/or NA
      seedtr.sp.plot.min <- NA      
    } else {
      seedtr.sp.plot.min <- min(unlist(seedtr.sp.plot),na.rm=TRUE)
    }
    
    regen.sp.plot.out <- data.frame(Regen_Plot=plot.id,regen.count=regen.sp.plot.tot,seed.tree=seedtr.sp.plot.min)
    regen.ret <- rbind.fill(regen.ret,regen.sp.plot.out)
  }
  return(regen.ret)
}



#### summarize regen counts for given ages, with all specified species as separate rows ####
summarize.regen.ind <- function(plot.df,regen.df,sp,years.regen,all.sp=FALSE) {
  
  plot.ids <- unique(plot.df$Regen_Plot)
  
  if(all.sp) sp <- unique(regen.df$species) # override the species list provided; use all species that were observed at least once
  regen.sp <- regen.df[regen.df$species %in% sp,]

  
  regen.cols <- paste("count.",years.regen,"yr",sep="")
  regen.sp$regen.count <- apply(regen.sp[,regen.cols],1,sum,na.rm=TRUE)
  regen.sp <- regen.sp[,c("Regen_Plot","species","regen.count","seedtr.dist","surviving.trees")] # could add in seed tree here if desired
  
  regen.all <- merge(plot.ids,sp,by=NULL) # list of all species for all plots
  names(regen.all) <- c("Regen_Plot","species")
  regen.tot <- merge(regen.all,regen.sp,by=c("Regen_Plot","species"),all.x=TRUE)
  regen.tot$regen.count[is.na(regen.tot$regen.count)] <- 0
  regen.tot$surviving.trees[is.na(regen.tot$surviving.trees)] <- 0
  
  ##create data frame of surviving counts per species per plot
  species.all <- regen.all
  species.all$species <- paste(species.all$species,"surviving",sep=".")
  species.surviving <- paste(regen.sp$species,"surviving",sep=".")
  surviving.df <- data.frame(Regen_Plot=regen.sp$Regen_Plot,species=species.surviving,regen.count=regen.sp$surviving.trees,seedtr.dist=NA)
  surviving.df2 <- merge(species.all,surviving.df,by=c("Regen_Plot","species"),all.x=TRUE)
  surviving.df2$regen.count[is.na(surviving.df2$regen.count)] <- 0
    #! need to make this include all species even if zeros (repeat what did for regen.all above)
  
  
  
  ## compute total across all species in each plot
  regen.allsp.tot <- aggregate(regen.tot$regen.count,by=list(regen.tot$Regen_Plot),FUN=sum)
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.allsp.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="ALL",regen.count=regen.allsp.tot$regen.count,seedtr.dist = NA)
  
  surviving.allsp.tot <- aggregate(regen.tot$surviving.trees,by=list(regen.tot$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot) <- c("Regen_Plot","surviving.count")
  surviving.allsp.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot$Regen_Plot,species="ALL.surviving",regen.count=surviving.allsp.tot$surviving.count,seedtr.dist = NA)
  
  
  

  ## compute total for all conifers
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("ABCO","PIPO","PSME","CADE27","PIJE","QUKE","PILA","PIAT","ABIES","PICO","PINUS","PSMA","TAXUS","TOCA","ABMA"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.con.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="CONIFER",regen.count=regen.allsp.tot$regen.count,seedtr.dist = NA)
  
  surviving.allsp.tot <- aggregate(regen.tot.focal$surviving.trees,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot) <- c("Regen_Plot","surviving.count")
  surviving.con.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot$Regen_Plot,species="CONIFER.surviving",regen.count=surviving.allsp.tot$surviving.count,seedtr.dist = NA)
  
  
  
  
  ## compute total for all hardwoods
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.hdwd.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="HARDWOOD",regen.count=regen.allsp.tot$regen.count,seedtr.dist = NA)

  surviving.allsp.tot <- aggregate(regen.tot.focal$surviving.trees,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot) <- c("Regen_Plot","surviving.count")
  surviving.hdwd.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot$Regen_Plot,species="HARDWOOD.surviving",regen.count=surviving.allsp.tot$surviving.count,seedtr.dist = NA)
  
  regen.tot <- regen.tot[,c("Regen_Plot","species","regen.count","seedtr.dist")] # remove the "surviving trees" column
  
  regen.tot <- rbind(regen.tot,regen.allsp.tot.df,regen.con.tot.df,regen.hdwd.tot.df,surviving.allsp.tot.df,surviving.con.tot.df,surviving.hdwd.tot.df,surviving.df2)
  
  regen.tot <- regen.tot[order(regen.tot$Regen_Plot),]
  

  return(regen.tot)
}



center.df <- function(df,leave.cols) {
  df.names <- names(df)
  new.df <- data.frame(SID=1:nrow(df))
  # get first species ID for thinning
  first.sp.id <- unique(df$species)[1]
  for(var.name in df.names) {
    col.vals <- df[,var.name]
    colnum <- which(df.names == var.name)
    if(!is.numeric(col.vals) | var.name %in% leave.cols) {
      new.df[,var.name] <- col.vals
    } else {
      new.var.name <- paste(var.name,"_c",sep="")
      col.vals.thinned <- col.vals[df$species == first.sp.id]
      new.df[,new.var.name] <- (col.vals - mean(col.vals.thinned,na.rm=TRUE)) / sd(col.vals.thinned,na.rm=TRUE)
    }
  }
  return(new.df)
}

mean.df <- function(df,leave.cols) {
  df.names <- names(df)
  means.array <- array()
  # get first species ID for thinning
  first.sp.id <- unique(df$species)[1]
  for(var.name in df.names) {
    col.vals <- df[,var.name]
    colnum <- which(df.names == var.name)
    if(is.numeric(col.vals) & !(var.name %in% leave.cols)) {
      col.vals.thinned <- col.vals[df$species == first.sp.id]
      means.array[var.name] <- mean(col.vals.thinned,na.rm=TRUE)
    }
  }
  return(means.array)
}

sd.df <- function(df,leave.cols) {
  df.names <- names(df)
  sds.array <- array()
  # get first species ID for thinning
  first.sp.id <- unique(df$species)[1]
  for(var.name in df.names) {
    col.vals <- df[,var.name]
    colnum <- which(df.names == var.name)
    if(is.numeric(col.vals) & !(var.name %in% leave.cols)) {
      col.vals.thinned <- col.vals[df$species == first.sp.id]
      sds.array[var.name] <- sd(col.vals.thinned,na.rm=TRUE)
    }
  }
  return(sds.array)
}

uncenter <- function(vect,mean,sd) {
  return((vect*sd)+mean)
}

recenter <- function (vect,mean,sd) {
  return((vect-mean)/sd)
}


deg2rad <- function(deg) {
  
  return(deg * (3.141593/180))

}
