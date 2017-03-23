extract.single <- function(layer, plots) {
  vals <- extract(layer,plots,method="bilinear")
  return(vals)
}

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
    def.cols <- grep("def.[0-9][0-9][0-9][0-9]",names(clim.plot))
    aet.cols <- grep("aet.[0-9][0-9][0-9][0-9]",names(clim.plot))
    
    tmean.DJF.plot <- unlist(clim.plot[,tmean.DJF.cols])
    tmean.JJA.plot <- unlist(clim.plot[,tmean.JJA.cols])
    tmean.plot <- unlist(clim.plot[,tmean.cols])
    ppt.plot <- unlist(clim.plot[,ppt.cols])
    def.plot <- unlist(clim.plot[,def.cols])
    aet.plot <- unlist(clim.plot[,aet.cols])
    
    tmean.DJF.mean <- mean(tmean.DJF.plot)
    tmean.JJA.mean <- mean(tmean.JJA.plot)
    tmean.mean <- mean(tmean.plot)
    ppt.mean <- mean(ppt.plot)
    def.mean <- mean(def.plot)
    aet.mean <- mean(aet.plot)
    
    tmean.DJF.sd <- sd(tmean.DJF.plot)
    tmean.JJA.sd <- sd(tmean.JJA.plot)
    tmean.sd <- sd(tmean.plot)
    ppt.sd <- sd(ppt.plot)
    def.sd <- sd(def.plot)
    aet.sd <- sd(aet.plot)
    
    
    ## get the post-fire values of each variable
    
    clim.years <- plot.plot$fire.year + years.clim
    
    tmean.cols <- paste("tmean.",clim.years,sep="")
    tmean.JJA.cols <- paste("tmean.JJA.",clim.years,sep="")
    tmean.DJF.cols <- paste("tmean.DJF.",clim.years,sep="")
    ppt.cols <- paste("ppt.",clim.years,sep="")
    def.cols <- paste("def.",clim.years,sep="")
    aet.cols <- paste("aet.",clim.years,sep="")
    
    snow.cols <- paste("snow.",clim.years,sep="")
    rain.cols <- paste("rain.",clim.years,sep="")
    
    tmean.plot <- unlist(clim.plot[,tmean.cols])
    tmean.DJF.plot <- unlist(clim.plot[,tmean.DJF.cols])
    tmean.JJA.plot <- unlist(clim.plot[,tmean.JJA.cols])
    ppt.plot <- unlist(clim.plot[,ppt.cols])
    snow.plot <- unlist(clim.plot[,snow.cols])
    rain.plot <- unlist(clim.plot[,rain.cols])
    def.plot <- unlist(clim.plot[,def.cols])
    aet.plot <- unlist(clim.plot[,aet.cols])
    
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
    
    aet.post.min <- min(aet.plot)
    aet.post.mean <- mean(aet.plot)
    def.post.max <- max(def.plot)
    def.post.mean <- mean(def.plot)
    
    tmean.normal <- clim.plot$tmean.normal.ann
    ppt.normal <- clim.plot$ppt.normal
    snow.normal <- clim.plot$snow.normal
    rain.normal <- clim.plot$rain.normal
    
    aet.normal <- aet.mean
    def.normal <- def.mean
    
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
    
    diff.norm.def.z <- (def.post.mean-def.mean)/def.sd
    diff.norm.def.max.z <- (def.post.max-def.mean)/def.sd
    
    diff.norm.aet.z <- (aet.post.mean-aet.mean)/aet.sd
    diff.norm.aet.min.z <- (aet.post.min-aet.mean)/aet.sd
    
    
    ### compute average negative (only) precipitation z-score departure
    # each years' departure
    ppt.anom.neg <- diff.norm.ppt.z[diff.norm.ppt.z < 0]
    tot.neg.ppt.anom <- sum(ppt.anom.neg)/length(years.clim)
    
    
    
    diff.norm.tmean <- tmean.post.mean - tmean.normal
    
    clim.plot.out <- data.frame(Regen_Plot=plot.id,tmean.post=tmean.post.mean,ppt.post=ppt.post.mean,def.post=def.post.mean,aet.post=aet.post.mean,ppt.post.min=ppt.post.min,def.post.max,aet.post.min,tmean.normal=tmean.normal,ppt.normal=ppt.normal,def.normal,aet.normal,perc.norm.ppt=perc.norm.ppt,perc.norm.ppt.min=perc.norm.ppt.min,diff.norm.ppt,diff.norm.ppt.min,diff.norm.tmean=diff.norm.tmean,
                                diff.norm.ppt.z,diff.norm.ppt.min.z,diff.norm.tmean.z,diff.norm.tmean.max.z,
                                diff.norm.def.z,diff.norm.def.max.z,diff.norm.aet.z,diff.norm.aet.min.z,
                                diff.norm.tmean.JJA.mean.z, diff.norm.tmean.DJF.mean.z,
                                diff.norm.tmean.JJA.max.z, diff.norm.tmean.DJF.min.z,
                                ppt.plot[1],ppt.plot[2],ppt.plot[3],ppt.plot[4],
                                def.plot[1],def.plot[2],def.plot[3],def.plot[4],
                                aet.plot[1],aet.plot[2],aet.plot[3],aet.plot[4],
                                rain.post.min,rain.post.mean,snow.post.min,snow.post.max,snow.post.mean,
                                snow.plot[1],snow.plot[2],snow.plot[3],snow.plot[4],
                                rain.plot[1],rain.plot[2],rain.plot[3],rain.plot[4],
                                snow.normal,rain.normal, tot.neg.ppt.anom)
    # removed these from above: tmean.plot,ppt.plot
    
    clim.plot.ret <- rbind.fill(clim.plot.ret,clim.plot.out)
  }
  return(clim.plot.ret)
}





#### summarize regen counts for given ages, with all specified species as separate rows ####
summarize.regen.ind <- function(plot.df,regen.df,sp,regen.ages,all.sp=FALSE,incl.unk.age.for.all=FALSE) {
  
  plot.ids <- unique(plot.df$Regen_Plot)

  if(all.sp) sp <- unique(regen.df$species) # if specified, override the species list provided; use all species that were observed at least once
  regen.sp <- regen.df[regen.df$species %in% sp,]
  
  
  ## function to determine which years of regen to use, based on plot age (number of years post-fire it was surveyed)
  ## then take the regen ages required for each plot and sum the counts of those ages in each plot

  get.regen.years <- function(plot.survey.years.post) {
    
    if(regen.ages == "young") {
      years.regen <- 1:2
    } else if(regen.ages == "old") {
      years.regen <- plot.survey.years.post - 1:-1 # the third-to-last and last year of the plot before it was surveyed (also the year after it was surveyed, in case of mis-aging)
    } else if(regen.ages == "all") {
      years.regen <- 1:(plot.survey.years.post+1) #+1 in case a seedling was over-aged by one year
    }

  }
  
  years.regen <- lapply(plot.df$survey.years.post,get.regen.years)
  nyears.regen <- lapply(years.regen,length)
  
  plot.df$years.regen <- years.regen
  
  regen.cols <- lapply(years.regen,function(x) paste("count.",x,"yr",sep=""))
  names(regen.cols) <- plot.df$Regen_Plot
  
  #for each row of regen.sp (plot-species combination), sum the counts of seedlings of the specified ages, and return the average number of seedlings per year
  regen.peryr.sp <- rep(NA,nrow(regen.sp))
  for(i in 1:nrow(regen.sp)) {
    
    regen.sp.row <- regen.sp[i,]
    
    regen.plot <- regen.sp.row["Regen_Plot"]
    
    # if there are no counts available for the plot, skip the plot-species combination
    if(!(regen.plot %in% names(regen.cols))) {
      #cat("No regen counts found for plot",regen.plot[[1]]," species", regen.sp.row$species)
      next()
    }
    regen.row.cols <- eval(parse(text=paste0("regen.cols$",regen.plot))) # get the names of the regen columns of the correct age seedlings for the current regen plot-species combination
    
    #include unknown age if examining all-aged trees
    if((regen.ages == "all")&(incl.unk.age.for.all==TRUE)) {
      regen.row.cols.complete <- c(regen.row.cols,"unk_yr") 
    } else {
      regen.row.cols.complete <- regen.row.cols
    }
    
    regen.tot.sp <- sum(regen.sp.row[,regen.row.cols.complete],na.rm=TRUE)
    
    #compute number of years
    if(regen.ages=="young") {
      regen.nyears <- length(regen.row.cols)
    } else if(regen.ages=="old") {
      regen.nyears <- length(regen.row.cols) -1 #-1 because added an extra year in case seedlings were over-aged, but in realiy they only came from the the years between the fire and the survey
    } else if(regen.ages=="all") {
      regen.nyears <- length(regen.row.cols) -1 #-1 because added an extra year in case seedlings were over-aged, but in realiy they only came from the the years between the fire and the survey
    }
      
    
    
    regen.peryr.sp[i] <- regen.tot.sp / regen.nyears #this is not "regen.row.cols.complete" because don't want to consider "unk_yr" as another year to have to divide counts by to get seedlings/yr (when computing all ages, the ages are restricted to the ages assigned to "all", which are all years between the fire and the survey)
    
  }
  
  regen.sp$regen.count <- regen.peryr.sp
  regen.sp <- regen.sp[,c("Regen_Plot","species","regen.count","surviving.trees.count","surviving.trees.ba")] 
  
  
  regen.all <- merge(plot.ids,sp,by=NULL) # list of all species for all plots
  names(regen.all) <- c("Regen_Plot","species")
  regen.tot <- merge(regen.all,regen.sp,by=c("Regen_Plot","species"),all.x=TRUE)
  regen.tot$regen.count[is.na(regen.tot$regen.count)] <- 0
  regen.tot$surviving.trees.count[is.na(regen.tot$surviving.trees.count)] <- 0
  regen.tot$surviving.trees.ba[is.na(regen.tot$surviving.trees.ba)] <- 0
  
  ##create data frame of surviving counts per species per plot (with all species--so 0 if none, rather than having no row)
  species.all <- regen.all #just a list of all species crossed with all plots
  surviving.df <- data.frame(Regen_Plot=regen.sp$Regen_Plot,species=regen.sp$species,adult.count=regen.sp$surviving.trees.count,adult.ba=regen.sp$surviving.trees.ba)
  surviving.df2 <- merge(species.all,surviving.df,by=c("Regen_Plot","species"),all.x=TRUE)
  surviving.df2$adult.count[is.na(surviving.df2$adult.count)] <- 0
  surviving.df2$adult.ba[is.na(surviving.df2$adult.ba)] <- 0
  
  
  ## compute total across all species in each plot
  regen.allsp.tot <- aggregate(regen.tot$regen.count,by=list(regen.tot$Regen_Plot),FUN=sum)
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.allsp.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="ALL",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot$surviving.trees.count,by=list(regen.tot$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot$surviving.trees.ba,by=list(regen.tot$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.allsp.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="ALL",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  
  ## compute total for all conifers
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("ABCO","PIPO","PSME","CADE27","PIJE","PILA","PIAT","ABIES","PICO","PINUS","PSMA","TAXUS","TOCA","ABMA","JUNIPERUS","JUOC","JUCA7","PISA","CONIFER"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)

  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.con.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="CONIF.ALLSP",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot.focal$surviving.trees.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot.focal$surviving.trees.ba,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.con.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="CONIF.ALLSP",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  ## compute total for all Pinus
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("PIPO","PIJE","PILA","PIAT","PICO","PINUS","PISA"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.pinus.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="PINUS.ALLSP",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot.focal$surviving.trees.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot.focal$surviving.trees.ba,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.pinus.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="PINUS.ALLSP",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  ## compute total for all yellow pine
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("PIPO","PIJE"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.pipj.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="PIPJ",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot.focal$surviving.trees.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot.focal$surviving.trees.ba,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.pipj.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="PIPJ",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  
  
  
  ## compute total for all shade sp
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("ABCO","CADE27","ABIES","TAXUS","TOCA","ABMA"),]
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.shade.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="SHADE.ALLSP",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot.focal$surviving.trees.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot.focal$surviving.trees.ba,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.shade.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="SHADE.ALLSP",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  
  
  
  ## compute total for all hardwoods
  regen.tot.focal <- regen.tot[regen.tot$species %in% c("ALNUS","QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","ACME","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA","ALRH"),] #! need to add some such as ALNUS
  regen.allsp.tot <- aggregate(regen.tot.focal$regen.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  
  names(regen.allsp.tot) <- c("Regen_Plot","regen.count")
  regen.hdwd.tot.df <- data.frame(Regen_Plot=regen.allsp.tot$Regen_Plot,species="HDWD.ALLSP",regen.count=regen.allsp.tot$regen.count)
  
  surviving.allsp.tot.count <- aggregate(regen.tot.focal$surviving.trees.count,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  surviving.allsp.tot.ba <- aggregate(regen.tot.focal$surviving.trees.ba,by=list(regen.tot.focal$Regen_Plot),FUN=sum)
  names(surviving.allsp.tot.count) <- c("Regen_Plot","adult.count")
  names(surviving.allsp.tot.ba) <- c("Regen_Plot","adult.ba")
  surviving.hdwd.tot.df <- data.frame(Regen_Plot=surviving.allsp.tot.count$Regen_Plot,species="HDWD.ALLSP",adult.count=surviving.allsp.tot.count$adult.count,adult.ba = surviving.allsp.tot.ba$adult.ba)
  
  
  regen.spgrps <- rbind.fill(regen.tot[,1:3],regen.allsp.tot.df,regen.con.tot.df,regen.hdwd.tot.df,regen.pinus.tot.df,regen.shade.tot.df,regen.pipj.tot.df)
  adult.spgrps <- rbind.fill(surviving.df2,surviving.allsp.tot.df,surviving.con.tot.df,surviving.hdwd.tot.df,surviving.pinus.tot.df,surviving.shade.tot.df,surviving.pipj.tot.df)
  
  trees.spgrps <- merge(regen.spgrps,adult.spgrps)

  trees.spgrps <- trees.spgrps[order(trees.spgrps$Regen_Plot),]
  
  
  return(trees.spgrps)
}


# remove specified columns from a data frame
remove.vars <- function(df,vars.remove) {
  
  df.names <- names(df)
  df.names.clean <- df.names[!(df.names %in% vars.remove)]
  df.clean <- df[,df.names.clean]
  
  return(df.clean)
  
}

