## Given a variable (vector of values) and a set of breakpoints, determine between which breakpoints each value of the variable falls
categorize <- function(x,breaks,name) { #name is used as the prefix for all categories output for this variable
  var.bounds <- c(min(x,na.rm=TRUE)-1,breaks,max(x,na.rm=TRUE)+1)
  var.bounds.shift <- c(var.bounds[-1],NA)
  var.bounds.pairs <- rbind(var.bounds,var.bounds.shift)[,-length(var.bounds)]
  colnames(var.bounds.pairs) <- 1:ncol(var.bounds.pairs)
  category <- rep(NA,length(x))
  for(i in 1:ncol(var.bounds.pairs)) { # for each set of precip bounds
    bounds <- var.bounds.pairs[,i]
    category <- ifelse((x > bounds[1]) & (x <= bounds[2]),i,category)
  }
  category <- paste(name,category,sep=".")
  return(category)
}

# remove specified columns from a data frame
remove.vars <- function(df,vars.remove) {
  
  df.names <- names(df)
  df.names.clean <- df.names[!(df.names %in% vars.remove)]
  df.clean <- df[,df.names.clean]
  
  return(df.clean)
  
}




cvfun <- function(formula,data) {
  
  errs <- NULL
  
  for(i in 1:nrow(data)) {
    
    data.train <- data[-i,]
    data.val <- data[i,]
    
    m <- betareg(formula,data=data.train)
    
    if(m$converged == FALSE) {
      message <- paste0("No convergence for ",sp," ",formula.name)
      print(message)
    }
    
    
    m.pred <- predict(m,newdat=data.val)
    obs <- data.val$response.var
    
    if((m.pred > 0.999)) { # if predicted probability identical to 1, rescale to the largest possible observed value (transformed)
      n.pts <- nrow(data.train)+nrow(data.val)
      one.t <- (1*(n.pts-1) + 0.5) / n.pts
      m.pred <- one.t
    }
    
    if((m.pred < 0.001)) { # if predicted probability identical to 0, rescale to the largest possible observed value (transformed)
      n.pts <- nrow(data.train)+nrow(data.val)
      zero.t <- (0*(n.pts-1) + 0.5) / n.pts
      m.pred <- zero.t
    }
    
    err <- (logit(m.pred)-logit(obs))
    
    errs[i] <- err
    
  }
  
  
  
  
  mae <- mean(abs(errs))
  stderr <- sd(abs(errs))/(sqrt(length(errs)))
  
  ret <- data.frame(mae,stderr)
  return(ret)
  
}







#plot-level cross-validation (withholding an entire fire at a time)
cvfun.fire <- function(formula,data) {
  
  by.fire <- TRUE
  
  data <- data[!is.na(data$response.var),]
  
  
  
  if(by.fire == TRUE) {
    for.in <- 1:length(unique(data$Fire))
  } else {
    for.in <- 1:nrow(data)
  }
  
  errs <- NULL
  
  for(i in for.in) {
    
  
    if(by.fire == TRUE) {
      Fire <- unique(data$Fire)[i]
      data.train <- data[data$Fire != Fire,]
      data.val <- data[data$Fire == Fire,]
    } else {
      data.train <- data[-i,]
      data.val <- data[i,]
    }
    
    if(sp %in% c(cover.opts,prop.opts)) {
      
      data.train <- data.train[!is.na(data.train$response.var),] 
      
      
      m <- try(betareg(formula,data=data.train),silent=TRUE)
      if(class(m) == "try-error") {
        message <- paste0("\nModel error during CV for ",sp," ",formula.name,"\n")
        print(message)
        next()
      }
      
      if(m$converged == FALSE) {
        message <- paste0("No convergence during CV for ",sp," ",formula.name,"\n")
        print(message)
        next()
      }
      
      m.pred <- predict(m,newdat=data.val)
      obs <- data.val$response.var
      
      # if predicted probability identical to 1, rescale to the largest possible observed value (transformed)
        n.pts <- nrow(data.train)+nrow(data.val)
        one.t <- (1*(n.pts-1) + 0.5) / n.pts
        m.pred[m.pred > 0.999] <- one.t
      
      # if predicted probability identical to 0, rescale to the largest possible observed value (transformed)
        n.pts <- nrow(data.train)+nrow(data.val)
        zero.t <- (0*(n.pts-1) + 0.5) / n.pts
        m.pred[m.pred < 0.001] <- zero.t
      
      
      err <- mean(abs(logit(m.pred)-logit(obs)),na.rm=TRUE)

      
    } else if (sp %in% htabs.opts) { # we're looking at absolute height
      
      
      m <- glm(formula,data=data.train,family="gaussian")
      
      if(m$converged == FALSE) {
        message <- paste0("No convergence during CV for ",sp," ",formula.name,"\n")
        print(message)
        next()
      }
      
      m.pred <- predict(m,newdat=data.val)
      obs <- data.val$response.var
      
      err <- mean(abs(m.pred-obs),na.rm=TRUE)  #! if we end up using a transformation like log we'll need to log these too
      
      
      
    } else { # it is not a cover or height response
      
      m <- glm(formula,data=data.train,family="binomial")
      
      if(m$converged == FALSE) {
        message <- paste0("No convergence during CV for ",sp," ",formula.name,"\n")
        print(message)
        next()
      }
      
      m.fit <- predict(m,type="response")
      
      ## find the cutoff to use for presence/absence, based on the proportion of plots that had presence in the training dataset
      prop.present <- mean(data.train$response.var,na.rm=TRUE)
      
      #sort plots by predicted probability
      m.fit.sort <- sort(m.fit,decreasing=TRUE)
      
      #what number of them should be presences?
      num.present <- round(prop.present*sum(!is.na(data.train$response.var)))
      
      #what is the probability at that point? use it as the cutoff
      cutoff <- m.fit.sort[num.present]
      
      # comp <- getCompSpecAndSens(data.train$response.var,m.fit)
      comp <- cutoff
      
      m.pred <- predict(m,newdat=data.val,type="response")
      m.pred.presab <- ifelse(m.pred > comp,1,0)
      
      obs <- data.val$response.var
      
      err <- 1-mean(m.pred.presab == obs,na.rm=TRUE)
      
      
    }
    


    
    
    errs[i] <- err
    
  }
  
  
  
  if(length(errs) < 3) {  #super high number if not enought of the fires converged
    errs <- c(9999,9999,9999)
  }
  
  mae <- mean(abs(errs))
  stderr <- sd(abs(errs))/(sqrt(length(errs)))
  
  ret <- data.frame(mae,stderr)
  return(ret)
  
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

inv.logit <- function(x) {
  exp(x)/(1+exp(x))
}






getCompSpecAndSens <- function(obs, fit)
{
  ascend <- order(fit)
  
  hosmer <- matrix(c(obs[ascend], fit[ascend]), ncol=2)
  #plot(hosmer, pch=16, col="#00000010")
  
  obs.pos <- which(hosmer[,1] == 1)
  obs.neg <- which(hosmer[,1] == 0)
  
  cutoffs <- seq(0, 1, .01)
  sens <- spec <- rep(NA, length(cutoffs))
  
  for(i in 1:length(cutoffs))
  {
    pred.pos <- which(hosmer[,2] >= cutoffs[i])
    pred.neg <- which(hosmer[,2] < cutoffs[i])
    
    oppp <- length(which(obs.pos %in% pred.pos))
    onpp <- length(which(obs.neg %in% pred.pos))
    oppn <- length(which(obs.pos %in% pred.neg))
    onpn <- length(which(obs.neg %in% pred.neg))
    
    sens[i] <- oppp / (oppp + oppn)
    spec[i] <- onpn / (onpp + onpn)
    
    # cat(paste("\ncutoff: ", cutoffs[i], "\n", sep=""))
    # cat(paste("\t\top\ton\n\tpp\t", oppp, "\t", onpp, "\n\tpn\t", oppn, "\t", onpn, "\n", sep=""))
    # cat(paste("\tsensitivity: ", sens[i], "\n\tspecificity: ", spec[i], "\n", sep=""))
  }
  
  compens <- which(spec > sens)[1]
  comp.cutoff <- (cutoffs[compens - 1] + cutoffs[compens]) / 2
  comp.spec <- (spec[compens - 1] + spec[compens]) / 2
  comp.sens <- (sens[compens - 1] + sens[compens]) / 2
  #return(list(comp.cutoff=comp.cutoff, comp.spec=comp.spec, comp.sens=comp.sens))
  return(comp.cutoff)
}






