# Supporting functions for the analysis script (Young_et_al_analysis.R) for "Post-fire forest regeneration shows limited climate tracking and potential for drought-induced type conversion"
# By Derek J. N. Young, Chhaya M. Werner, Kevin R. Welch, Truman P. Young, Hugh D. Safford, Andrew M. Latimer

# test if a number is odd
is.odd <- function(x) x %% 2 != 0

# summarize a vector of numbers
summary.string <- function(x,decimals) {
  
  med <- round(median(x),digits=decimals)
  
  lims <- round(quantile(x,probs=c(0.0,1)),digits=decimals)
  range <- paste(lims,collapse=", ")
  sum.str <- paste0(med," (",range,")")
  
  return(sum.str)
}


# capitalize the first letter
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}



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



# Fit a model and perform fire-level cross-validation (withholding an entire fire at a time)
cvfun.fire <- function(formula,data) {
  
  by.fire <- TRUE
  
  data <- data[complete.cases(data$response.var),]

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
    
    if(sp %in% c(cover.opts)) { # 
      
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

      
    }  else { # it is not a cover or height dominance response
      
      m <- glm(formula,data=data.train,family="binomial")
      
      if(m$converged == FALSE) {
        message <- paste0("No convergence during CV for ",sp," ",formula.name,"\n")
        print(message)
        next()
      }
      
      m.fit <- predict(m,type="response")
      
      
    
      
      ## find the cutoff to use for presence/absence, based on the proportion of plots that had presence in the training dataset
      prop.present <- mean(data.train$response.var,na.rm=TRUE)

      m.pred <- predict(m,newdat=data.val,type="response")
      obs <- data.val$response.var
      
      pred <- mean(m.pred)
      obs <- mean(obs)
      
      err <- abs(pred-obs)     
      
      
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





### Functions to find the center (mean), lower (20th percentile) and upper (80th percentile) values of a "fixed" predictor variable within the species distribution
mid.val.fun <- function(var) {
  mid.val <- mean(d.c.modfit[,var],na.rm=TRUE)
  return(mid.val)
}

low.val.fun <- function(var) {
  a <- quantile(d.c.modfit[,var],0.20,na.rm=TRUE)
  #a <- min(d.c[,var])
  return(a)
}

high.val.fun <- function(var) {
  a <- quantile(d.c.modfit[,var],0.80,na.rm=TRUE)
  #a <- max(d.c[,var])
  return(a)
}



## Center all columns in a data frame except the specified columns (leave.cols)
center.df <- function(df,leave.cols) {
  df.names <- names(df)
  center.data <- data.frame()
  new.df <- data.frame(SID=1:nrow(df))
  
  for(var.name in df.names) {
    col.vals <- df[,var.name]
    colnum <- which(df.names == var.name)
    if(!is.numeric(col.vals) | var.name %in% leave.cols) { # if it's not a column to be centered
      new.df[,var.name] <- col.vals
    } else { # it is a column to be centered
      new.var.name <- paste(var.name,"_c",sep="")
      col.vals.thinned <- col.vals   # <- col.vals[df$species == first.sp.id]
      new.df[,new.var.name] <- (col.vals - mean(col.vals.thinned,na.rm=TRUE)) / sd(col.vals.thinned,na.rm=TRUE)
      
      
      
      ## store a DF with centering data
      var.mean <- mean(col.vals.thinned,na.rm=TRUE)
      var.sd <- sd(col.vals.thinned,na.rm=TRUE)
      var.center.dat <- data.frame(var=new.var.name,var.mean,var.sd)
      center.data <- rbind(center.data,var.center.dat)
      
    }
  }
  return(list(centered.df=new.df,center.data=center.data))
}



# Inverse logit transformation
inv.logit <- function(x) {
  exp(x)/(1+exp(x))
}





