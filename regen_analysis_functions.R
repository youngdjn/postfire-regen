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
      m.pred <- one.t
    }
    
    err <- (logit(m.pred)-logit(obs))
    
    errs[i] <- err
    
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
