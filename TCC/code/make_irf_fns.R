##########################################################################################
#
# make_irf_fn.R
#
# Code to define the functions that make the irfs for the "Inflation Sentiments" project
# Philip Barrett, Washington DC
# First version: 28jan2022
#
##########################################################################################

make.irf.rf <- function( l.var, n.pds=20 ){
  # Makes the reduced form IRFs for the var and returns them as a dataframe
  
  ## Calculations
  n.vars <- l.var$n.vars # var.names %>% length()
  # The number of variables
  var.names <- l.var$B %>% rownames %>% .[1:n.vars]
  # The variable names
  l.irfs <- lapply( 1:n.pds, function(x) out <- l.var$B %^% x )
  # The list of IRF vectors
  
  ## Output: wrangling the initial period
  m.init <- diag(n.vars)
  rownames(m.init) <- colnames(m.init) <- var.names
  # The initial period response
  df.init <- m.init %>%
    as.data.frame() %>%
    rownames_to_column('outcome') %>%
    gather( shock, value, -outcome) %>%
    mutate( period=0 ) %>%
    select( shock, outcome, period, value )
  
  ## Output: wrangling the later periods
  df.out <- df.init
  # Initialization
  for(i in 1:n.pds){
    m.add <- l.irfs[[i]][1:n.vars,1:n.vars]
    rownames(m.add) <- colnames(m.add) <- var.names
    # The response
    df.add <- m.add %>%
      as.data.frame() %>%
      rownames_to_column('outcome') %>%
      gather( shock, value, -outcome) %>%
      mutate( period=i ) %>%
      select( shock, outcome, period, value )
    # In data frame form
    df.out <- df.out %>%
      rbind( df.add )
      # full_join( df.add, by = c("shock", "outcome", "period", "value") )
  }
  return( df.out %>% arrange(shock, outcome, period) )
}

make.irf.struct <- function( l.var, A, fcast.horiz, inf.name, n.pds=20, fcast=fcast,
                             fcast.cumul=TRUE ){
# Makes the structural IRFs using the shock decomp in A
  
  ## Calculations
  n.vars <- l.var$n.vars # var.names %>% length()
  # The number of variables
  p <- l.var$lags
  # The number of lags
  var.names <- l.var$B %>% rownames %>% .[1:n.vars]
  # The variable names
  A.extended <- 0*diag( n.vars * p ) ; A.extended[1:n.vars, 1:n.vars] <- A
  # The extended version of A
  l.irfs <- lapply( 1:n.pds, function(x) out <- (l.var$B %^% x) %*% A.extended )
  # The list of IRF vectors
  
  ## Output: wrangling the initial period
  m.init <- A
  rownames(m.init) <- var.names
  n.fcast <- inf.name %>% length()
  shk.names <- c( paste0( 'Sentiment #', 1:n.fcast), paste0( 'Fundamental #', 1:(n.vars-n.fcast) ) )
  colnames(m.init) <- shk.names
    # The initial period response
  df.init <- m.init %>%
    as.data.frame() %>%
    rownames_to_column('outcome') %>%
    gather( shock, value, -outcome) %>%
    mutate( period=0 ) %>%
    select( shock, outcome, period, value )
  
  ## Output: wrangling the later periods
  df.out <- df.init
    # Initialization
  for(i in 1:n.pds){
    m.add <- l.irfs[[i]][1:n.vars,1:n.vars]
    rownames(m.add) <- var.names
    colnames(m.add) <- shk.names
    # The response
    df.add <- m.add %>%
      as.data.frame() %>%
      rownames_to_column('outcome') %>%
      gather( shock, value, -outcome) %>%
      mutate( period=i ) %>%
      select( shock, outcome, period, value )
    # In data frame form
    df.out <- df.out %>%
      rbind(df.add)
      # full_join( df.add )
  }
  df.out <- df.out %>% arrange(shock, outcome, period)
  this.fcast <- fcast
    # Need to do this to stop dplyr getting confused about names
  make.fire <- function(df, i){
  # Make the FIRE part
    fire.name <- paste0( inf.name[i], '.fire' )
    senti.name <- paste0( inf.name[i], '.sentiment')
    df.out <- df %>%
      group_by(shock) %>%
      mutate( !!(fire.name):= if(fcast.cumul[i]) c(rollsum(get(inf.name[i]), 
          fcast.horiz, fill = NA, align='l' )[-1], NA) else c( get(inf.name[i])[-(1:fcast.horiz)], 
                   rep(NA,fcast.horiz) ) ) %>%
      mutate( !!(senti.name) := get(this.fcast[i]) - get(fire.name) )
    return(df.out)
  }
  
  df.re.fcast.irf <- df.out %>% 
    filter(outcome %in% c( inf.name, fcast ) ) %>% 
    spread( outcome, value)
  for( i in 1:n.fcast){
    df.re.fcast.irf <- make.fire( df.re.fcast.irf, i )
  } 
    # group_by(shock) %>%
    # mutate( inf.fire=c(rollsum(get(inf.name), fcast.horiz, fill = NA, align='l' )[-1], NA),
    #         inf.sentiment = get(this.fcast) - inf.fire )
    # Create the rational expectations part of the IRFs.
  df.out.full <- full_join( df.out, 
                            df.re.fcast.irf %>% select(-c(fcast, inf.name)) %>% 
                              gather(outcome, value, -shock, -period ),
                            by=c('shock', 'outcome', 'period', 'value'))
    # Merge the RE decomp back in
  return( df.out.full )
}

make.bootstrap <- function( l.var, A, fcast.horiz, inf.name, n.pds=20, 
                            n.boot=100, seed=42, do.var.decomp=TRUE, print.iter=TRUE,
                            make.phi.fn=make.phi.h, do.nD=FALSE, fcast.cumul=TRUE,
                            reorder.to.idx=0 ){
  # Makes the bootstrapped IRFs and variance decomp
  
  ## 1. Set up
  set.seed(seed)
  # Set the seed
  est.resid <- l.var$var %>% resid
  # The primary residuals from the main estimation
  n.vars <- l.var$n.vars
  # The number of variables
  est.m.Y <- l.var$data
  # The real data
  n.pds.data <- nrow(est.m.Y)
  # The number of periods of the data
  p <- l.var$lags
  n.pds.boot <- n.pds.data - p
  # The number of periods to create the bootstrapped sample
  est.m.Y.x <- est.m.Y
  n.fcast <- inf.name %>% length()
  
  if(p>1) for( j in 2:p ) est.m.Y.x <- cbind( est.m.Y.x, rbind( matrix(NA, j-1, n.vars ), 
                                                                est.m.Y[1:(n.pds.data-j+1),] ))
  # The extended form of Y
  B <- l.var$B
  # The estimated transition matrix
  l.Sigma <- l.est.phi.h <- l.est.phi.h.theta <- l.var.decomp <- l.A.est <- l.B.est <- list()
  # The containers for the estimated persistence matrix B and the structural decomp A
  resid.x <- rep(0, n.vars*(p-1) )
  # The matrix to extend the residuals
  anc.0 <- A[-1,1]
  # The initial guess for the decomposition
  
  ## 2. Main loop to estimate the coefficients
  for( iter in 1:n.boot ){
    
    if(print.iter) message('Bootstrap, simulation=', iter)
    # if( iter==2 ) browser()
    
    this.resid <- est.resid %>% 
      as.data.frame %>% 
      sample_n( n.pds.boot, replace=TRUE) %>% 
      as.matrix()
    # the new sample of residuals
    this.m.Y.x <- 0 * est.m.Y.x
    this.m.Y.x[1:p, ] <- est.m.Y.x[1:p,]
    # Initialize the new data sample
    for( i in 1:n.pds.boot ) this.m.Y.x[p+i, ] <- B %*% this.m.Y.x[p+i-1,] + c( this.resid[i,], resid.x )
    # Fill in the data with the bootstrapped residuals
    this.l.var <- make.var( NA, NA, NA, NA, NA, NA, lags=p, m.Y=this.m.Y.x[,1:n.vars] )
    # Estimate the VAR on the new data
    l.B.est[[iter]] <- this.l.var$B
    # Record the estimate of the coefficient matrix
    l.Sigma[[iter]] <- this.l.var$Sigma
      # Record the variance matrix
    # if(do.nD){
      # browser()
      l.est.phi.h[[iter]] <- sapply( 1:n.fcast,
                           function(j) make.phi.fn( this.l.var, horiz = fcast.horiz, 
                                 inf.idx = n.fcast + j, cumul=fcast.cumul[j] )$phi.h ) %>%
        set_colnames( rownames(A)[1:n.fcast] ) %>% t
      l.A.est[[iter]] <- anc.analytic( l.Sigma[[iter]], l.est.phi.h[[iter]] ) %>%
        A.reorder(., l.Sigma[[iter]], n.fcast, reorder.to.idx )
    # }else{
    #   this.phi.h.out <- make.phi.fn( this.l.var, horiz=fcast.horiz )
    #   l.est.phi.h[[iter]] <- this.phi.h.out$phi.h
    #   l.est.phi.h.theta[[iter]] <- this.phi.h.out$theta
    #   # l.A.est[[iter]] <- anc.nlslv.robust( l.Sigma[[iter]], l.est.phi.h[[iter]], anc=anc.0 ) #, maxit=1000, gain=.05 )
    #   l.A.est[[iter]] <- anc.analytic( l.Sigma[[iter]], l.est.phi.h[[iter]] )
    #     # The structural decomp.  Want to be numerically conservative because some
    #     # guesses might not be PD.
    #   
    # }
    l.var.decomp[[iter]] <- var.decomp( l.A.est[[iter]]$A, this.l.var, n.pds=n.pds, n.fcast=n.fcast )
    # The variance decomposition
  }
  
  # 2.1 Check the structural decomp
    #### NEED TO RECOMPUTE FOR THE nD CASE ####
  
  A.status <- sapply(l.A.est,function(x) x$status )
  A.Sigma.err <- sapply(l.A.est,function(x) x$sigma.err )
  A.iter <- A.err <- n.fail <- 0
  # if(do.nD){
  # }else{
  #   A.err <- sapply(l.A.est,function(x) x$err )
  #   A.iter <- sapply(l.A.est, function(x) x$iter )
  #   n.fail <- ( A.status == 'failure' ) %>% sum
  #   message( 'Dropping ', n.fail, ' bootstrap trials due to failure of structural decomposition' )
  #   l.A.est <- l.A.est[ A.status == 'success' ]
  #   l.B.est <- l.B.est[ A.status == 'success' ]
  #   l.var.decomp <- l.var.decomp[ A.status == 'success' ]
  # }
  n.boot <- n.boot - n.fail
  
  ## 3. Compute the IRFs 
  
  ## 3.1 Reduced form
  message('Creating reduced form IRFs')
  df.boot.irfs.rf <- lapply( 1:n.boot, function(x){
    if(print.iter) message('Bootstrap, reduced form IRF, n=', x)
    make.irf.rf( list(n.vars=n.vars, B=l.B.est[[x]] ), n.pds ) %>% mutate(i.boot=x) 
  } ) %>%
    reduce(rbind)
    # The dataframe of all the RF IRFs
  df.irf.quantiles.rf <- df.boot.irfs.rf %>% group_by( shock, outcome, period ) %>%
    summarise( #q.005=quantile(value, .005), q.025=quantile(value, .025), 
               q.05=quantile(value, .05), q.50=quantile(value, .50), 
               q.95=quantile(value, .95), 
               # q.975=quantile(value, .975), q.995=quantile(value, .995),
               .groups = 'drop' )
    # The quantiles of the RF IRF
  # browser()
  
  ## 3.2 Structural 
  message('Creating structural IRFs')
  eval(fcast)
    # Force instantiation of this variable
  df.boot.irfs.struct <- lapply( 1:n.boot, function(i){
    if(print.iter) message('Bootstrap, structural IRF, n=', i) 
    make.irf.struct( list(n.vars=n.vars, B=l.B.est[[i]], lags=p ), 
                     l.A.est[[i]]$A, fcast.horiz, inf.name, 
                     n.pds, fcast, fcast.cumul ) %>%
      mutate(i.boot=i)
    } ) %>% 
    reduce(rbind)
  # The dataframe of all the IRFs
  df.irf.quantiles.struct <- df.boot.irfs.struct %>% group_by( shock, outcome, period ) %>%
    summarise( #q.005=quantile(value, .005, na.rm=T), q.025=quantile(value, .025, na.rm=T), 
               q.05=quantile(value, .05, na.rm=T), q.50=quantile(value, .50, na.rm=T), 
               q.95=quantile(value, .95, na.rm=T), 
               # q.975=quantile(value, .975, na.rm=T), q.995=quantile(value, .995, na.rm=T),
               .groups = 'drop' )
  # The quantiles of the structural IRF
  
  ## 4. Do rhw variance decomposition
  message('Variance decomposition')
  if(do.var.decomp){
    df.var.decomp <- 
      lapply( 1:n.boot, function(x){
        if(print.iter) message('Bootstrap, variance decomposition, n=', x)
        lapply( 1:length( l.var.decomp[[x]]$decomp.ratio), 
                function(i) l.var.decomp[[x]]$decomp.ratio[[i]] %>% 
                  as.data.frame %>%
                  rownames_to_column('outcome') %>%
                  gather( shock, value, -outcome ) %>%
                  mutate( shock = ifelse(str_detect(shock, 'Sentiment'), shock, 'Fundamental')) %>%
                  group_by( outcome, shock ) %>%
                  summarise(value=sum(value), .groups='drop') %>%
                  mutate( horizon=i ) ) %>%
          reduce(rbind) %>% 
          mutate(i.boot=x) %>%
          as_tibble() %>%
          select( outcome, shock, i.boot, horizon, value ) %>%
          arrange( outcome, shock, i.boot, horizon ) 
      } ) %>% 
      reduce(rbind)
    # The dataframe of the bootstrapped varianc decomp
    df.var.decomp.quantiles <- df.var.decomp %>% 
      group_by( shock, outcome, horizon ) %>%
      summarise( #q.005=quantile(value, .005), q.025=quantile(value, .025), 
                 q.05=quantile(value, .05), q.50=quantile(value, .50), 
                 q.95=quantile(value, .95), 
                 # q.975=quantile(value, .975), q.995=quantile(value, .995),
                 .groups = 'drop' )
    # The quantiles of the structural IRF
  }else{
    df.var.decomp.quantiles <- NA
  }
  
  return( list( reduced.form=df.irf.quantiles.rf, structural=df.irf.quantiles.struct, 
                var.decomp=df.var.decomp.quantiles, theta=unlist(l.est.phi.h.theta),
                report=list( status=A.status, err=A.err, sigma.err=A.Sigma.err, 
                             iter=A.iter, n.fail=n.fail ) ) )
  
}
