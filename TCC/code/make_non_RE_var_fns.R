##########################################################################################
#
# make_non_RE_var_fns.R
#
# Code to define the functions to make non-RE VAR decompositions
# Philip Barrett, Washington DC
# First version: 29aug2022
#
##########################################################################################

library(tidyverse)
library(expm)
select <- dplyr::select
  # Otherwise picks up the MASS usage

make.phi.h.DO <- function( l.var, horiz=12, inf.idx=2, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for delayed observations
  
  n.vars <- l.var$n.vars
    # Number of variables
  out <- list()
  out$phi.h <- rep(0, n.vars) ; names(out$phi.h) <- l.var$B %>% rownames() %>% head(n.vars)
  return( out )
  # Select the inflation response
}

make.phi.h.PO <- function( l.var, this.df=df, z.var=c('fed.funds.rate'), horiz=12, inf.idx=2, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for partially
# observed data.
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  rf.resids <- l.var$var %>% resid()
    # The reduced form residuals
  n.z <- length(z.var)
    # Number of instruments
  z <- matrix( l.var$data[, z.var], nrow=nrow(l.var$data), ncol=n.z ) %>%
    set_colnames(z.var)
    # this.df %>% select(-value) %>% spread(variable, value.deseas) %>% select(!!z.var) %>%
    # as.matrix()
    # The z variables
  m.rhs <- lapply( 1:l.var$lags, function(i) rbind( matrix( NA, nrow=i, ncol=n.vars ), 
                                                    l.var$data[-(nrow(l.var$data)+1-(1:i)),] ) ) %>%
    reduce(cbind)
    # The VAR information set
  z.tilde <- sapply( 1:n.z, function(i) lm( z[,i] ~ m.rhs ) %>% resid )
    # The residual information in z
  B.uz <- sapply( 1:n.vars, function(i) lm( rf.resids[,i] ~ 0 + z.tilde ) %>% coef )
  B.zu <- sapply( 1:n.z, function(i) lm( z.tilde[,i] ~ 0 + rf.resids ) %>% coef )
    # The regressions of instruments on residuals and vice versa
  if( n.z==1 ) B.uz <- t(B.uz)
  phi.h.PO <- ( t(B.uz) %*% t(B.zu) %*% phi.h.re ) %>% c
  names(phi.h.PO) <- l.var$B %>% rownames() %>% head(n.vars)
    # The new phi
  out <- list( phi.h=phi.h.PO )
  return( out )
}

make.phi.h.PO.fed.funds <- function( l.var, this.df=df, horiz=12, inf.idx=2, cumul=NA ){
  return(make.phi.h.PO( l.var, this.df=this.df, z.var=c('fed.funds.rate'), horiz=12, inf.idx=2 ) )
}

make.phi.h.PO.cdty <- function( l.var, this.df=df, horiz=12, inf.idx=2, cumul=NA ){
  return(make.phi.h.PO( l.var, this.df=this.df, z.var=c('fed.funds.rate', 'inf.commod'), horiz=12, inf.idx=2 ) )
}

make.phi.h.PO.big <- function( l.var, this.df=df, horiz=12, inf.idx=2, cumul=NA ){
  return(make.phi.h.PO( l.var, this.df=this.df, 
          z.var=c('fed.funds.rate', 'i.3.mo', 'i.1.yr', 'i.2.yr', 'i.5.yr', 'i.10.yr', 'i.30.yr',
                  'log.will.idx', 'log.fx.jpn', 'log.fx.gbr', 'log.fx.cad' ), 
          horiz=12, inf.idx=2 ) )
}

make.phi.h.PO.small <- function( l.var, this.df=df, horiz=12, inf.idx=2, cumul=NA ){
  return(make.phi.h.PO( l.var, this.df=this.df, z.var=c('fed.funds.rate', 'i.1.yr', 'log.will.idx' ),  
                        horiz=12, inf.idx=2 ) )
}

make.phi.h.CD <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for cognitive
# discounting. Assuming that forecasts are ordered first
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  rf.resids <- l.var$var %>% resid()
    # The reduced form residuals
  fcast.resid <- rf.resids[,fcast.idx]
    # The residuals for the forecast
  reg.rhs <- rf.resids %*% phi.h.re
    # The rational responses to the residuals
  lm.CD <- lm( fcast.resid ~ 0 + reg.rhs )
    # The auxiliary regression
  theta.CD <- lm.CD$coef[1]
    # The attenuation from the cognitive discounting
  phi.h.CD <- phi.h.re * theta.CD
  names(phi.h.CD) <- l.var$B %>% rownames() %>% head(n.vars)
    # The new phi
  out <- list( phi.h=phi.h.CD, theta=theta.CD, aux.reg=lm.CD )
  return( out )
    # The output
}

make.phi.h.CD.reg <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for cognitive
# discounting via explicit regression on rational forecasts
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  var.fcast <- make.var.fcast.fit( l.var, horiz, TRUE )
    # The forecast from the VAR
  pi.e <- l.var$data[,fcast.idx]
    # The residuals for the forecast
  lm.CD <- lm( pi.e ~ 0 + var.fcast[,inf.idx] )
    # The auxiliary regression
  theta.CD <- lm.CD$coef[1]
    # The attenuation from the cognitive discounting
  phi.h.CD <- phi.h.re * theta.CD
  names(phi.h.CD) <- l.var$B %>% rownames() %>% head(n.vars)
    # The new phi
  out <- list( phi.h=phi.h.CD, theta=theta.CD, aux.reg=lm.CD )
  return( out )
  # The output
}

make.phi.h.AE <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for adaptive
# expectations.
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  rf.resids <- l.var$var %>% resid()
    # The reduced form residuals
  fcast.resid <- rf.resids[,fcast.idx]
    # The residuals for the forecast
  reg.rhs <- rf.resids[,inf.idx]
    # The inflation residuals
  lm.AE <- lm( fcast.resid ~ 0 + reg.rhs )
    # The auxiliary regression
  theta.AE <- 1 - lm.AE$coef[1]
    # The attenuation from adaptive expectations
  phi.h.AE <- 0 * phi.h.re
  phi.h.AE[inf.idx] <- 1-theta.AE
  names(phi.h.AE) <- l.var$B %>% rownames() %>% head(n.vars)
    # The new phi
  out <- list( phi.h=phi.h.AE, theta=theta.AE, aux.reg=lm.AE )
  return( out )
    # The output
}

make.phi.h.AE.reg <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, n.lags.ae=6, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for adaptive
# expectations via regression
  
  n.vars <- l.var$n.vars
      # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
      # The rational phi
  pi.e <- l.var$data[,fcast.idx]
      # The forecast
  pi <- l.var$data[,inf.idx]
      # The inflation series
  df.nls <- sapply( 0:n.lags.ae, function(i) if(i==0) pi else c( rep(NA, i), pi[-(1:i)] ) ) %>%
    set_colnames( paste0( 'pi.l', 0:n.lags.ae) ) %>%
    as_tibble() %>%
    mutate( pi.e=pi.e )
      # The dataset of lags and forecasts
  nls.fla <- paste0( 'pi.e ~ (1-theta) * ( ', paste0( 'theta^', 0:n.lags.ae, 
                         '*pi.l', 0:n.lags.ae, collapse = ' + ' ) , ' )' )
  nls.AE <- nls( nls.fla %>% as.formula(), data = df.nls, start=c(theta=.6), 
                 upper=c(theta=1), lower=c(theta=0) )
      # The nonlinear regression
  theta.AE <- nls.AE %>% coef
      # The attenuation from adaptive expectations
  phi.h.AE <- 0 * phi.h.re
  phi.h.AE[inf.idx] <- 1-theta.AE
  names(phi.h.AE) <- l.var$B %>% rownames() %>% head(n.vars)
      # The new phi
  out <- list( phi.h=phi.h.AE, theta=theta.AE, aux.reg=nls.AE )
  return( out )
    # The output
}

make.phi.h.AE.recursive <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for adaptive
# expectations via regression on the recursive formulation
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  pi.e <- l.var$data[,fcast.idx]
    # The forecast
  pi <- l.var$data[,inf.idx]
    # The inflation series
  df.nls <- data.frame( pi.e=pi.e, pi=pi, pi.e.l1=c(NA, pi.e[-length(pi.e)]) )
    # The dataset of lags and forecasts
  nls.AE <- nls( pi.e ~ (1-theta) * pi + theta * pi.e.l1, data = df.nls, start=list(theta=.6) )
    # The nonlinear regression
  theta.AE <- nls.AE %>% coef
    # The attenuation from adaptive expectations
  phi.h.AE <- 0 * phi.h.re
  phi.h.AE[inf.idx] <- 1-theta.AE
  names(phi.h.AE) <- l.var$B %>% rownames() %>% head(n.vars)
  # The new phi
    out <- list( phi.h=phi.h.AE, theta=theta.AE, aux.reg=nls.AE )
  return( out )
    # The output
}

make.phi.h.AE.z <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1,
                                       z=c( 'fed.funds.rate'), cumul=NA ){
  # Makes the cumulative response for inflation at horizon h for adaptive
  # expectations via regression on the recursive formulation
  
  n.vars <- l.var$n.vars
    # Number of variables
  n.z <- z %>% length
    # Number of other variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  pi.e <- l.var$data[,fcast.idx]
    # The forecast
  df.nls <- data.frame( pi.e=pi.e, l.var$data[,c(inf, z)], pi.e.l1=c(NA, pi.e[-length(pi.e)]) )
    # The dataset of lags and forecasts
  fla.nls <- paste0( 'pi.e ~ (1-theta) * ( ', inf, ' + ', paste0( 'beta.', 1:n.z, ' * ', z, collapse=' + ' ), ' ) + theta * pi.e.l1' )
    # The restricted regression formula
  l.start <- lapply( 1:n.z, function(i) 0 )
  names(l.start) <- paste0( 'beta.', 1:n.z)
  l.start$theta <- .9
    # Initial values for the parameters
  nls.AE <- nls( fla.nls %>% as.formula, data = df.nls, start=l.start )
    # The nonlinear regression
  theta.AE <- (nls.AE %>% coef)['theta']
  beta.AE <- (nls.AE %>% coef)[paste0('beta.',1:n.z)]
    # The attenuation from adaptive expectations
  phi.h.AE <- 0 * phi.h.re
  phi.h.AE[inf.idx] <- (1-theta.AE)
  phi.h.AE[z] <- (1-theta.AE) * (beta.AE)
  names(phi.h.AE) <- l.var$B %>% rownames() %>% head(n.vars)
    # The new phi
  out <- list( phi.h=phi.h.AE, theta=theta.AE, beta=beta.AE, aux.reg=nls.AE )
  return( out )
  # The output
}

make.phi.h.AE.fix <- function( l.var, theta, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the adaptive expectations for inflation at horizon h for diagnostic
# discounting for fixed theta
  
  n.vars <- l.var$n.vars
      # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
      # The rational phi
  theta.AE <- theta
      # The attenuation from adaptive expectations
  phi.h.AE <- 0 * phi.h.re
  phi.h.AE[inf.idx] <- 1-theta.AE
  names(phi.h.AE) <- l.var$B %>% rownames() %>% head(n.vars)
      # The new phi
  out <- list( phi.h=phi.h.AE, theta=theta.AE )
  return(out)
}

make.phi.h.DE <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for diagnostic
# discounting. Assuming that forecasts are ordered first
  
  n.vars <- l.var$n.vars
  # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
  # The rational phi
  rf.resids <- l.var$var %>% resid()
  # The reduced form residuals
  fcast.resid <- rf.resids[,fcast.idx]
  # The residuals for the forecast
  reg.rhs <- rf.resids %*% phi.h.re
  # The rational responses to the residuals
  lm.DE <- lm( fcast.resid ~ 0 + reg.rhs )
  # The auxiliary regression
  theta.DE <- lm.DE$coef[1] - 1
  # The diagnostic coefficient
  B.h <- l.var$B %^% horiz
  phi.h.DE <- phi.h.re + theta.DE * B.h[ inf.idx, 1:n.vars ] 
  names(phi.h.DE) <- l.var$B %>% rownames() %>% head(n.vars)
  # The new phi
  out <- list( phi.h=phi.h.DE, theta=theta.DE, aux.reg=lm.DE )
  return( out )
  # The output
}

make.phi.h.DE.reg <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for diagnostic
# discounting. Assuming that forecasts are ordered first
  
  n.vars <- l.var$n.vars
    # Number of variables
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  var.fcast <- make.var.fcast.fit( l.var, horiz, TRUE )
  var.fcast.1 <- make.var.fcast.fit( l.var, 1, TRUE )
  var.fcast.13 <- make.var.fcast.fit( l.var, horiz+1, TRUE )
    # The forecast from the VAR
  
  pi.e <- l.var$data[,fcast.idx]
    # The residuals for the forecast
  df.nls <- data.frame( pi.e=pi.e, re.fcast.l0=var.fcast[,inf.idx], 
                    re.fcast.l1=c(NA,(var.fcast.13-var.fcast.1)[-nrow(var.fcast),inf.idx]) )
  # y <- pi.e - var.fcast[,inf.idx]
  # x <- c(NA, diff( var.fcast[,inf.idx] ) )
  #   # The difference form of the regression
  # # lm.DE <- lm( pi.e ~ 0 + var.fcast[,inf.idx] )
  # lm.DE <- lm( y ~ 0 + x )
    ### OLD, WRONG REGRESSIONS
  nls.DE <- nls( pi.e ~ (1+theta) * re.fcast.l0 - theta * re.fcast.l1, df.nls, start=list(theta=-.2) )
    # The auxiliary regression
  
  theta.DE <- nls.DE %>% coef
    # The diagnostic coefficient
  # B.h <- l.var$B %^% horiz
  # phi.h.DE <- phi.h.re + theta.DE * B.h[ inf.idx, 1:n.vars ] 
    ### OLD, WRONG FROMULATION FOR PHI
  phi.h.DE <- (1+theta.DE) * phi.h.re
    # The new phi
  out <- list( phi.h=phi.h.DE, theta=theta.DE, aux.reg=nls.DE )
  return( out )
    # The output
}

make.phi.h.PO.DE.fix <- function( l.var, theta, z.var=c('fed.funds.rate'), horiz=12, inf.idx=2, cumul=NA ){
  # Makes the cumulative response for inflation at horizon h for partially
  # observed data.
  
  n.vars <- l.var$n.vars
  # Number of variables
  theta.DE <- theta
  # The diagnostic coefficient
  phi.h.PO <- make.phi.h.PO( l.var, z.var = z.var,horiz = horiz, inf.idx=inf.idx )$phi.h
  # The rational phi
  phi.h.DE.PO <- (1+theta.DE) * phi.h.PO
  # The new phi
  out <- list( phi.h=phi.h.DE.PO, theta=theta.DE )
  return( out )
}

make.phi.h.DE.fix <- function( l.var, theta, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
# Makes the cumulative response for inflation at horizon h for diagnostic
# discounting for fixed theta
  
  n.vars <- l.var$n.vars
    # Number of variables
  theta.DE <- theta
    # The diagnostic coefficient
  phi.h.re <- make.phi.h( l.var, horiz = horiz, inf.idx=inf.idx )$phi.h
    # The rational phi
  phi.h.DE <- (1+theta.DE) * phi.h.re
    # The new phi
  out <- list( phi.h=phi.h.DE, theta=theta.DE )
  return( out )
    # The output
}

## NB: Fixed theta versions need global values to be set.
make.phi.h.DE.cg <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.DE.fix( l.var, theta.cg, horiz=horiz, inf.idx=inf.idx, fcast.idx=fcast.idx ) )
}
make.phi.h.DE.bord <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.DE.fix( l.var, theta.bord, horiz=horiz, inf.idx=inf.idx, fcast.idx=fcast.idx ) )
}
make.phi.h.DE.bord.indiv <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.DE.fix( l.var, theta.bord.indiv, horiz=horiz, inf.idx=inf.idx, fcast.idx=fcast.idx ) )
}
make.phi.h.AE.christiano <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.AE.fix( l.var, theta.christiano, horiz=horiz, inf.idx=inf.idx, fcast.idx=fcast.idx ) )
}
make.phi.h.AE.gelain <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.AE.fix( l.var, theta.gelain, horiz=horiz, inf.idx=inf.idx, fcast.idx=fcast.idx ) )
}

make.phi.h.PO.DE.bord <- function( l.var, horiz=12, inf.idx=2, fcast.idx=1, cumul=NA ){
  return( make.phi.h.PO.DE.fix( l.var, theta.bord, 'fed.funds.rate', horiz=horiz, inf.idx=inf.idx ) )
}