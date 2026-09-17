##########################################################################################
#
# lp_irfs.R
#
# Code to make local projection IRFs in the  "Inflation Sentiments" project
# Philip Barrett, Washington DC
# First version: 14apr2022
#
##########################################################################################

library(lpirfs)

make.lp.irf <- function( df.in, est.A.in, var.coefs.in, inf.name, this.fcast, 
                         irf.pds, p.conf=.95, struct.idx=1 ){
# Makes the IRFs for the structural shock and returns in the appropriate format
# for plotting with everything else
  
  ### Make the shocks
  struct.resid <- solve( est.A.in$A, var.coefs.in$var %>% resid() %>% t ) %>% t
  struct.shock <- c( rep(NA, var.coefs.in$lags), struct.resid[,struct.idx] ) %>%
    as.data.frame()
  
  ### Make the LPs
  confint <- qnorm(1-(1-p.conf)/2)
  results_lin_iv <- lp_lin_iv(endog_data = l.var.coefs$data %>% 
                                as.data.frame(), 
                              lags_endog_lin = var.coefs.in$lags,
                              shock = struct.shock, trend = 0,
                              confint = confint, hor = irf.pds+1 )
  
  ### Assemble into a dataframe that is useful
  df.out <- results_lin_iv$irf_lin_mean %>% t %>%
    set_colnames(var.coefs.in$data %>% colnames ) %>%
    data.frame( period=0:irf.pds, . ) %>%
    mutate( inf.fire=c(rollsum(get(inf.name), fcast.horiz, fill = NA, align='l' )[-1], NA),
            inf.sentiment = get(this.fcast) - inf.fire ) %>%
    gather(outcome, value, -period ) %>%
    as_tibble()
  # Create the rational expectations part of the IRFs.
  
  return(df.out)  
}
