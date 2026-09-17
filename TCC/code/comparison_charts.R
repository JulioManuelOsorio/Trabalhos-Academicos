##########################################################################################
#
# comparison_charts.R
#
# Code to make charts comparing the IRFs in the  "Inflation Sentiments" project
# Philip Barrett, Washington DC
# First version: 13feb2022
#
##########################################################################################

#### 0. Set up ####
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(zoo)
library(xtable)
select <- dplyr::select
  # Otherwise picks up the MASS usage 

#### 1. Controls ####
v.st.switch <- c('rnr_2', 'lags', 'oil', 'fcasts' )
root.dir <- 'replication'

for(st.switch in v.st.switch){
  
  if(st.switch=='fcasts'){
    v.irfs <- paste0( 'model_solutions/', c( 'michigan_cpi_1982_m_3', 'cleveland_cpi_1982_m_4', 
                                  'spf_pgdp_1982_q_2', 'fed_pgdp_1982_q_4'  ), '.rdata')
      # The locations of the data
    l.compare <- list( 
      c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ),
      c( fcast='cleveland.fcast', inf='inf.cpi',  inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ),
      c( fcast='spf.fcast', inf='inf.gdp', inf.fire='inf.gdp.fire', inf.sentiment='inf.gdp.sentiment',
         act='log.gdp', fed.funds.rate='fed.funds.rate' ),
      c( fcast='fed.fcast', inf='inf.gdp', inf.fire='inf.gdp.fire', inf.sentiment='inf.gdp.sentiment',
         act='log.gdp', fed.funds.rate='fed.funds.rate' ) )
      # The details for comparing the different types of estimation
    v.freq <- c( ' Michigan'='m', Cleveland='m', SPF='q', 'Fed Greenbook'='q' )
      # The frequency
    do.lp <- FALSE
  }
  
  if(st.switch=='oil'){
    v.irfs <- paste0( 'model_solutions/', c( 'michigan_cpi_1982_m_3', 'michigan_cpi_oil_inf_2_1982_m_3' ), '.rdata')
    # The locations of the data
    l.compare <- list( 
      c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ),
      c( fcast='michigan.fcast', inf='inf.cpi',  inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ) )
    # The details for comparing the different types of estimation
    v.freq <- c( 'Baseline'='m', 'Including oil price inflation'='m' )
    # The frequency
    do.lp <- FALSE
  }
  
  
  if(st.switch=='lags'){
    v.irfs <- paste0( 'model_solutions/', c( 'michigan_cpi_1982_m_3', 'michigan_cpi_1lag_1982_m_1',
                                  'michigan_cpi_2lag_1982_m_2', #'michigan_cpi_4lag_1982_m_4',
                                  'michigan_cpi_6lag_1982_m_6','michigan_cpi_8lag_1982_m_8',
                                  'michigan_cpi_12lag_1982_m_12' ), '.rdata')
    # The locations of the data
    l.compare <- list( 
      c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ) )
    l.compare[[2]] <- l.compare[[3]] <- l.compare[[4]] <- l.compare[[5]] <- l.compare[[6]] <- #l.compare[[7]] <- 
      l.compare[[1]]
    # The details for comparing the different types of estimation
    v.freq <- c( '  Baseline (3 lags)'='m', '1 lag'='m', '2 lags'='m', # '4 lags'='m', 
                 '6 lags'='m', '8 lags'='m', '12 lags'='m')
    # The frequency
    do.lp <- TRUE
  } 
  
  if(st.switch=='rnr_2'){
    v.irfs <- paste0( 'model_solutions/', c( 'michigan_cpi_1982_m_3', 'michigan_cpi_1st_1982_m_3', 
                                  'michigan_cpi_trend_1982_m_3'), '.rdata')
    # , 'michigan_cpi_2010_1982_m_3', 'michigan_cpi_2015_1982_m_3'
    # The locations of the data
    l.compare <- list( 
      c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ),
      c( fcast='michigan.fcast', inf='inf.cpi.first', inf.fire='inf.cpi.first.fire', inf.sentiment='inf.cpi.first.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' ),
      c( fcast='michigan.fcast.trend', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
         act='log.ip', fed.funds.rate='fed.funds.rate' )
      # c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
      #    act='log.ip', fed.funds.rate='fed.funds.rate' ),
      # c( fcast='michigan.fcast', inf='inf.cpi', inf.fire='inf.cpi.fire', inf.sentiment='inf.cpi.sentiment',
      #    act='log.ip', fed.funds.rate='fed.funds.rate' )
      )
    # The details for comparing the different types of estimation
    v.freq <- c( ' Baseline'='m', 'Using initial CPI release'='m', 'Expectations detrended using inflation trend'='m',
                 'Until 2010'='m', 'Until 2015'='m')
    # The frequency
    do.lp <- FALSE
  }
  
  var.labs <- c( fcast='Inflation, one year ahead inflation cummulative,\n forecast', 
                 inf.a=' Inflation, annual rate',
                 act='Log Industrial Production * 100', 
                 fed.funds.rate='Federal Funds Rate', 
                 inf.fire='Inflation, one year ahead cummulative', 
                 inf.sentiment='Inflation, one year ahead cummulative,\n sentiment' )
  # The variable labels
  var.labs.tab <- c( fcast='Year-ahead inf. exp.', inf='Realized inflation',
                     act='100 * Industrial Production', fed.funds.rate='Federal Funds Rate' )
  # The variable labels for the table
  n.compare <- l.compare %>% length()
  save.file <- paste0( root.dir, '/compare/', st.switch )
  
  #### 2. Assemble the data ####
  
  l.sent.ar <- list()
  df.var.decomp.all <- NULL
  
  for( i.compare in 1:n.compare ){
    
    load( v.irfs[i.compare] )
      # Load the file
    freq <- v.freq[i.compare]
    inf.ann <- if(freq=='m') function(x) x*12 else function(x) x*4
    # inf.ann <- if(freq=='m') function(x) rollsum( c( rep(0,11), x ), 12, align='l' ) else function(x) rollsum( c( rep(0,3), x ), 4, align='l' )
    #   # Alternative specification: show annual inflation rate response
    pd.mth <- if(freq=='m') identity else function(x) x*3+1
      # Functions to adjust inflation to be annualized and periods to be monthly only
    this.df.nf.chart <- df.irf.struct %>% 
      # filter(shock=='Non-fundamental') %>% 
      filter(shock=='Sentiment #1') %>% 
      select(-outcome.lab) %>% 
      spread(outcome,value) %>%
      rename(l.compare[[i.compare]]) %>% 
      mutate( period.m=pd.mth(period), inf.a=inf.ann(inf) ) %>%
      select( period.m, fcast, inf.a, act, fed.funds.rate, inf.fire, inf.sentiment ) %>%
      gather( outcome, value, -period.m ) %>%
      mutate( type=names(v.freq)[i.compare],
              outcome=fct_relevel( outcome, c('fcast', 'inf.fire', 'inf.sentiment', 'inf.a', 'fed.funds.rate', 'act') ) )
    df.nf.chart <- if(i.compare==1) this.df.nf.chart else full_join( df.nf.chart, this.df.nf.chart )
    
    if(do.lp){
      this.est.lp <- est.lp %>%
        spread(outcome,value) %>%
        mutate( inf.a=inf.ann(inf.cpi), fcast=michigan.fcast, act=log.ip ) %>%
        select( period, fcast, inf.a, act, fed.funds.rate, inf.fire, inf.sentiment ) %>%
        gather( outcome, value, -period ) %>%
        mutate( type=names(v.freq)[i.compare],
                outcome=fct_relevel( outcome, c('fcast', 'inf.fire', 'inf.sentiment', 'inf.a', 'fed.funds.rate', 'act') ) )
      df.lp.chart <- if(i.compare==1) this.est.lp else full_join( this.est.lp, df.lp.chart )
    }
    
    l.sent.ar[[names(v.freq)[i.compare]]] <- this.df.nf.chart %>% filter( outcome == 'inf.sentiment', !is.na(value), value >=.01 ) %>% 
      pull(value) %>% ar(order.max=1, aic=FALSE)
    
    this.var.decomp <- 
      # full_join( 
        est.bootstrap$var.decomp %>% 
                 filter( horizon==max(horizon), shock!='Fundamental') %>% 
                 select(outcome, horizon, q.05, q.95, pt.est=q.50) %>% #,
                 # df.var.decomp %>% filter( horizon==max(horizon), shock!='Fundamental') %>%
                 # select( outcome, horizon, pt.est=value )
                 # ) 
      filter( outcome %in% l.compare[[i.compare]] ) %>%
      mutate( type=names(v.freq)[i.compare],
              outcome.tab=names(l.compare[[i.compare]])[match( outcome, l.compare[[i.compare]])],
              outcome.name=var.labs.tab[outcome.tab] ) %>%
      select(type, outcome.name, everything() )
    df.var.decomp.all <- if( i.compare > 1 ) full_join(df.var.decomp.all, this.var.decomp) else this.var.decomp
    
  }
  
  
  gg.compare <- 
  ggplot(df.nf.chart, aes( x=period.m, y=value, group=type, color=type, shape=type)) +
    geom_line(lwd=.9) +
    geom_point(size=3) +
    facet_wrap(~outcome, label=labeller( outcome=var.labs ), scales = 'free_y' ) +
    theme_minimal() + 
    theme(legend.position = 'bottom') +
    labs( color='Measure', lty='Measure', shape='Measure' ) +
    xlab('Months') + ylab('') +
    geom_hline(yintercept=0) +
    # xlim(c(0,24)) +
    scale_x_continuous( breaks=seq(0,24,6), limits = c(0,24) )
  gg.compare %>% print()
  ggsave(gg.compare, filename = paste0( 'graphs/', save.file, '_irf_compare.pdf'), height=6, width=10 )
  
  if(do.lp){
    gg.compare.lp <- 
      ggplot(df.lp.chart, aes( x=period, y=value, group=type, color=type, shape=type)) +
      geom_line(lwd=.9) +
      geom_point(size=3) +
      facet_wrap(~outcome, label=labeller( outcome=var.labs ), scales = 'free_y' ) +
      theme_minimal() + 
      theme(legend.position = 'bottom') +
      labs( color='Measure', lty='Measure', shape='Measure' ) +
      xlab('Months') + ylab('') +
      geom_hline(yintercept=0) +
      # xlim(c(0,24)) +
      scale_x_continuous( breaks=seq(0,24,6), limits = c(0,24) )
    gg.compare.lp %>% print()
    ggsave(gg.compare.lp, filename = paste0( 'graphs/', save.file, '_lp_compare.pdf'), height=6, width=10 )
  }
  
  
  
  if(st.switch=='fcasts'){
    
    load('data/baseline_time_series_update.rdta')  
    
    df.all.fcast <- df.m %>% filter( variable %in% c('michigan.fcast', 'cleveland.fcast')) %>%
      full_join( ., df.q %>% filter( variable %in% c('spf.fcast', 'fed.fcast') ) %>% mutate(value=as.numeric(value)) ) %>%
      mutate( plt.var=case_when( variable=='michigan.fcast'~' Michigan',
                                 variable=='cleveland.fcast'~'Cleveland Fed',
                                 variable=='spf.fcast'~'SPF',
                                 variable=='fed.fcast'~'Fed Greenbook' ) ) %>%
      ungroup() %>%
      select(date, plt.var, value, value.deseas ) %>%
      gather( type, value, -date, -plt.var )
    gg.fcast.compare <- ggplot( df.all.fcast, 
                                aes(x=date, y=value, color=plt.var, type=plt.var, shape=plt.var) ) +
      geom_line(lwd=.9) +
      # geom_point(size=3) +
      facet_wrap(~type, label=labeller( type=c( value='Raw data', value.deseas='Deseasonalized and detrended') )) +
      theme_minimal() +
      theme(legend.position = 'bottom') +
      geom_hline(yintercept=0, lwd=.5) +
      labs( color='Measure') +
      xlab('Date')
    gg.fcast.compare %>% print()
    ggsave(gg.fcast.compare, filename = paste0('graphs/', save.file, '_fcast_compare.pdf'), height=6, width=10 )
    
    v.sent.ar <- l.sent.ar %>% sapply( ., function(x) x$ar )
      # The sentiment persistences
    df.var.decomp.tab <- df.var.decomp.all %>%
      ungroup() %>%
      mutate( range = paste0( '(', formatC(q.05, 2, format = 'f'), ', ', 
                              formatC(q.95, 2, format='f'), ')' ),
              pt.est=formatC(pt.est,2,format='f')) %>% 
      select(type, outcome.name, pt.est, range ) %>%
      gather( key, value, -type, -outcome.name ) %>%
      spread( type, value ) %>%
      select(-key) %>%
      mutate( outcome.name=ifelse( outcome.name==c(outcome.name[-1],' '), outcome.name, '')) %>%
      select( everything(), -'Fed Greenbook', 'Fed Greenbook') %>%
      rename( ' '=outcome.name)
    
    table.string <- 'Long-run share of variance, sentiment shock. 
              Bootstrapped 90 percent confidence interval in parentheses.'
    print( xtable( df.var.decomp.tab, label='t:var_decomp',
                   caption = table.string, align = c('l','l','c','c','c','c') ),
          include.rownames=FALSE, hline.after=c(-1,-1,0,8), floating = FALSE, 
          file = paste0('tables/', st.switch, '_var_decomp.tex' ) ) #,
          # floating.environment = 'sidewaystable' )
  }
  
}  
  
  
  
  
  
  
