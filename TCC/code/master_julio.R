##########################################################################################
#
# master.R  (Versão adaptada para ler Excel/CSV em vez de .rdta)
#
##########################################################################################

#### 0. Set up global options ####

## 0.1 Housekeeping ##
rm(list=ls())
library(tidyverse)
# instalar pacote apenas se necessário
if(!requireNamespace("lpirfs", quietly = TRUE)) install.packages('lpirfs')

# Ajuste o root.dir se quiser um caminho absoluto; por padrão uso o diretório corrente
root.dir <- 'C:\\Users\\julio\\OneDrive\\Desktop\\USP\\TCC'
setwd(root.dir)

# --- CONFIGURE AQUI: caminho para o seu arquivo Excel/CSV ---
# Você pode colocar: 'data/baseline_time_series_update.rdta' (original) 
# ou 'data/final_data.xlsx' (seu Excel). 
data.file <- 'data_tcc/new_final_data_trimmed.xlsx'     # <<-- troque se necessário

xl.in.dta <- 'data_tcc'
#### TO ADD: CLEAR THE OLD GRAPHS AND TABLES ####

## 0.2 Define functions ####
source('code/R/functions/make_var_fns.R')
source('code/R/functions/make_irf_fns.R')
source('code/R/functions/make_lp_fns.R')
source('code/R/functions/make_non_RE_var_fns.R')

## 0.3 Controls ####
make.data <- FALSE # FALSE # Remake the data from scratch?
out_prefix <- 'replication/'               # Where to store the outputs
n.boot <- 50 # 100 # 500 # Set a small number for testing

# ----- CASES: ATENÇÃO: ajustei os nomes de fcast/x/y para nomes plausíveis do seu Excel ----
# Se preferir, edite o bloco abaixo para usar as variáveis que você tem no Excel.
l.cases <- list( 
  basecase_consensus <- list(
    freq='m', fcast=c('consensus_expec_inflation'),
    x='inflation',
    fcast.cumul=c(TRUE),
    y=c('output_gap', 'selic_rate', 'delta_er', 'credit_growth', 'embi'),
    exog = c('commodities_index', 'fed_funds_rate', 'vix'),
    init.yr=NA, robust=TRUE,
    do.lp=TRUE, do.non.re=TRUE, non.re.ci=TRUE),
  basecase_survey <- list(
    freq='m', fcast=c('survey_expec_inflation'),
    x='inflation',
    fcast.cumul=c(TRUE),
    y=c('output_gap', 'selic_rate', 'delta_er', 'credit_growth', 'embi'),
    exog = c('commodities_index', 'fed_funds_rate', 'vix'),
    init.yr=NA, robust=TRUE,
    do.lp=TRUE, do.non.re=TRUE, non.re.ci=TRUE),
    selec1_consensus <- list(
      freq='m', fcast=c('consensus_expec_inflation'),
      x='inflation',
      fcast.cumul=c(TRUE),
      y=c('output_gap', 'selic_rate', 'delta_er', 'embi'),
      exog = c('commodities_index', 'fed_funds_rate', 'vix'),
      init.yr=NA, robust=TRUE,
      do.lp=TRUE, do.non.re=TRUE, non.re.ci=TRUE),
  selec1_survey <- list(
    freq='m', fcast=c('survey_expec_inflation'),
    x='inflation',
    fcast.cumul=c(TRUE),
    y=c('output_gap', 'selic_rate', 'delta_er', 'embi'),
    exog = c('commodities_index', 'fed_funds_rate', 'vix'),
    init.yr=NA, robust=TRUE,
    do.lp=TRUE, do.non.re=TRUE, non.re.ci=TRUE)
  
  
      
    )

## 0.4 Settings for the non-RE functions ##
theta.bord <- .55 ; theta.bord.indiv <- - 0.15
theta.christiano <- .83 ; theta.gelain <- .91
l.non.re.fns <- c('Delayed Observation'=make.phi.h.DO, 
                  'Adaptive expectations, Gelain et al. 2019'=make.phi.h.AE.gelain,
                  'Diagnostic Expectations, over-reaction'=make.phi.h.DE.bord,
                  'Diagnostic Expectations, under-reaction'=make.phi.h.DE.bord.indiv)
v.non.re.idx.show <- 1:length(l.non.re.fns)
l.non.re.show <- l.non.re.fns[v.non.re.idx.show] %>% names()

## 0.5 Derived settings ##
n.cases <- l.cases %>% length()           
make.dir.if.needed <- function(dir.path) if(!file.exists(dir.path)) dir.create(dir.path,recursive = TRUE)
for( this.vars in c( 'compare', names(l.cases)) ) make.dir.if.needed(paste0('graphs/',out_prefix, this.vars))

########## Main Computational loop
for(i.case in 1:n.cases ){
  
  #### 1. Set up ####
  freq <- l.cases[[i.case]]$freq
  fcast <- l.cases[[i.case]]$fcast
  fcast.cumul <- l.cases[[i.case]]$fcast.cumul
  init.yr <- l.cases[[i.case]]$init.yr
  final.yr <- if( is.null(l.cases[[i.case]]$final.yr) ) NA else l.cases[[i.case]]$final.yr
  x <- l.cases[[i.case]]$x
  y <- l.cases[[i.case]]$y
  do.lp <- l.cases[[i.case]]$do.lp
  non.re.ci <- l.cases[[i.case]]$non.re.ci
  do.non.re <- l.cases[[i.case]]$do.non.re
  robust <- if( is.null(l.cases[[i.case]]$robust) ) FALSE else l.cases[[i.case]]$robust
  irf.pds <- if(freq=='m') 36 else 16
  irf.plot.pds <- irf.pds - if(freq=='m') 12 else 4
  n.fcast <- fcast %>% length()
  n.vars <- 2*n.fcast + length(y)
  n.non.re <- l.non.re.fns %>% length()
  this.case.name <- names(l.cases)[i.case]
  order.first <- if( is.null(l.cases[[i.case]]$order.first) ) 'none' else l.cases[[i.case]]$order.first
  lags <- if(is.null(l.cases[[i.case]]$lags)) NA else l.cases[[i.case]]$lags
  
  message( '\n\n******** CASE # ', i.case, ' ********')
  message( 'Forecast series = ', paste0( fcast, collapse=', '),
           ', Initial year = ', init.yr, '\n\n' )
  
  
  #### 2. Read and process the data ####
  # -------------------------------------------------
  # Objetivo: suportar 3 formatos:
  #  - arquivo .RData/.rdata/.rdta (mantém compatibilidade)
  #  - arquivo .csv
  #  - arquivo .xlsx / .xls  (Excel)
  # O resultado deve criar df.m (monthly long) e df.q (quarterly long).
  # Formato esperado: colunas = date, variable, value.deseas
  # -------------------------------------------------
  ext <- tools::file_ext(data.file) %>% tolower()
  
  if(ext %in% c('rdata','rdta','rda')){
    # Mantém comportamento original
    load(data.file)   # espera que este load crie df.m e/ou df.q
  } else if(ext %in% c('csv','xlsx','xls')){
    # bibliotecas necessárias
    library(readxl)
    library(lubridate)
    library(zoo)
    
    if(ext %in% c('xlsx','xls')){
      df_wide <- readxl::read_excel(data.file)
    } else {
      df_wide <- read_csv(data.file, show_col_types = FALSE)
    }
    
    # Normalizar nomes de colunas: remover espaços e transformar para snake_case simples
    names(df_wide) <- names(df_wide) %>% 
      str_replace_all('\\s+','_') %>% str_replace_all('\\.','_') %>% tolower()
    
    # Tentar detectar coluna de data automaticamente (colunas com 'date', 'data', 'mes', 'period', 'ano')
    date.candidates <- names(df_wide)[grepl('date|data|mes|period|ano', names(df_wide))]
    if(length(date.candidates) >= 1){
      date.col <- date.candidates[1]
    } else {
      stop("Não detectei automaticamente a coluna de datas no seu Excel. Renomeie a coluna de datas para 'Data' ou 'date' e tente novamente.")
    }
    
    # Converter para Date (tenta formatos usuais)
    df_wide[[date.col]] <- tryCatch(as.Date(df_wide[[date.col]]),
                                    error = function(e){
                                      # tentar dd/mm/YYYY
                                      try(as.Date(df_wide[[date.col]], format = "%d/%m/%Y"), silent = TRUE)
                                    })
    if(all(is.na(df_wide[[date.col]]))){
      stop("Falha ao converter a coluna de data. Verifique o formato de datas no Excel (ex.: YYYY-MM-DD ou DD/MM/YYYY).")
    }
    # renomear para 'date'
    df_wide <- df_wide %>% rename(date = !!date.col)
    
    # Arrange por data
    df_wide <- df_wide %>% arrange(date)
    
    # Escolher colunas numéricas (todas exceto 'date') — cuidado: mantenha as que deseja usar
    value_cols <- names(df_wide)[ names(df_wide) != 'date' ]
    
    # Pivotar para formato longo com coluna 'value.deseas'
    df_long <- df_wide %>%
      pivot_longer(cols = all_of(value_cols), names_to = 'variable', values_to = 'value.deseas') %>%
      mutate(variable = as.character(variable))
    
    # Criar df.m (monthly long)
    # Garantir que a frequência seja mensal; se tiver duplicatas de data/variable, deixamos como está
    df.m <- df_long
    
    # Criar df.q (quarterly) agregando média por trimestre (último dia do trimestre como date)
    df.q <- df_long %>%
      mutate(yearq = as.yearqtr(date)) %>%
      group_by(yearq, variable) %>%
      summarise(value.deseas = mean(value.deseas, na.rm = TRUE), .groups='drop') %>%
      mutate(date = as.Date(as.yearqtr(yearq, format = "%Y Q%q"), frac = 1)) %>%
      select(date, variable, value.deseas)
    
    # Exibir resumo para checagem
    message("Leitura feita de: ", data.file)
    message("Observações (datas) carregadas: ", min(df_wide$date), " até ", max(df_wide$date))
    message("Variáveis detectadas: ", paste(unique(df_long$variable)[1:20], collapse=', '), if(length(unique(df_long$variable))>20) " ...")
    
  } else {
    stop("Formato de arquivo não suportado: ", ext, ". Use .RData/.rda/.rdta ou .csv ou .xlsx.")
  }
  
  # Definir janelas inicial/final com base nas séries carregadas
  df <- if(freq=='m') df.m else df.q
  init.yr <- if(is.na(init.yr)) df$date %>% min %>% year else init.yr
  final.yr <- if(is.na(final.yr)) df$date %>% max %>% year else final.yr
  df <- df %>% filter( year( date ) >= init.yr, year( date ) <= final.yr ) 
  
  #### 3. Estimate the reduced forms ####
  if(!is.na(lags)){
    l.var.coefs <- make.var( na.omit(df), value.name='value.deseas', 
                             var.name='variable', fcast=fcast, 
                             inf=x, y=y, lags=lags )
  }else{
    l.var.coefs <- make.var( na.omit(df), value.name='value.deseas', 
                             var.name='variable', fcast=fcast, 
                             inf=x, y=y, lag.max=8 )
  }
  
  #### 4. Compute the structural decomposition ####
  fcast.horiz <- if(freq=='m') 12 else 4
  est.phi.h <- sapply( 1:n.fcast,
                       function(i) make.phi.h( l.var.coefs, horiz = fcast.horiz, 
                                               inf.idx = n.fcast + i, cumul=fcast.cumul[i] )$phi.h ) %>%
    set_colnames( fcast ) %>% t
  reorder.to.idx <- if(order.first=='none') 0 else n.fcast * 2 + which(y==order.first)
  est.A <- anc.analytic( l.var.coefs$Sigma, est.phi.h ) %>%
    A.reorder(., l.var.coefs$Sigma, n.fcast, reorder.to.idx )
  
  est.var.decomp <- var.decomp( est.A$A, l.var.coefs, n.pds = irf.pds, n.fcast=n.fcast )
  if(do.non.re){
    l.non.re.phi <- lapply( l.non.re.fns, function(fn) fn(l.var.coefs, horiz = fcast.horiz) )
    l.est.A.non.re <- lapply( l.non.re.phi, function(x) anc.analytic( l.var.coefs$Sigma, x$phi.h ) )
  }
  if(do.lp) est.lp <- make.lp.irf( df, est.A, l.var.coefs, inf.name=x, this.fcast=fcast, irf.pds=irf.pds )
  
  #### 5. Make the IRFs ####
  est.bootstrap <- make.bootstrap( l.var.coefs, est.A$A, fcast.horiz=fcast.horiz,
                                   inf.name=x, n.pds = irf.pds, n.boot = n.boot,
                                   fcast.cumul=fcast.cumul, print.iter = FALSE,
                                   reorder.to.idx=reorder.to.idx)
  if(non.re.ci){
    l.est.bootstrap.non.re <- lapply( 1:n.non.re, function(i){
      message( paste0( '### Bootstrap, ', names(l.non.re.fns)[i], ' ###' ) )
      make.bootstrap( l.var.coefs, l.est.A.non.re[[i]]$A, fcast.horiz=fcast.horiz, inf.name=x,
                      n.pds = irf.pds, n.boot = n.boot, fcast.cumul=fcast.cumul,
                      make.phi.fn = l.non.re.fns[[i]], print.iter = FALSE, do.var.decomp = FALSE,
                      reorder.to.idx=reorder.to.idx )
    } ) %>% set_names( names(l.non.re.fns) )
  }
  
  # ... (o resto do seu script permanece igual) ...
  
  #### 7. Make some charts ####
  figure.location <- paste0(out_prefix, names(l.cases)[[i.case]] )
  df.irf.struct <- make.irf.struct( l.var.coefs, est.A$A, fcast.horiz, x, 
                                    irf.pds, fcast, fcast.cumul)
  df.quick.irfs <- lapply( 0:irf.pds, 
                           function(i) ((l.var.coefs$B %^% i)[1:n.vars,1:n.vars] %*% est.A$A) %>%
                             as.data.frame %>% rownames_to_column('variable') %>%
                             gather(shock, value, -variable) %>% mutate(pd=i) ) %>%
    reduce(full_join)
  
  ggplot( df.irf.struct %>% filter( grepl('Sentiment', shock ) ),
          aes(x=period, y=value, color=shock, shape=shock) ) +
    geom_hline(aes(yintercept=0)) +
    geom_line(lwd=1.2) +
    geom_point(size=3) +
    facet_wrap(~outcome, scales='free_y') +
    theme_minimal() + 
    theme(legend.position = 'bottom')
  ggsave(paste0('graphs/replication/',names(l.cases)[[i.case]], '/', i.case, 'case_quick_irfs.pdf'))
  df.var.decomp <- suppressMessages(
    lapply(1:length(est.var.decomp$decomp.ratio), function(i){
      est.var.decomp$decomp.ratio[[i]] %>%
        as.data.frame() %>%
        rownames_to_column("outcome") %>%
        gather(shock, value, -outcome) %>%
        mutate(shock = ifelse(str_detect(shock, "Fundamental"), "Fundamental", shock)) %>%
        group_by(outcome, shock) %>%
        summarise(value = sum(value), .groups = "drop") %>%
        mutate(horizon = i)
    }) %>% reduce(full_join) %>% as_tibble() %>%
      mutate(outcome.f = fct_relevel(outcome %>% as.factor, c(fcast,x, y ))) %>%
      arrange(outcome, shock, horizon)
  )
  ## 7.2 Now make lots more charts ##
  #source('code/R/sessions/charts.R')
  
  #### 8. Save the data ####
  save.string <- c( 'est.A', 'l.var.coefs', 'df', 'df.irf.struct', 'df.var.decomp', 
                    'est.bootstrap', if(non.re.ci) 'l.est.bootstrap.non.re', if(do.lp) 'est.lp' )
  save( list=save.string, file=paste0( 'model_solutions/', i.case, '_', init.yr, '_', freq, '_', l.var.coefs$lags, '.rdata' ) )
  
  df.shk.ar <- df.irf.struct %>% 
    filter( grepl( 'Sentiment', shock ), grepl('sentiment', outcome) ) %>% 
    group_by(shock,outcome) %>% 
    summarise( ar.1=ar(na.omit(value),order=1, aic=FALSE)$ar )
  # write_csv( df.shk.ar, file = paste0( 'graphs/', figure.location, '/shk_ar.csv') )
  
  struct.resid <- solve( est.A$A, l.var.coefs$var %>% resid() %>% t ) %>% t
  #df.struct.resid <- data.frame( date=(df %>% select(-value.deseas) %>% 
  #                                       spread( variable, value ) %>% 
  #                                       pull(date))[-(1:l.var.coefs$lags)],
  #                               struct.resid)
  # write_csv( df.struct.resid, file=paste0( this.case.name, '_resid.csv' ) )
} # fim do loop de casos

#### 9. Tidying up ####
source('code/R/sessions/comparison_charts.R')
source('code/R/sessions/simulation_nD.R')
source('code/R/sessions/narrative_compare.R')
source('code/R/sessions/renaming.R')
