##########################################################################################
#
# make_var_fns.R  (adaptado para suportar variáveis exógenas)
#
# Philip Barrett, adaptado por ChatGPT
#
##########################################################################################

library(tidyverse)
library(vars)
library(expm)   
library(lubridate)
library(magrittr)
select <- dplyr::select

# ----------------------------
# make.var
# ----------------------------
# Parâmetros principais:
#  df: data.frame long (date, variable, value.deseas)
#  value.name: 'value.deseas'
#  var.name: 'variable'
#  fcast, inf, y : vectores de nomes (nomes das variáveis contidas em df)
#  lags, lag.max, lag.select: seleção de defasagens (lag.select padrão 'SC')
#  m.Y: alternativa: matriz já montada de dados (linhas = períodos, colunas = variáveis)
#  exog: NA (padrão) OU um vetor de nomes das variáveis exógenas (presentes em df) OU uma matriz de exógenas já alinhada
# Retorno inclui elementos: B, Sigma, n.vars, lags, mu, data (m.Y), var (objeto VAR), exog (matriz ou NULL), exog.names
make.var <- function( df, value.name, var.name, fcast, inf, y,
                      lags=NA, lag.max=24, lag.select='SC', m.Y=NA, exog=NA ){
  # Se m.Y não for fornecido (é o NA escalar), constrói a matriz a partir de df
  if( length(m.Y) == 1 && is.na(m.Y) ){
    tmp <- df %>%
      filter( !!as.symbol(var.name) %in% c( fcast, inf, y ) ) %>%
      select(date, !!var.name, !!value.name) %>%
      spread( !!as.symbol(var.name), !!as.symbol(value.name) ) %>%
      arrange(date)
    # Guardamos as datas para alinhamento
    m.Y.dates <- tmp$date
    m.Y <- tmp %>% select(-date) %>% as.matrix()
    rownames(m.Y) <- as.character(m.Y.dates)
  } else {
    # Se m.Y foi fornecido diretamente, garantir que é matrix e sem "date"
    if(is.data.frame(m.Y)) m.Y <- as.matrix(m.Y)
    if(!is.matrix(m.Y)) stop("make.var: m.Y deve ser uma matrix quando providenciado.")
    # Se tiver rownames com datas, ok; senão não fazemos nada
    m.Y.dates <- if(!is.null(rownames(m.Y))) rownames(m.Y) else NULL
  }
  
  # Construir matriz de exógenas (m.exog) se exog fornecido como nomes
  m.exog <- NULL; exog.names <- NULL
  if( length(exog) == 1 && is.na(exog) ){
    # sem exógenas
    m.exog <- NULL
  } else if( is.character(exog) ){
    # exog é vetor de nomes em df; extraímos e alinhamos por date
    tmp.exog <- df %>%
      filter( !!as.symbol(var.name) %in% exog ) %>%
      select(date, !!var.name, !!value.name) %>%
      spread( !!as.symbol(var.name), !!as.symbol(value.name) ) %>%
      arrange(date)
    # alinhar datas entre m.Y e tmp.exog
    if(is.null(m.Y.dates)){
      # m.Y sem rownames: assumimos mesma ordem; definimos m.Y.dates a partir do tmp.exog se for necessário
      m.exog.dates <- tmp.exog$date
      m.exog <- tmp.exog %>% select(-date) %>% as.matrix()
      rownames(m.exog) <- as.character(m.exog.dates)
    } else {
      m.exog.dates <- tmp.exog$date
      m.exog.mat <- tmp.exog %>% select(-date) %>% as.matrix()
      rownames(m.exog.mat) <- as.character(m.exog.dates)
      # Interseção de datas
      common.dates <- intersect(rownames(m.Y), rownames(m.exog.mat))
      if(length(common.dates) == 0) stop("make.var: nenhuma data em comum entre variáveis endógenas e exógenas.")
      # Subset para as common.dates e preservar a ordem de m.Y
      m.Y <- m.Y[common.dates, , drop=FALSE]
      m.exog <- m.exog.mat[common.dates, , drop=FALSE]
      m.Y.dates <- common.dates
    }
    exog.names <- colnames(m.exog)
  } else if( is.matrix(exog) || is.data.frame(exog) ){
    # exog já é uma matrix/data.frame alinhada: converte e usa
    if(is.data.frame(exog)) exog <- as.matrix(exog)
    m.exog <- exog
    # tentar setar nomes
    exog.names <- colnames(m.exog)
    # se o número de linhas for igual às de m.Y, ok; se não, tentar alinhar se rownames existir
    if(nrow(m.exog) != nrow(m.Y)){
      if(!is.null(rownames(m.exog)) && !is.null(rownames(m.Y))){
        common.dates <- intersect(rownames(m.Y), rownames(m.exog))
        if(length(common.dates)==0) stop("make.var: exog matrix não alinha com m.Y (linhas).")
        m.Y <- m.Y[common.dates, , drop=FALSE]
        m.exog <- m.exog[common.dates, , drop=FALSE]
      } else {
        stop("make.var: exog matrix tem número de linhas diferente de m.Y e não tem rownames para alinhar.")
      }
    }
  } else {
    stop("make.var: argumento 'exog' inválido. Deve ser NA, vetor de nomes (character) ou matrix/data.frame.")
  }
  
  n.vars <- ncol(m.Y)
  y.means <- m.Y %>% colMeans()
  var.names <- colnames(m.Y)
  
  # -------------- Seleção de lags --------------
  if( lags %>% is.na ){
    # VARselect não aceita exog nos seus argumentos antigos; fazemos seleção com base em m.Y apenas.
    lagselect <- VARselect( na.omit(m.Y), lag.max = lag.max )
    # usa SC por padrão (ou outro se passado)
    lags <- lagselect$selection[paste0(lag.select, '(n)')]
  }
  p <- lags
  
  # -------------- Estima o VAR --------------
  # passamos exog se existir
  var.est <- tryCatch({
    if(!is.null(m.exog)) {
      VAR(y = na.omit(m.Y), p = p, type = 'none', exogen = m.exog)
    } else {
      VAR(y = na.omit(m.Y), p = p, type = 'none')
    }
  }, error = function(e) stop("Erro em VAR(): ", e$message))
  
  # Extrair matriz de coeficientes (mesma lógica original)
  m.B.1 <- var.est$varresult %>% sapply(., coef) %>% t
  if( p == 1 ){
    m.B <- m.B.1
  } else {
    lagmat <- cbind( diag(n.vars*(p-1)), matrix(0, n.vars*(p-1), n.vars) )
    rownames(lagmat) <- colnames(m.B.1)[1:(n.vars*(p-1))]
    m.B <- rbind(m.B.1, lagmat)
  }
  
  m.Sigma <- var.est %>% resid() %>% var()
  
  return( list( B = m.B,
                Sigma = m.Sigma,
                n.vars = n.vars,
                lags = p,
                mu = y.means,
                data = m.Y,
                var = var.est,
                exog = m.exog,
                exog.names = exog.names ) )
}

# ----------------------------
# make.phi.h  (sem alterações importantes - usa B)
# ----------------------------
make.phi.h <- function( l.var, horiz=12, inf.idx=2, cumul=TRUE ){
  n.vars <- l.var$n.vars
  phi.k <- if(cumul) lapply(1:horiz, function(x) l.var$B %^% x) %>% Reduce('+', .) else l.var$B %^% horiz
  out <- list()
  out$phi.h <- phi.k[inf.idx, 1:n.vars]
  names(out$phi.h) <- rownames(phi.k)[1:n.vars]
  return(out)
}

# ----------------------------
# make.var.fcast  (adaptado para exógenas / fallback ao companion)
# ----------------------------
# Versão melhorada:
# - se 'var.obj' (o objecto retornado por make.var) for passado e contiver var (objeto VAR),
#   usa predict(var.obj$var, n.ahead = n.fcast.pds, dumvar = exog.future) quando exog.future for fornecido.
# - se não houver var.obj/exog.future, faz o forecast clássico via matriz companion B (ignora exog).
make.var.fcast <- function( m.B, m.Y, n.fcast.pds=1, var.obj = NULL, exog.future = NULL ){
  # Se temos um objeto VAR e exog.future, use predict()
  if(!is.null(var.obj) && !is.null(var.obj$var) && !is.null(var.obj$exog) && !is.null(exog.future)){
    # exog.future deve ser matrix com nrow = n.fcast.pds e colnames na mesma ordem que var.obj$exog
    if(is.data.frame(exog.future)) exog.future <- as.matrix(exog.future)
    if(nrow(exog.future) != n.fcast.pds) stop("make.var.fcast: exog.future deve ter nrow == n.fcast.pds")
    # predict
    fc <- tryCatch({
      predict(var.obj$var, n.ahead = n.fcast.pds, dumvar = exog.future)
    }, error = function(e){
      warning("predict(VAR) com exógenas falhou: ", e$message, " — fallback para método por companion (ignora exog).")
      NULL
    })
    if(!is.null(fc)){
      # fc$fcst é uma lista; cada elemento é uma matrix com coluna 'fcst' (primeira)
      varnames <- names(fc$fcst)
      n.vars <- length(varnames)
      out.mat <- matrix(NA, nrow = n.fcast.pds, ncol = n.vars)
      colnames(out.mat) <- varnames
      for(i in seq_along(varnames)){
        # A coluna com forecasts (primeira coluna)
        out.mat[, i] <- fc$fcst[[ varnames[i] ]][, 1]
      }
      return(out.mat)
    }
    # se fc == NULL, fallback mais abaixo
  }
  
  # Fallback clássico (sem exógenas): companion matrix B
  n.vars <- ncol(m.Y)
  n.pds <- nrow(m.Y)
  n.lags <- (ncol(m.B)) / n.vars
  m.fcast <- (( m.B %^% n.fcast.pds ) %>% t()) %>% .[,1:n.vars]
  m.out <- matrix( NA, n.pds, n.vars ) %>% set_colnames( colnames(m.Y) )
  m.Y.use <- matrix(NA, n.pds, n.vars*n.lags ) %>% set_colnames( rownames(m.B) )
  for(j in 1:n.lags) m.Y.use[j:n.pds, (j-1)*n.vars + 1:n.vars] <- m.Y[1:(n.pds+1-j), ]
  for(i in n.lags:n.pds){
    m.out[i, ] <- m.Y.use[i,] %*% m.fcast
  }
  return(m.out)
}

# ----------------------------
# make.var.fcast.fit  (usa l.var quando possível)
# ----------------------------
make.var.fcast.fit <- function( l.var, n.fcast.pds = 1, cumul = FALSE, exog.future = NULL ){
  # Se cumul TRUE, soma os forecasts 1..n
  if(!is.null(l.var$var) && !is.null(l.var$exog) && !is.null(exog.future)){
    # tenta usar predict com exog.future
    fc.mat <- make.var.fcast(NULL, NULL, n.fcast.pds = n.fcast.pds, var.obj = l.var, exog.future = exog.future)
    if(cumul){
      out <- do.call(cbind, lapply(1:n.fcast.pds, function(i) fc.mat[i, , drop = FALSE])) # não faz exatamente soma por horizonte - fallback
      # Simpler: para cumul, soma ao longo dos horizontes:
      out <- sapply(1:ncol(fc.mat), function(j) sapply(1:nf <- nrow(fc.mat), function(h) sum(fc.mat[1:h, j, drop = TRUE]))) # but ugly
    } else {
      return(fc.mat)
    }
  } else {
    # Fallback: chamar make.var.fcast passando m.B/l.var$data
    return(make.var.fcast(l.var$B, l.var$data, n.fcast.pds = n.fcast.pds, var.obj = NULL, exog.future = NULL))
  }
}

# ----------------------------
# var.decomp, anc.analytic, A.reorder (mantidos, sem alteração)
# ----------------------------

var.decomp <- function( A, l.var, n.pds=20, n.fcast=1 ){
  l.var.decomp.ratio <- l.var.decomp <- list()
  n.vars <- l.var$n.vars
  p <- l.var$lags
  big.sigma <- l.var$Sigma
  B <- l.var$B
  
  A.sq.extended <- lapply( 1:n.vars, function(x){ 
    out <- 0*diag( n.vars * p )
    out[1:n.vars,1:n.vars] <- A[,x] %*% t(A[,x])
    return(out)
  } )
  big.sigma.extended <- 0*diag( n.vars * p ) ; big.sigma.extended[1:n.vars,1:n.vars] <- big.sigma
  big.sigma.denominator <- 0
  
  for( i in 1:n.pds ){
    if( i==1 ){
      l.var.decomp[[i]] <- lapply( 1:n.vars, function(x) ((B %^% i) %*% A.sq.extended[[x]] %*% t(B %^% i))[1:n.vars,1:n.vars] )
    } else {
      l.var.decomp[[i]] <- lapply( 1:n.vars, function(x) l.var.decomp[[i-1]][[x]] + 
                                     ((B %^% i) %*% A.sq.extended[[x]] %*% t(B %^% i))[1:n.vars,1:n.vars] )
    }
    big.sigma.denominator <- big.sigma.denominator + ((B %^% i) %*% big.sigma.extended %*% t(B %^% i))[1:n.vars,1:n.vars]
    l.var.decomp.ratio[[i]] <- sapply( 1:n.vars, function(x) diag(l.var.decomp[[i]][[x]]) / 
                                         diag( big.sigma.denominator) )
    shk.names <- c( paste0( 'Sentiment #', 1:n.fcast), paste0( 'Fundamental #', 1:(n.vars-n.fcast) ) )
    colnames(l.var.decomp.ratio[[i]]) <- shk.names 
  }
  
  return( list( var.decomp = l.var.decomp, decomp.ratio = l.var.decomp.ratio ) )
}

anc.analytic <- function( big.sigma, phi.h ){
  if(is.null(phi.h %>% nrow())) phi.h <- phi.h %>% t
  nn <- dim( big.sigma )[2]
  n.fcast <- phi.h %>% nrow()
  phi.h.fcast <- phi.h[1:n.fcast, 1:n.fcast]
  phi.h.c <- phi.h[,-(1:n.fcast)]
  if(is.null(phi.h.c %>% nrow())) phi.h.c <- phi.h.c %>% t
  gamma <- solve( diag(n.fcast) - phi.h.fcast, phi.h.c )
  
  sigma.11 <- big.sigma[ 1:n.fcast, 1:n.fcast  ]
  sigma.12 <- big.sigma[ 1:n.fcast, (1+n.fcast):nn  ]
  sigma.22 <- big.sigma[ (1+n.fcast):nn, (1+n.fcast):nn  ]
  
  V.S <- cbind( diag(n.fcast) - phi.h.fcast, -phi.h.c ) %*%
    big.sigma %*% t(cbind( diag(n.fcast) - phi.h.fcast, -phi.h.c ))
  m.Lambda <- chol(V.S) %>% t
  
  anc <- solve( m.Lambda, ( diag(n.fcast) - phi.h.fcast) %*% ( sigma.12 - gamma %*% sigma.22 ) ) %>% t
  anf <- solve( diag(n.fcast) - phi.h.fcast, m.Lambda ) + gamma %*% anc
  afc <- chol( sigma.22 - anc %*% t(anc) ) %>% t
  aff <- gamma %*% afc
  m.A <- cbind( rbind( anf, anc ),
                rbind( aff, afc ) )
  rownames(m.A) <- rownames(big.sigma)
  colnames(m.A) <- c( paste0( 'Sentiment ', 1:n.fcast),
                      paste0( 'Fundamental ', 1:(nn-n.fcast) ) )
  sigma.err <- max( abs( big.sigma - (m.A %*% t(m.A)) ) )
  status <- if(sigma.err<1e-08) 'success' else 'failure'
  return( list(A = m.A, sigma.err = sigma.err) )
}

A.reorder <- function( l.A, Sigma, n.fcast, fundamental.idx ){
  if(fundamental.idx==0) return(l.A)
  m.A <- l.A$A
  reorder.from.idx <- n.fcast + 1
  reorder.to.idx <- fundamental.idx
  m.trans <- diag(nrow(m.A))
  m.trans[reorder.from.idx,reorder.from.idx] <- m.trans[reorder.to.idx,reorder.to.idx] <- 0
  m.trans[reorder.from.idx,reorder.to.idx] <- m.trans[reorder.to.idx,reorder.from.idx] <- 1
  m.trans.sub <- m.trans[-(1:n.fcast),-(1:n.fcast)]
  m.Sigma <- Sigma
  m.Sigma.trans <- m.trans %*% m.Sigma %*% t(m.trans)
  m.A.trans <- m.trans %*% m.A
  m.A.trans.check <- m.A.trans %*% t(m.A.trans) - m.Sigma.trans
  m.D <- m.A.trans[-(1:n.fcast),-(1:n.fcast)]
  m.D.sq <- m.D %*% t(m.D)
  m.D.tilde <- chol( m.D.sq ) %>% t()
  m.D.tilde.check <- m.D.sq - ( m.D.tilde %*% t(m.D.tilde) )
  m.B <- m.A.trans[1,-1]
  m.B.tilde <- (solve( m.D.tilde ) %*% m.D %*% m.B) %>% t
  m.B.tilde.check <- m.B.tilde ^2 %>% sum() - m.B ^2 %>% sum()
  m.A.alt <- m.A.trans ; m.A.alt[-(1:n.fcast),-(1:n.fcast)] <- m.D.tilde ; m.A.alt[1,-1] <- m.B.tilde
  m.Sigma.trans.check <- m.A.alt %*% t(m.A.alt) - m.Sigma.trans
  m.A.final <- m.trans %*% m.A.alt %>% set_colnames(colnames(m.A)) %>% set_rownames(rownames(m.A))
  m.A.final.check <- (m.A.final %*% t(m.A.final) - m.Sigma) %>% abs %>% max
  out <- list( A = m.A.final, sigma.err = m.A.final.check )
  return(out)
}

# FIM
