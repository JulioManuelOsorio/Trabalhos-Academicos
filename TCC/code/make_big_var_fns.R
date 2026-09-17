##########################################################################################
#
# make_big_var_fns.R  (modificado: escolher lag sempre por SC / BIC + suporte a exógenas)
#
##########################################################################################

library(tidyverse)
library(BigVAR)
library(expm)
library(magrittr)
library(vars)   # necessário para VARselect e VAR
select <- dplyr::select

make.big.var <- function( df, value.name, var.name, fcast, inf, y, lag.max=24, m.Y=NA, exog=NA, ... ){
  # Faz VAR/BigVAR. Seleciona ordem por SC (BIC).
  # - Se exog == NA -> usa BigVAR (cv.BigVAR)
  # - Se exog fornecido (vector de nomes ou matrix/data.frame) -> usa vars::VAR(..., exogen = m.exog)
  #
  # Args:
  #  df, value.name, var.name, fcast, inf, y : como antes
  #  lag.max : máximo a considerar em VARselect
  #  m.Y : opcional, matrix já preparada (linhas = períodos, colnames = variáveis)
  #  exog: NA (padrão) OU vector de nomes das variáveis exógenas (presentes em df) OU matrix/data.frame alinhada
  #  ... : parâmetros adicionais para constructModel() quando BigVAR for usado
  #
  # Returns: lista com B, Sigma, n.vars, lags, mu, data, bigvar (ou NULL), var (obj VAR ou NULL),
  #          exog (matrix ou NULL), exog.names (vector ou NULL)
  
  # -----------------------------
  # montar m.Y se necessário
  # -----------------------------
  if( length(m.Y) == 1 && is.na(m.Y) ){
    tmp <- df %>%
      filter( !!as.symbol(var.name) %in% c( fcast, inf, y ) ) %>%
      select(date, !!as.symbol(var.name), !!as.symbol(value.name)) %>%
      spread( !!as.symbol(var.name), !!as.symbol(value.name) ) %>%
      arrange(date)
    m.Y.dates <- tmp$date
    m.Y <- tmp %>% select(-date) %>% as.matrix()
    rownames(m.Y) <- as.character(m.Y.dates)
  } else {
    if(is.data.frame(m.Y)) m.Y <- as.matrix(m.Y)
    if(!is.matrix(m.Y)) stop("make.big.var: m.Y precisa ser uma matrix (ou forneça df no formato esperado).")
    m.Y.dates <- if(!is.null(rownames(m.Y))) rownames(m.Y) else NULL
  }
  
  # -----------------------------
  # montar m.exog (se exog fornecido)
  # -----------------------------
  m.exog <- NULL; exog.names <- NULL
  if( !(length(exog)==1 && is.na(exog)) ){
    # exog foi fornecido: pode ser vector de nomes (character) ou matrix/data.frame
    if(is.character(exog)){
      tmp.exog <- df %>%
        filter( !!as.symbol(var.name) %in% exog ) %>%
        select(date, !!as.symbol(var.name), !!as.symbol(value.name)) %>%
        spread( !!as.symbol(var.name), !!as.symbol(value.name) ) %>%
        arrange(date)
      m.exog.dates <- tmp.exog$date
      m.exog.mat <- tmp.exog %>% select(-date) %>% as.matrix()
      rownames(m.exog.mat) <- as.character(m.exog.dates)
      # alinhar com m.Y (preservar ordem de m.Y)
      if(is.null(m.Y.dates)){
        # se m.Y n tem rownames, assumimos que tmp.exog tem as mesmas linhas na mesma ordem
        m.exog <- m.exog.mat
      } else {
        common.dates <- intersect(rownames(m.Y), rownames(m.exog.mat))
        if(length(common.dates) == 0) stop("make.big.var: nenhuma data em comum entre m.Y e exógenas.")
        m.Y <- m.Y[common.dates, , drop=FALSE]
        m.exog <- m.exog.mat[common.dates, , drop=FALSE]
        rownames(m.Y) <- rownames(m.exog)
      }
      exog.names <- colnames(m.exog)
    } else if(is.matrix(exog) || is.data.frame(exog)){
      if(is.data.frame(exog)) exog <- as.matrix(exog)
      m.exog <- exog
      exog.names <- colnames(m.exog)
      # alinhar se possível se rownames existem
      if(nrow(m.exog) != nrow(m.Y)){
        if(!is.null(rownames(m.exog)) && !is.null(rownames(m.Y))){
          common.dates <- intersect(rownames(m.Y), rownames(m.exog))
          if(length(common.dates)==0) stop("make.big.var: exog matrix não alinha com m.Y (linhas).")
          m.Y <- m.Y[common.dates, , drop=FALSE]
          m.exog <- m.exog[common.dates, , drop=FALSE]
        } else {
          stop("make.big.var: exog matrix tem número de linhas diferente de m.Y e não tem rownames para alinhar.")
        }
      }
    } else {
      stop("make.big.var: argumento 'exog' inválido. Deve ser NA, vetor de nomes (character) ou matrix/data.frame.")
    }
  } # fim montagem m.exog
  
  n.vars <- ncol(m.Y)
  y.means <- m.Y %>% colMeans()
  var.names <- colnames(m.Y)
  
  # -----------------------------
  # seleção de lag por SC (BIC)
  # -----------------------------
  vs <- VARselect(m.Y, lag.max = lag.max)
  chosen_p <- tryCatch({ as.integer(vs$selection['SC(n)']) }, error = function(e) NA)
  if(is.na(chosen_p) || chosen_p < 1){
    warning("make.big.var: seleção por SC falhou ou devolveu <1; usando p = 1")
    chosen_p <- 1
  }
  chosen_p <- min(chosen_p, lag.max)
  message("make.big.var: lag selecionado por SC (BIC) = ", chosen_p)
  
  # -----------------------------
  # Se houver exógenas -> estimar com vars::VAR (com exogen)
  # Caso contrário -> usar BigVAR (cv.BigVAR)
  # -----------------------------
  if(!is.null(m.exog)){
    message("make.big.var: exógenas detectadas -> usando vars::VAR com exogen (BigVAR ignorado para incluir exógenas).")
    var.est <- tryCatch({
      VAR(y = m.Y, p = chosen_p, type = 'none', exogen = m.exog)
    }, error = function(e){
      stop("make.big.var: erro ao estimar VAR com exógenas: ", e$message)
    })
    # extrair coeficientes (mesma lógica)
    m.B.1 <- var.est$varresult %>% sapply(., coef) %>% t
    if(chosen_p == 1){
      m.B <- m.B.1
    } else {
      lagmat <- cbind( diag(n.vars*(chosen_p-1)), matrix(0, n.vars*(chosen_p-1), n.vars) )
      rownames(lagmat) <- colnames(m.B.1)[1:(n.vars*(chosen_p-1))]
      m.B <- rbind(m.B.1, lagmat)
    }
    m.Sigma <- var.est %>% resid() %>% var()
    return( list( B = m.B,
                  Sigma = m.Sigma,
                  n.vars = n.vars,
                  lags = chosen_p,
                  mu = y.means,
                  data = m.Y,
                  bigvar = NULL,
                  var = var.est,
                  exog = m.exog,
                  exog.names = exog.names ) )
  } else {
    # sem exógenas: usar BigVAR (cv.BigVAR)
    message("make.big.var: sem exógenas -> usando BigVAR (cv.BigVAR).")
    mod <- constructModel( m.Y, p = chosen_p, intercept = FALSE, ... )
    res <- cv.BigVAR(mod)
    m.B.1 <- res@betaPred[,-1] %>%
      set_rownames(var.names) %>%
      set_colnames( paste0( rep(var.names, chosen_p), '.l', rep(1:chosen_p, each = n.vars) ) )
    if(chosen_p == 1){
      m.B <- m.B.1
    } else {
      lagmat <- cbind( diag(n.vars * (chosen_p - 1)), matrix(0, n.vars * (chosen_p - 1), n.vars) )
      rownames(lagmat) <- colnames(m.B.1)[1:(n.vars * (chosen_p - 1))]
      m.B <- rbind(m.B.1, lagmat)
    }
    m.Sigma <- (res@resids) %>% set_colnames(var.names) %>% var()
    return( list( B = m.B,
                  Sigma = m.Sigma,
                  n.vars = n.vars,
                  lags = chosen_p,
                  mu = y.means,
                  data = m.Y,
                  bigvar = res,
                  var = NULL,
                  exog = NULL,
                  exog.names = NULL ) )
  }
}
