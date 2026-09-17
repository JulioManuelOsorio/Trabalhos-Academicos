##########################################################################################
#
# charts_adaptado_BR_v2.R
# Versão adaptada para o projeto com dados do Brasil (Excel/CSV)
# (Usa as variáveis definidas pelo usuário)
#
##########################################################################################

library(tidyverse)
library(readxl)
library(zoo)    # rollmean
library(lubridate)

# -------------------------
# 0. Configurações / rótulos
# -------------------------
# Mapeie os nomes exatos das colunas do seu Excel (fornecidos por você)
# Variáveis que você declarou:
# consensus_expec_inflation, survey_expec_inflation, inflation, delta_er,
# exchange_rate, embi, selic_rate, ibc_br, commodities_index, fed_funds_rate, output_gap

setwd("C:\\Users\\julio\\OneDrive\\Desktop\\USP\\TCC")
load("model_solutions\\2_2005_m_2.rdata")

source('code/R/functions/make_var_fns.R')
source('code/R/functions/make_irf_fns.R')
source('code/R/functions/make_lp_fns.R')
source('code/R/functions/make_non_RE_var_fns.R')


var.labs <- c(
  consensus_expec_inflation = "Expectativa (consenso)",
  survey_expec_inflation    = "Expectativa (survey FGV)",
  inflation                 = "Inflação (anual)",
  delta_er                  = "Variação ER (delta)",
  exchange_rate             = "Taxa de câmbio (nível)",
  embi                      = "EMBI (spread)",
  selic_rate                = "SELIC (%)",
  ibc_br                    = "IBC-BR (atividade)",
  commodities_index         = "Índice Commodities",
  fed_funds_rate            = "Fed Funds Rate",
  output_gap                = "Hiato do Produto (output gap)",
  log_ip                    = "100*log(Produção Industrial)",
  log_dol                   = "100*log(USD/BRL)"
  
)

# Função helper para pegar rótulo ou retornar nome se não houver rótulo
lab_or_name <- function(x){
  # x pode ser vetor; retorna vetor de rótulos (ou nomes caso não tenha rótulo definido)
  sapply(x, function(xx){
    lab <- var.labs[xx]
    if(is.na(lab) || length(lab)==0) return(xx) else return(lab)
  }, USE.NAMES = FALSE)
}

# -------------------------
# 0.5 Variáveis de fallback (caso não existam no ambiente)
# -------------------------
# O master.R já define fcast, x, fcast.cumul, n.fcast, irf.pds, irf.plot.pds, figure.location etc.
# Aqui, se algum desses não existir, deixamos defaults sensatos para testar os gráficos.
if(!exists("fcast"))       fcast <- c("survey_expec_inflation")   # default: série de expectativa survey
if(!exists("x"))           x <- "inflation"
if(!exists("fcast.cumul")) fcast.cumul <- rep(TRUE, length(fcast))
if(!exists("n.fcast"))     n.fcast <- length(fcast)
if(!exists("irf.pds"))     irf.pds <- 36
if(!exists("irf.plot.pds")) irf.plot.pds <- ifelse(irf.pds==36, 24, irf.pds)
if(!exists("fcast.horiz")) fcast.horiz <- ifelse(irf.pds==36, 12, 4)

# plot_folder: usa figure.location (do master) se existir, senão fallback
if(exists("figure.location")){
  plot_folder <- paste0("graphs/", figure.location)
} else {
  plot_folder <- "graphs/replication"
}
if(!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)

# Annualize helper (se necessário)
annualize.vars <- c()  # se quiser anualizar respostas de inflação
#annualize <- function(x){ if(!exists("freq") || freq=="m") x*12 else x*4 }
irf.ylim <- if(!exists("freq") || freq=="m") c(-.8, .4) else c(-2, 2)


est.var.decomp <- var.decomp( est.A$A, l.var.coefs, n.pds = irf.pds, n.fcast=n.fcast)

# -------------------------
# 1. Reduced-form IRFs plot (robusto a faltas)
# -------------------------
if(exists("l.var.coefs") && exists("irf.pds")){
  df.irf.rf <- tryCatch(make.irf.rf(l.var.coefs, n.pds = irf.pds),
                        error = function(e){ message("make.irf.rf erro: ", e$message); NULL })
  if(!is.null(df.irf.rf)){
    gg.irfs.rf <- ggplot(df.irf.rf) +
      geom_path(aes(x=period, y=value), lwd = .9) +
      { if(exists("est.bootstrap") && !is.null(est.bootstrap$reduced.form))
        geom_ribbon(data = est.bootstrap$reduced.form, aes(x=period, ymin=q.05, ymax=q.95), alpha=.25)
      } +
      facet_grid(rows = vars(shock), cols = vars(outcome), labeller = labeller(shock = lab_or_name, outcome = lab_or_name)) +
      theme_minimal() +
      labs(title = paste0("Reduced-form IRFs — ", l.var.coefs$lags, " lags"), x = "Period", y = "") +
      xlim(c(0, irf.plot.pds))
    ggsave(filename = file.path(plot_folder, paste0("irf_rf_", ifelse(exists('init.yr'), init.yr, ''), "_", ifelse(exists('freq'),freq,''), "_", l.var.coefs$lags, "_lags.pdf")), plot = gg.irfs.rf, height = 10, width = 10)
  } else {
    message("df.irf.rf não disponível — pulando reduced-form IRF.")
  }
} else {
  message("l.var.coefs ou irf.pds não encontrado — pulando reduced-form IRF.")
}

# -------------------------
# 2. Structural IRFs
# -------------------------
# Helper: regex para capturar ambos "Sentiment 1" e "Sentiment #1" etc.
sent_regex <- function(i) paste0("Sentiment\\s*#?\\s*", i)
fund_regex <- function(i) paste0("Fundamental\\s*#?\\s*", i)

if(exists("l.var.coefs") && exists("est.A")){
    df.irf.struct <- tryCatch(make.irf.struct(l.var.coefs, est.A$A, fcast.horiz = fcast.horiz,
                                            inf.name = x, n.pds = irf.pds, fcast = fcast,
                                            fcast.cumul = fcast.cumul),
                            error = function(e){ message("make.irf.struct erro: ", e$message); NULL })
  if(!is.null(df.irf.struct)){
    # garantir coluna outcome.lab com rótulos
    df.irf.struct <- df.irf.struct %>% mutate(outcome.lab = lab_or_name(outcome))
    gg.irfs.struct <- ggplot(df.irf.struct) +
      geom_path(aes(x=period, y=value), lwd = .9) +
      { if(exists("est.bootstrap") && !is.null(est.bootstrap$structural))
        geom_ribbon(data = est.bootstrap$structural %>% mutate(outcome.lab = lab_or_name(outcome)),
                    aes(x=period, ymin=q.05, ymax=q.95), alpha = .25)
      } +
      facet_grid(rows = vars(shock), cols = vars(outcome.lab)) +
      theme_minimal() +
      xlim(c(0, irf.plot.pds)) +
      labs(title = paste0("Structural IRFs — ", l.var.coefs$lags, " lags"), x = "Período", y = "")
    #ggsave(filename = file.path(plot_folder, paste0("irf_struct_", ifelse(exists('init.yr'), init.yr, ''), "_", ifelse(exists('freq'),freq,''), "_", l.var.coefs$lags, "_lags.pdf")), plot = gg.irfs.struct, height = 10, width = 10)
  } else {
    message("df.irf.struct é NULL — pulando structural IRF.")
  }
} else {
  message("l.var.coefs ou est.A ausente — pulando structural IRF.")
}

# 2.1 Plots por forecast shock (non-fundamental)
if(exists("df.irf.struct") && !is.null(df.irf.struct)){

  for(this.fcast in seq_len(n.fcast)){
    # filtro robusto que casa com "Sentiment 1" e "Sentiment #1"
    pat <- sent_regex(this.fcast)
    df.irf.nf <- df.irf.struct %>%
      filter(str_detect(shock, regex(pat))) %>%
      mutate(outcome.lab = lab_or_name(outcome),
             value = ifelse(outcome %in% annualize.vars, annualize(value), value))
    if(nrow(df.irf.nf) == 0) {
      message("Nenhuma resposta encontrada para ", pat, " — pulando.")
      next
    }
    ribbon.df <- NULL
    if(exists("est.bootstrap") && !is.null(est.bootstrap$structural)){
      ribbon.df <- est.bootstrap$structural %>%
        filter(str_detect(shock, regex(pat))) %>%
        mutate(outcome.lab = lab_or_name(outcome),
               q.05 = ifelse(outcome %in% annualize.vars, annualize(q.05), q.05),
               q.95 = ifelse(outcome %in% annualize.vars, annualize(q.95), q.95))
    }
    gg.irfs.struct.nf <- ggplot(df.irf.nf) +
      geom_hline(yintercept = 0, lwd = .5) +
      geom_path(aes(x = period, y = value), lwd = .9) +
      { if(!is.null(ribbon.df)) geom_ribbon(data = ribbon.df, aes(x=period, ymin=q.05, ymax=q.95), alpha=.25) } +
      facet_wrap(~outcome.lab, scales = 'free_y') +
      theme_minimal() +
      xlab('Horizonte') + ylab('Resposta acumulada') +
      scale_x_continuous(breaks = seq(0, irf.plot.pds, 6), limits = c(0, irf.plot.pds))
    print(gg.irfs.struct.nf)
    ggsave(filename = file.path(plot_folder, paste0("irf_struct_nf_", this.fcast, "_", ifelse(exists('init.yr'), init.yr, ''), "_", ifelse(exists('freq'),freq,''), "_", l.var.coefs$lags, "_lags.pdf")), plot = gg.irfs.struct.nf, height = 8, width = 10)
  }
}

# 2.2 Non-RE comparisons (se disponíveis)
if(exists("l.est.A.non.re") && length(l.est.A.non.re) > 0 && exists("df.irf.struct")){
  for(this.fcast in seq_len(n.fcast)){
    pat <- sent_regex(this.fcast)
    df.irf.struct.non.re <- lapply(names(l.est.A.non.re), function(nn){
      tryCatch({
        make.irf.struct(l.var.coefs, l.est.A.non.re[[nn]]$A, fcast.horiz=fcast.horiz, inf.name=x, n.pds = irf.pds, fcast = fcast) %>%
          mutate(expectations = nn)
      }, error = function(e){ NULL })
    }) %>% discard(is.null) %>% reduce(full_join)
    if(nrow(df.irf.struct.non.re) == 0) next
    df.irf.struct.non.re <- df.irf.struct.non.re %>%
      filter(str_detect(shock, regex(pat))) %>%
      mutate(outcome.lab = lab_or_name(outcome),
             value = ifelse(outcome %in% annualize.vars, annualize(value), value)) %>%
      full_join(df.irf.struct %>% filter(str_detect(shock, regex(pat))) %>% mutate(expectations = "ER")) %>%
      mutate(expectations.f = as.factor(expectations))
    gg.irfs.non.re <- ggplot(df.irf.struct.non.re, aes(x=period, y=value, group=expectations.f, color=expectations.f, shape=expectations.f)) +
      geom_path(lwd = .9) + geom_hline(yintercept = 0, lwd = .5) + geom_point(size = 3) +
      facet_wrap(~outcome.lab, scales = 'free_y') +
      theme_minimal() + xlab("Horizonte") + ylab("Pontos percentuais") +
      theme(legend.position = "bottom") +
      scale_x_continuous(breaks = seq(0, irf.plot.pds, 6), limits = c(0, irf.plot.pds))
    print(gg.irfs.non.re)
    ggsave(filename = file.path(plot_folder, paste0("irf_struct_nf_", this.fcast, "_non_RE.pdf")), plot = gg.irfs.non.re, height = 8, width = 10)
  }
}

# -------------------------
# 3. Variance decomposition
# -------------------------
if(exists("est.var.decomp") && !is.null(est.var.decomp$decomp.ratio)){
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
      mutate(outcome.f = fct_relevel(outcome %>% as.factor, c(fcast,x,y))) %>%
      arrange(outcome, shock, horizon)
  )
  gg.var.decomp <- ggplot(df.var.decomp, aes(group=shock)) +
    geom_line(aes(x=horizon, y=value, color=shock), lwd = .9) +
    geom_point(aes(x=horizon, y=value, color=shock), size = 2) +
    { if(exists("est.bootstrap") && !is.null(est.bootstrap$var.decomp))
      geom_ribbon(data = est.bootstrap$var.decomp %>% mutate(outcome.f = outcome %>% as_factor()),
                  aes(x=horizon, ymin=q.05, ymax=q.95, fill=shock), alpha=.25)
    } +
    facet_wrap(~outcome.f, labeller = labeller(outcome.f = lab_or_name), scales = "free_y") +
    theme_minimal() + theme(legend.position="bottom") +
    labs(color="Shocks", shape="Shocks") + xlab("Horizonte") + guides(fill="none") +
    scale_x_continuous(breaks = seq(0, irf.plot.pds, 6), limits = c(0, irf.plot.pds)) +
    ylim(c(0,1)) + ylab("Fração da Variância")
  print(gg.var.decomp)
  ggsave(filename = file.path(plot_folder, paste0("var_decomp_", l.var.coefs$lags, "_lags.pdf")), plot = gg.var.decomp, height = 6, width = 10)
} else {
  message("est.var.decomp não disponível — pulando variance decomposition.")
}

# -------------------------
# 4. Série do choque de sentimento estimado
# -------------------------
if(exists("est.A") && exists("l.var.coefs")){
  struct.resid <- tryCatch({
    solve(est.A$A, (l.var.coefs$var %>% resid()) %>% t) %>% t
  }, error = function(e){ message("Erro em calcular struct.resid: ", e$message); NULL })
  if(!is.null(struct.resid)){
    # construir df com datas (usar df do ambiente global se existir)
    if(exists("df")){
      dates.all <- df %>% distinct(date) %>% arrange(date) %>% pull(date)
    } else if(!is.null(rownames(l.var.coefs$data))){
      dates.all <- as.Date(rownames(l.var.coefs$data))
    } else {
      dates.all <- seq_len(nrow(struct.resid))  # fallback
    }
    df.struct.resid <- data.frame(
      date = dates.all,
      rbind(matrix(NA, nrow = l.var.coefs$lags, ncol = ncol(est.A$A)), struct.resid)
    )
    names(df.struct.resid) <- c("date", paste0("Choque de Sentimento", 1:n.fcast), paste0("Choque Fundamental", 1:(ncol(est.A$A)-n.fcast)))
    df.struct.resid.long <- df.struct.resid %>%
      pivot_longer(cols = -date, names_to = "shock", values_to = "value") %>%
      group_by(shock) %>%
      mutate(value.12m = zoo::rollmean(value, k = 12, fill = NA, align = "right"))
    for(i in seq_len(n.fcast)){
      pattern <- paste0("Choque de Sentimento\\s*#?\\s*", i)
      sub.df <- df.struct.resid.long %>% filter(str_detect(shock, regex(pattern)))
      if(nrow(sub.df)==0) next
      gg.sentiment <- ggplot(sub.df, aes(x=date, y=value.12m)) +
        geom_line(lwd = .9) + geom_hline(yintercept = 0) +
        theme_minimal() + xlab("") + ylab(paste0("Sentiment #", i)) +
        theme(legend.position = "bottom")
      print(gg.sentiment)
      ggsave(filename = file.path(plot_folder, paste0("sentiment_", i, "_shk_", l.var.coefs$lags, "_lags.pdf")), plot = gg.sentiment, height = 5, width = 10)
    }
  } else {
    message("Não foi possível calcular struct.resid.")
  }
} else {
  message("est.A ou l.var.coefs ausente — pulando gráficos de sentimento.")
}

wri

df_struct <- df.struct.resid.long

df_struct <- df_struct %>%
  filter(shock == "Choque de Sentimento1")

df_struct

df_struct <- df_struct[3:nrow(df_struct),]

df_struct <- df_struct %>%
  mutate(value.12m = zoo::rollmean(value, k = 12, fill = NA, align = "right"))
ggplot(df_struct, aes(x = date)) +
  geom_line(aes(y = value, color = "Choque de Sentimento"), linewidth = 0.9, alpha = 0.7) +
  geom_line(aes(y = value.12m, color = "Média Móvel (12 meses)"), linewidth = 1.1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.9) +
  
  scale_color_manual(
    name = NULL,
    values = c("Choque de Sentimento" = "gray40", 
               "Média Móvel (12 meses)" = "black")
  ) +
  labs(
    x = "Data",
    y = "Sentimento"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),  # Aumenta o tamanho da fonte da legenda
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    legend.key.width = unit(1.5, "cm"),
    axis.title = element_text(size = 14),      # Títulos dos eixos (x e y)
    axis.text = element_text(size = 12),
    panel.border = element_blank()
  )
var(df_struct$value)
