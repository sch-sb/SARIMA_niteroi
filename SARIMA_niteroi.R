install.packages("pacman")

pacman::p_load(tidyverse, forecast, tseries, fable, readxl, zoo, vroom, foreign, scales, lubridate)

# ==== 1. População ====
POP <- read.csv2("C:/Users/sacha/Downloads/popBR2000-2024.long.csv") %>%
  mutate(MUNCOD = as.integer(MUNCOD),
         ANO    = as.integer(ANO),
         POP    = as.numeric(POP))

# ==== 2. Dengue ====
files_dengue <- paste0("F:/Relatorios DC e Nowcasting/Dengue/",
                       list.files(path="F:/Relatorios DC e Nowcasting/Dengue/", pattern = ".csv$"))

dengue <- map_dfr(files_dengue,
                  ~vroom(., num_threads = 5,
                         col_select=c(DT_SIN_PRI, SEM_PRI, SG_UF, ID_MN_RESI, CLASSI_FIN, DT_DIGITA)))

file_dengue_atual <- paste0("F:/Relatorios DC e Nowcasting/Dengue atual/",
                            list.files(path="F:/Relatorios DC e Nowcasting/Dengue atual/", pattern = ".csv$"))

dengue_atual <- map_dfr(file_dengue_atual,
                        ~vroom(., num_threads = 5,
                               col_select=c(DT_SIN_PRI, SEM_PRI, SG_UF, ID_MN_RESI, CLASSI_FIN, DT_DIGITA)))

dengue_atual$DT_SIN_PRI <- as.Date(as.character(dengue_atual$DT_SIN_PRI), format = "%Y-%m-%d")
dengue_atual$DT_DIGITA  <- as.Date(as.character(dengue_atual$DT_DIGITA), format = "%Y-%m-%d")

dengue_total <- rbind(dengue, dengue_atual)

rm(dengue); rm(dengue_atual)

dengue_total <- dengue_total %>%
  mutate(DT_SIN_PRI = as.Date(as.character(DT_SIN_PRI), format = "%Y-%m-%d")) %>%
  filter(DT_SIN_PRI >= as.Date("2008-01-01") & DT_SIN_PRI <= as.Date("2025-12-31"))

# Excluir os casos descartados
dengue_total$class <- as.character(as.factor(dengue_total$CLASSI_FIN))

dengue_total$class <- ifelse(is.na(dengue_total$class), 99, dengue_total$class)

dengue_total <- dengue_total %>%
  filter(class != 5) %>%
  filter(class != 13)

# ==== 3. Niterói ====
dengue_niteroi <- dengue_total %>%
  filter(ID_MN_RESI == 330330) %>%
  mutate(ano_mes = floor_date(DT_SIN_PRI, "month"),
         ano     = year(ano_mes)) %>%
  count(ano_mes, ano, name = "casos") %>%
  arrange(ano_mes) %>%
  left_join(POP %>% filter(MUNCOD == 330330) %>% select(ano = ANO, pop = POP),
            by = "ano") %>%
  mutate(incidencia = (casos / pop) * 100000)

# ==== 4. Itaboraí ====
dengue_itaborai <- dengue_total %>%
  filter(ID_MN_RESI == 330190) %>%
  mutate(ano_mes = floor_date(DT_SIN_PRI, "month"),
         ano     = year(ano_mes)) %>%
  count(ano_mes, ano, name = "casos") %>%
  arrange(ano_mes) %>%
  left_join(POP %>% filter(MUNCOD == 330190) %>% select(ano = ANO, pop = POP),
            by = "ano") %>%
  mutate(incidencia = (casos / pop) * 100000)

# ==== 5. Completar meses ====
seq_meses <- tibble(
  ano_mes = seq.Date(from = floor_date(min(dengue_total$DT_SIN_PRI, na.rm = TRUE), "month"),
                     to   = floor_date(max(dengue_total$DT_SIN_PRI, na.rm = TRUE), "month"),
                     by   = "1 month")
)

dengue_niteroi <- seq_meses %>%
  left_join(dengue_niteroi, by="ano_mes") %>%
  mutate(casos = replace_na(casos, 0),
         incidencia = replace_na(incidencia, 0))

dengue_itaborai <- seq_meses %>%
  left_join(dengue_itaborai, by="ano_mes") %>%
  mutate(casos = replace_na(casos, 0),
         incidencia = replace_na(incidencia, 0))

# ==== 6. Criar objetos ts ====
ts_full <- ts(dengue_niteroi$incidencia,
              start = c(year(min(dengue_niteroi$ano_mes)),
                        month(min(dengue_niteroi$ano_mes))),
              frequency = 12)

# ==== 7. Ajustar SARIMA ====

interv <- as.Date("2017-01-01")

df_pre <- dengue_niteroi %>% filter(ano_mes < interv)
df_pos <- dengue_niteroi %>% filter(ano_mes >= interv)

ts_pre <- ts(df_pre$incidencia,
             start = c(year(min(df_pre$ano_mes)), month(min(df_pre$ano_mes))),
             frequency = 12)

ts_pos <- ts(df_pos$incidencia,
             start = c(year(min(df_pos$ano_mes)), month(min(df_pos$ano_mes))),
             frequency = 12)

# Ajuste SARIMA apenas no pré-intervenção 
lambda <- BoxCox.lambda(ts_pre)
fit_pre <- auto.arima(ts_pre,
                      seasonal = TRUE,
                      lambda = lambda, biasadj = TRUE,
                      stepwise = FALSE, approximation = FALSE,
                      max.P = 2, max.Q = 2, max.p = 5, max.q = 5,
                      max.order = 10, D = 1)
summary(fit_pre)
checkresiduals(fit_pre)


# ==== 8. Contrafactual ====
h_cf <- min(48, length(ts_pos))
fc_cf <- forecast(fit_pre, h = h_cf)

fc_cf_tibble <- tibble(
  ano_mes = seq(from = interv, by = "month", length.out = h_cf),
  mean    = as.numeric(fc_cf$mean),
  lower   = as.numeric(fc_cf$lower[,2]),
  upper   = as.numeric(fc_cf$upper[,2])
)

impacto <- tibble(
  ano_mes     = fc_cf_tibble$ano_mes,
  obs_niteroi = as.numeric(ts_pos[1:h_cf]),
  cf_mean     = fc_cf_tibble$mean,
  cf_lwr      = fc_cf_tibble$lower,
  cf_upr      = fc_cf_tibble$upper,
  evitados    = cf_mean - obs_niteroi
)

# ==== 9. Dados para o gráfico ====
df_comp <- dengue_niteroi %>%
  select(ano_mes, incid_nit = incidencia) %>%
  left_join(dengue_itaborai %>% select(ano_mes, incid_sg = incidencia), by = "ano_mes")

df_plot  <- df_comp %>% mutate(tipo = "Observado")

df_cf <- tibble(
  ano_mes   = fc_cf_tibble$ano_mes,
  incid_nit = fc_cf_tibble$mean,
  incid_sg  = NA_real_,
  tipo      = "Contrafactual"
)

df_cf_ic <- tibble(
  ano_mes = fc_cf_tibble$ano_mes,
  lwr     = fc_cf_tibble$lower,
  upr     = fc_cf_tibble$upper
)

ultimo_mes_cf <- max(fc_cf_tibble$ano_mes)

# ==== 10. Gráfico final ====

ggplot() +
  geom_line(data = df_plot %>% filter(ano_mes <= ultimo_mes_cf),
            aes(x = ano_mes, y = incid_sg, color = "Itaboraí"), size = 0.8) +
  geom_line(data = df_plot %>% filter(ano_mes <= ultimo_mes_cf),
            aes(x = ano_mes, y = incid_nit, color = "Niterói observado"), size = 0.8) +
  geom_line(data = df_cf,
            aes(x = ano_mes, y = incid_nit, color = "Niterói contrafactual"),
            linetype = "dashed", size = 0.9) +
  # Área pós-intervenção pintada
  geom_rect(aes(xmin = interv, xmax = ultimo_mes_cf,
                ymin = -Inf, ymax = Inf),
            fill = "blue", alpha = 0.08, inherit.aes = FALSE) +
  geom_vline(xintercept = as.numeric(interv), linetype = "dashed", color = "black") +
  labs(title = "Impacto da Wolbachia",
       y = "Incidência de dengue por 100 mil habitantes",
       x = "Ano") +
  scale_color_manual(values = c("Niterói observado" = "blue",
                                "Niterói contrafactual" = "blue",
                                "Itaboraí" = "red")) +
  scale_x_date(limits = c(as.Date("2008-01-01"), ultimo_mes_cf)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())


#=============== RESULTADOS / DISCUSSÃO ==================


modelos <- list()
ic <- tibble()

for(p in 0:2) for(d in 0:1) for(q in 0:2) {
  for(P in 0:2) for(D in 0:1) for(Q in 0:2) {
    fit <- try(Arima(ts_pre, order=c(p,d,q), seasonal=c(P,D,Q)), silent=TRUE)
    if(!inherits(fit,"try-error")){
      ic <- bind_rows(ic, tibble(
        p=p,d=d,q=q,P=P,D=D,Q=Q,
        AIC=fit$aic, AICc=fit$aicc, BIC=fit$bic
      ))
    }
  }
}
top5 <- ic %>% arrange(AICc) %>% slice(1:5)
print(top5)


checkresiduals(fit_pre, lag=48)

# holdout
h <- 12
treino <- window(ts_pre, end=c(2019,12)) 
teste  <- window(ts_pre, start=c(2016,1))

fit_hold <- auto.arima(treino, seasonal=TRUE)
fc_hold  <- forecast(fit_hold, h=length(teste))

acc <- accuracy(fc_hold, teste)
acc_baseline <- accuracy(snaive(treino, h=length(teste)), teste)
print(acc)
print(acc_baseline)


fc_48 <- forecast(fit_pre, h=48)
autoplot(fc_48)




impacto_sum <- impacto %>%
  summarise(obs_total = sum(obs_niteroi),
            cf_total = sum(cf_mean),
            evitados = sum(evitados),
            reducao_pct = 100 * sum(evitados) / sum(cf_mean))
print(head(impacto,12)) 
print(impacto_sum)



