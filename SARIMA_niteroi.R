install.packages("pacman")

pacman::p_load(tidyverse, forecast, tseries, fable, readxl, zoo, vroom, foreign)




files_dengue<-paste0("C:/Relatorios DC e Nowcasting/DC e Nowcasting_CSVs/Dengue/",list.files(path="C:/Relatorios DC e Nowcasting/DC e Nowcasting_CSVs/Dengue/", pattern = ".csv$"))

dengue <- map_dfr(files_dengue, ~vroom(.,num_threads = 5, col_select=c(DT_SIN_PRI, SEM_PRI, SG_UF, ID_MN_RESI, CLASSI_FIN, DT_DIGITA)))

file_dengue_atual<-paste0("C:/Relatorios DC e Nowcasting/DC e Nowcasting_CSVs/Dengue atual/",list.files(path="C:/Relatorios DC e Nowcasting/DC e Nowcasting_CSVs/Dengue atual/", pattern = ".csv$"))

dengue_atual <- map_dfr(file_dengue_atual, ~vroom(.,num_threads = 5, col_select=c(DT_SIN_PRI, SEM_PRI, SG_UF, ID_MN_RESI, CLASSI_FIN, DT_DIGITA)))


dengue_atual$DT_SIN_PRI<-as.Date(as.character(dengue_atual$DT_SIN_PRI), format = "%Y-%m-%d")
dengue_atual$DT_DIGITA<-as.Date(as.character(dengue_atual$DT_DIGITA), format = "%Y-%m-%d")

dengue_total <- rbind(dengue, dengue_atual)

rm(dengue)
rm(dengue_atual)


dengue_niteroi <- dengue_total %>%
  filter(ID_MN_RESI == 330330) 

dengue_niteroi <- dengue_niteroi %>%
  mutate(ano_mes = floor_date(DT_SIN_PRI, "month"))

dengue_niteroi_mensal <- dengue_niteroi %>%
  count(ano_mes, name = "casos")


ts_niteroi <- ts(dengue_niteroi_mensal$casos,
                 start = c(2018, 1),   
                 frequency = 12)

plot(ts_niteroi, main = "Casos mensais de dengue em Niterói",
     ylab = "Casos", xlab = "Ano")

