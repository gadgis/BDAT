library(dplyr)
library(ggplot2)
library(tidyr)
# metrics_INLA_full <- readRDS("~/bdat/stage_bdat/output/metrics_INLA_full.rds")
# metrics_RF_full <- readRDS("~/bdat/stage_bdat/output/metrics_RF_full.rds")
# 
# 
# tt <- bind_rows(pred_INLA_full , pred_RF_full ) %>%
#   mutate( diffR = pred - obs ) %>%
#   group_by(sample_size, method, approach, type_val) %>%
#   reframe(res = mean(diffR))
# 
# 
# 
# results_summary <- bind_rows(metrics_INLA_full , metrics_RF_full) %>%
#   mutate(NSE = if_else(NSE<0,0,NSE)) %>%
#   group_by(sample_size, method, approach, type_val) %>%
#   summarise(across(c(ME, MAE, RMSE, r, r2, NSE, Cb), 
#                    mean, na.rm = TRUE),
#             .groups = "drop"
#   ) %>%
#   mutate(across(c(ME, MAE, RMSE, r, r2, NSE, Cb), round, digits = 4))
# 
# results_summaryV <- 
#   bind_rows(
#     metrics_INLA_full , metrics_RF_full
#     ) %>%
#   mutate(NSE = if_else(NSE<0,0,NSE)) %>%
#   
#   group_by(sample_size, method) %>%
#   summarise(across(c(ME, MAE, RMSE, r, r2, NSE, Cb),
#                    list(mean,sd), na.rm = TRUE),
#             .groups = "drop" 
#   ) 
# 
# 
# 
# pred_INLA_full <- readRDS("~/bdat/stage_bdat/output/pred_INLA_full.rds")
# pred_RF_full <- readRDS("~/bdat/stage_bdat/output/pred_RF_full.rds")

Myeval <- function(x, y){
  ME <- mean(y - x, na.rm = TRUE)
  RMSE <- sqrt(mean((y - x)^2, na.rm = TRUE))
  MAE <- mean(abs(y - x), na.rm = TRUE)
  r2 <- (cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2)
  SSE <- sum((y - x)^2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T))^2, na.rm = T)
  NSE <- 1 - SSE/SST
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  v <- sd(x, na.rm = T) / sd(y, na.rm = T)
  sx2 <- var(x, na.rm = T) * (length(x) - 1) / length(x)
  sy2 <- var(y, na.rm = T) * (length(y) - 1) / length(y)
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1/v + u^2)/2)^-1
  rCb <- r * Cb
  CCC <- rCb
  data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, CCC = CCC, Cb = Cb)
}

pattern <- "Xval_pH*"
pattern <- "Xval_arg*"

lf <- list.files(path = "output/", 
                 pattern = pattern,
                full.names = T
                 )



pred_INLA_full <- bind_rows(
  lapply(lf, readRDS)
  ) 


tt <- pred_INLA_full %>%
  pivot_longer(cols =pred:predKED ,
               names_to = "method",
               values_to = "pred",
  ) %>%
group_by(sample_size, method, approach, type_val,rep) %>%
  
  group_modify(~Myeval(.$pred,.$obs)) %>%
  ungroup()


results_summary <- tt %>%
  mutate(NSE = if_else(NSE<0,0,NSE)) %>%
  group_by(sample_size, method, approach, type_val) %>%
  summarise(across(c(ME, MAE, RMSE, r, r2, NSE, Cb), 
                   mean, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(across(c(ME, MAE, RMSE, r, r2, NSE, Cb), round, digits = 4))

results_summaryV <- tt %>%
  mutate(NSE = if_else(NSE<0,0,NSE)) %>%
  
  group_by(sample_size, method, approach, type_val) %>%
  dplyr::summarise(across(c(ME, MAE, RMSE, r, r2, NSE, Cb),
                   list(mean,sd), na.rm = TRUE),
            .groups = "drop" 
  ) 

results_summaryN <- tt %>%

  group_by(sample_size, method, approach, type_val) %>%
  dplyr::summarise(n = n(),
                   .groups = "drop" 
  ) %>%
  distinct(n)


# Visualisation des résultats
library(tidyr)
library(stringr)



results_summaryV %>%
  pivot_longer(ME_1:Cb_2, names_to = "Indice", values_to = "valeur") |>
  mutate(type = if_else(  grepl( "[1]", Indice ) ,"mean" , "sd" ) ,
         indice2 = str_split(Indice,"_",simplify = T)[,1]
         
  ) %>%
  
  pivot_wider(id_cols = !Indice,names_from = type, values_from = valeur) %>%
  filter(indice2 %in% c("RMSE","NSE","Cb")) %>%
  mutate(sample_size = as.numeric(sample_size)) %>%
  
  ggplot(
    aes(x = sample_size,
        y = mean
    )
  ) +
  
  geom_ribbon(aes(ymin = mean + 1.96 * sd / sqrt(results_summaryN$n),
                  ymax = mean - 1.96 * sd / sqrt(results_summaryN$n), 
                  fill = method,
                  linetype=type_val), alpha = 0.3)+
  
  geom_line( aes(
    color = method,
    linetype=type_val
  ) , size = 1) +
  
  geom_point(aes(color = method,
                 shape=type_val) , 
             size = 2) +
  
  geom_vline(xintercept = 2000, 
             
             color = "black", 
             linetype = "solid",
             size = 1
  ) +
  
  
  facet_grid(indice2~approach, scales= "free") +
  labs(
    x = "Taille de l'échantillon (calibration)",
    y = "Indice",
    color = "Méthode",
    title = paste("Évolution des indicateurs de la validation croisée (10) \nselon la taille d'échantillon")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )



