##Data preparation.

#Load the database and the libraries needed for the analysis.
load("C:/Users/Victoria/Desktop/TFM/Datos ATP/ATP.RData")
WOCAUSES <- read_excel("WOCAUSES.xlsx")
library(readxl)
library(dplyr)
library(nortest)
library(SmartEDA)
library(openxlsx)
library(ggplot2)
library(compareGroups)
library(mice)
library(MASS)
library(epiDisplay)
library(epitools)
library(epiR)
library(naniar)
library(howManyImputations)
library(ResourceSelection)
library(car)

##Selection, creation and recoding of variables.
df_ATP <- ATP %>%
  select(-Columna1, -Columna2, -sex)
df_ATP$WalkOver <- ifelse(df_ATP$Final_Partit == "W/O", TRUE, FALSE)
df_ATP$Default <- ifelse(df_ATP$Final_Partit == "Default", TRUE, FALSE)
df_ATP$Default<-factor(df_ATP$Default, levels = c("TRUE", "FALSE"), labels = c("YES", "NO"))
df_ATP$WalkOver<-factor(df_ATP$WalkOver, levels = c("TRUE", "FALSE"), labels = c("YES", "NO"))
df_ATP$Default_Type <- ifelse(df_ATP$Default == "YES" & grepl("-", df_ATP$score), "Durante", ifelse(df_ATP$Default == "YES", "Antes", NA))
df_ATP$year <- as.numeric(df_ATP$year)
df_ATP$tourney_level<- as.factor(df_ATP$tourney_level)
df_ATP$winner_hand<- as.factor(df_ATP$winner_hand)
df_ATP$surface<- as.factor(df_ATP$surface)
df_ATP$loser_hand<- as.factor(df_ATP$loser_hand)
df_ATP$round<- as.factor(df_ATP$round)
df_ATP$Final_Partit<- as.factor(df_ATP$Final_Partit)
df_ATP$sets<- as.factor(df_ATP$sets)
df_ATP$winner_rank<-as.numeric(df_ATP$winner_rank)
df_ATP$loser_rank<- as.numeric(df_ATP$loser_rank)
df_ATP$dif_rank<- as.numeric(df_ATP$dif_rank)
df_ATP$Default_Type<- as.factor(df_ATP$Default_Type)
str(df_ATP)

##Exploratory analysis.
##Analysis of the normality. 
qqnorm(df_ATP$winner_age, main = "Gr?fico QQ de Edad del Ganador")
qqline(df_ATP$winner_age)
lillie.test(x =  df_ATP$winner_age)

qqnorm(df_ATP$loser_age, main = "Gr?fico QQ de Edad del Perdedor")
qqline(df_ATP$loser_age)
lillie.test(x =  df_ATP$loser_age)

qqnorm(df_ATP$dif_age, main = "Gr?fico QQ de diferencia de edad")
qqline(df_ATP$dif_age)
lillie.test(x =  df_ATP$dif_age)

qqnorm(df_ATP$winner_rank, main = "Gr?fico QQ de Ranking del Ganador")
qqline(df_ATP$winner_rank)
lillie.test(x =  df_ATP$winner_rank)

qqnorm(df_ATP$loser_rank, main = "Gr?fico QQ de Ranking del perdedor")
qqline(df_ATP$loser_rank)
lillie.test(x =  df_ATP$loser_rank)

qqnorm(df_ATP$dif_rank, main = "Gr?fico QQ de diferencia de rankig")
qqline(df_ATP$dif_rank)
lillie.test(x =  df_ATP$dif_rank)

qqnorm(df_ATP$games, main = "Gr?fico QQ de juegos")
qqline(df_ATP$games)
lillie.test(x =  df_ATP$games)

#Creation of an exploratory report.
ExpReport(df_ATP, Template = NULL, Target = NULL, label = NULL, theme = "Default", op_file = "report.html", op_dir = getwd(), sc = NULL, sn = NULL, Rc = NULL)

Summary_Cat<- ExpCTable(df_ATP, Target = NULL, margin = 1, clim = 10, nlim = 10, round = 2, bin = 3, per = FALSE, weight = NULL)
write.xlsx(Summary_Cat, "C:/Users/Victoria/Desktop/TFM/Datos ATP/Summary_Cat.csv", rowNames = FALSE)

Summary_Num<-ExpNumStat(df_ATP,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2)
selected_columns <- c("Vname", "mean", "median", "SD", "IQR")
Summary_numeric <- Summary_Num[, selected_columns]
write.xlsx(Summary_Num, "C:/Users/Victoria/Desktop/TFM/Datos ATP/Summary_Num.csv", rowNames = FALSE)

Summary_Cat<- ExpCTable(WOCAUSES, Target = NULL, margin = 1, clim = 10, nlim = 10, round = 2, bin = 3, per = FALSE, weight = NULL)
write.xlsx(Summary_Cat, "C:/Users/Victoria/Desktop/TFM/Datos ATP/Summary_CatWO.csv", rowNames = FALSE)


##Epidemiological analysis. 

#Calculation of PI and 95%CI
prop_walkover <- sum(df_ATP$WalkOver == "YES") / nrow(df_ATP)
ci_walkover <- prop.test(sum(df_ATP$WalkOver == "YES"), nrow(df_ATP))$conf.int
prop_per_1000_walkover <- prop_walkover * 1000
ci_per_1000_walkover <- prop.test(sum(df_ATP$WalkOver == "YES"), nrow(df_ATP), conf.level = 0.95)$conf.int * 1000
prop_default <- sum(df_ATP$Default == "YES") / nrow(df_ATP)
ci_default <- prop.test(sum(df_ATP$Default == "YES"), nrow(df_ATP))$conf.int
prop_per_1000_default <- prop_default * 1000
ci_per_1000_default <- prop.test(sum(df_ATP$Default == "YES"), nrow(df_ATP), conf.level = 0.95)$conf.int * 1000
cat("Proporci?n de incidencia de WalkOver:", prop_walkover, "\n")
cat("Intervalo de confianza para WalkOver:", ci_walkover[1], "-", ci_walkover[2], "\n\n")
cat("Proporci?n de incidencia de WalkOver por cada 1000 partidos:", prop_per_1000_walkover, "%\n\n")
cat("Intervalo de confianza para WalkOver por cada 1000 partidos:", ci_per_1000_walkover[1], "-", ci_per_1000_walkover[2], "%\n\n")
cat("Proporci?n de incidencia de Default por cada 1000 partidos:", prop_per_1000_default, "%\n")
cat("Intervalo de confianza para Default por cada 1000 partidos:", ci_per_1000_default[1], "-", ci_per_1000_default[2], "%\n")
cat("Proporci?n de incidencia de Default:", prop_default, "\n")
cat("Intervalo de confianza para Default:", ci_default[1], "-", ci_default[2], "\n")

#Analysis of the time trend of the Default and WalkOver.

walkover_cases_by_year <- df_ATP %>%
  filter(WalkOver == "YES") %>%
  group_by(year) %>%
  summarize(WalkOver = n())

default_cases_by_year <- df_ATP %>%
  filter(Default == "YES") %>%
  group_by(year) %>%
  summarize(Default = n())

walkover_cases_by_year <- walkover_cases_by_year %>%
  mutate(Proportion_WalkOver = WalkOver / sum(walkover_cases_by_year$WalkOver))

default_cases_by_year <- default_cases_by_year %>%
  mutate(Proportion_Default = Default / sum(default_cases_by_year$Default))

ggplot(walkover_cases_by_year, aes(x = year, y = Proportion_WalkOver)) +
  geom_line() +
  geom_smooth(method = "loess", se = TRUE) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(1973, 2019, by = 2), limits = c(1973, 2019)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Tendencia temporal WalkOver", y = "PI", x = "Year")

ggplot(default_cases_by_year, aes(x = year, y = Proportion_Default)) +
  geom_line() +
  geom_smooth(method = "loess", se = TRUE) +  
  theme_minimal() +
  scale_x_continuous(breaks = seq(1973, 2019, by = 2), limits = c(1973, 2019)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = "Tendencia temporal Default", y = "PI", x = "Year")



##Bivariate analysis.

#Summary of covariates by groups
Summary_walkover<-compareGroups(WalkOver~tourney_level+surface+round+games+sets+winner_hand+loser_hand+winner_age+loser_age+dif_age+winner_rank+loser_rank+dif_rank, df_ATP, byrow=TRUE)
df_Summary_walkover<- createTable(Summary_walkover)
export2word(df_Summary_walkover, file='Summary_walkover.docx')

summary_default<-compareGroups(Default~tourney_level+surface+round+sets+games+winner_hand+loser_hand+winner_age+loser_age+dif_age+winner_rank+loser_rank+dif_rank, df_ATP, byrow = TRUE)
df_summary_default<- createTable(summary_default)
export2word(df_summary_default, file='summary_default.docx')

#Epidemiological analysis by groups 

 #Tourney level and WO
(tabla_tourney_walkover <- table(df_ATP$tourney_level, df_ATP$WalkOver))
valores_M_GS <- c( tabla_tourney_walkover["Masters", "YES"], tabla_tourney_walkover["Masters", "NO"], tabla_tourney_walkover["Grand Slams", "YES"], tabla_tourney_walkover["Grand Slams", "NO"])
epi.2by2(dat = valores_M_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

valores_250.500_GS <- c( tabla_tourney_walkover["250 or 500", "YES"], tabla_tourney_walkover["250 or 500", "NO"], tabla_tourney_walkover["Grand Slams", "YES"], tabla_tourney_walkover["Grand Slams", "NO"])
epi.2by2(dat = valores_250.500_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

valores_TF_GS<- c( tabla_tourney_walkover["Tour Finals", "YES"], tabla_tourney_walkover["Tour Finals", "NO"], tabla_tourney_walkover["Grand Slams", "YES"], tabla_tourney_walkover["Grand Slams", "NO"])
epi.2by2(dat = valores_TF_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

 #Surface and WO
(tabla_surface_walkover <- table(df_ATP$surface, df_ATP$WalkOver))

valores_G_C<- c(tabla_surface_walkover["Grass", "YES"], tabla_surface_walkover["Grass", "NO"],tabla_surface_walkover["Clay", "YES"], tabla_surface_walkover["Clay", "NO"])
epi.2by2(dat = valores_G_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

valores_H_C<- c(tabla_surface_walkover["Hard", "YES"], tabla_surface_walkover["Hard", "NO"], tabla_surface_walkover["Clay", "YES"], tabla_surface_walkover["Clay", "NO"])
epi.2by2(dat = valores_H_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

valores_C_C<- c(tabla_surface_walkover["Carpet", "YES"], tabla_surface_walkover["Carpet", "NO"], tabla_surface_walkover["Clay", "YES"], tabla_surface_walkover["Clay", "NO"])
epi.2by2(dat = valores_C_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

 #Sets and WO
(tabla_sets_walkover <- table(df_ATP$sets, df_ATP$WalkOver))
valores_3_5<- c(tabla_sets_walkover["3", "YES"], tabla_sets_walkover["3", "NO"], tabla_sets_walkover["5", "YES"], tabla_sets_walkover["5", "NO"])
epi.2by2(dat = valores_3_5, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")
 
#Round and WO
(tabla_round_walkover <- table(df_ATP$round, df_ATP$WalkOver))
valores_F_Q<- c(tabla_round_walkover["Final", "YES"], tabla_round_walkover["Final", "NO"], tabla_round_walkover["Qualifying", "YES"], tabla_round_walkover["Qualifying", "NO"])
epi.2by2(dat = valores_F_Q, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

valores_P_Q<- c(tabla_round_walkover["Preliminary", "YES"], tabla_round_walkover["Preliminary", "NO"], tabla_round_walkover["Qualifying", "YES"], tabla_round_walkover["Qualifying", "NO"])
epi.2by2(dat = valores_P_Q, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")


 #Tourney level and def
tabla_tourney_Default <- table(df_ATP$tourney_level, df_ATP$Default)
Valores_M_GS <- c(tabla_tourney_Default["Masters", "YES"], tabla_tourney_Default["Masters", "NO"], tabla_tourney_Default["Grand Slams", "YES"], tabla_tourney_Default["Grand Slams", "NO"])
epi.2by2(dat = Valores_M_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

Valores_250.500_GS <- c(tabla_tourney_Default["250 or 500", "YES"], tabla_tourney_Default["250 or 500", "NO"], tabla_tourney_Default["Grand Slams", "YES"], tabla_tourney_Default["Grand Slams", "NO"])
epi.2by2(dat = Valores_250.500_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

Valores_TF_GS<- c(tabla_tourney_Default["Tour Finals", "YES"], tabla_tourney_Default["Tour Finals", "NO"], tabla_tourney_Default["Grand Slams", "YES"], tabla_tourney_Default["Grand Slams", "NO"])
epi.2by2(dat = Valores_TF_GS, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

 #Surface and def
(tabla_surface_Default <- table(df_ATP$surface, df_ATP$Default))

Valores_G_C<- c(tabla_surface_Default["Grass", "YES"], tabla_surface_Default["Grass", "NO"], tabla_surface_Default["Clay", "YES"], tabla_surface_Default["Clay", "NO"])

epi.2by2(dat = Valores_G_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

Valores_H_C<- c(tabla_surface_Default["Hard", "YES"], tabla_surface_Default["Hard", "NO"], tabla_surface_Default["Clay", "YES"], tabla_surface_Default["Clay", "NO"])

epi.2by2(dat = Valores_H_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

Valores_C_C<- c( tabla_surface_Default["Carpet", "YES"], tabla_surface_Default["Carpet", "NO"], tabla_surface_Default["Clay", "YES"], tabla_surface_Default["Clay", "NO"])

epi.2by2(dat = Valores_C_C, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

 #Sets and def
(tabla_sets_Default <- table(df_ATP$sets, df_ATP$Default))

epi.2by2(dat = tabla_sets_Default, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

 #Round and def
(tabla_round_Default <- table(df_ATP$round, df_ATP$Default))

Valores_F_Q<- c(
  tabla_round_Default["Final", "YES"],
  tabla_round_Default["Final", "NO"],
  tabla_round_Default["Qualifying", "YES"],
  tabla_round_Default["Qualifying", "NO"]
)

epi.2by2(dat = Valores_F_Q, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

Valores_P_Q<- c(
  tabla_round_Default["Preliminary", "YES"],
  tabla_round_Default["Preliminary", "NO"],
  tabla_round_Default["Qualifying", "YES"],
  tabla_round_Default["Qualifying", "NO"]
)

epi.2by2(dat = Valores_P_Q, method = "cohort.count", digits = 2,
         conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")



##Logistic Regression Analysis
#Imputation of missing data.
variables_regression <- c("tourney_level", "surface", "sets", "round", "year", "dif_age", "dif_rank","WalkOver", "Default")

df_regression <- df_ATP[variables_regression]

pct_miss(df_regression)
pct_miss_case(df_regression)
gg_miss_var(df_regression)
(prop_miss_case(df_regression))

imputed<- mice(df_regression, method = "pmm", m = 20, maxit = 10, seed = 123)
model_impt <- with(imputed, glm(WalkOver ~ tourney_level + surface + sets + round + dif_age + dif_rank + year,  family = "binomial"))
(how_many_imputations(model_impt))
imputed<- mice(df_regression, method = "pmm", m = 29, maxit = 10, seed = 123)
df_regression <- complete(imputed)

#Variable response transformation and selection of reference dummmies
df_regression$WalkOver <- ifelse(df_regression$WalkOver == "YES", 1, 0)
df_regression$Default <- ifelse(df_regression$Default == "YES", 1, 0)
df_regression$tourney_level <- relevel(df_regression$tourney_level, ref = "Grand Slams")
df_regression$surface <- relevel(df_regression$surface, ref = "Clay")
df_regression$sets <- relevel(df_regression$sets, ref = "5")
df_regression$round <- relevel(df_regression$round, ref = "Qualifying")

#Logistic regression model for WalKover.
modelo_wo_inicial <- glm(WalkOver ~ tourney_level + surface + round + sets + dif_age + dif_rank + year,
                              data = df_regression, family = binomial)
summary(modelo_wo_inicial)

modelo_wo_forward <- stepAIC(modelo_wo_inicial, direction = "forward", trace = FALSE)

summary(modelo_wo_forward)

modelo_wo_backward <- stepAIC(modelo_wo_inicial, direction = "backward", trace = FALSE)

summary(modelo_wo_backward)

logistic.display(modelo_wo_forward, crude.p.value=FALSE, decimal = 2)
logistic.display(modelo_wo_backward, crude.p.value=FALSE, decimal = 2)

coeficientes <- summary(modelo_wo_backward)$coefficients
or <- exp(coeficientes[, "Estimate"])
lower_ci <- exp(coeficientes[, "Estimate"] - 1.96 * coeficientes[, "Std. Error"])
upper_ci <- exp(coeficientes[, "Estimate"] + 1.96 * coeficientes[, "Std. Error"])

forest_data <- data.frame(Variable = rownames(coeficientes), OR = or, LowerCI = lower_ci, UpperCI = upper_ci)
forest_data <- forest_data[-1, ]

ggplot(forest_data, aes(x = Variable, y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  coord_flip() +  # Voltea el gráfico para una mejor visualización
  xlab("Variables") +
  ylab("Odds Ratios (IC 95%)") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  geom_text(aes(label = paste("OR:", round(OR, 2), "\n[", round(LowerCI, 2), "-", round(UpperCI, 2), "]")), 
            hjust = 1.1, size = 2, position = position_dodge(width = 0.8))


aic_inicial <- AIC(modelo_wo_inicial)
aic_forward <- AIC(modelo_wo_forward)
aic_backward <- AIC(modelo_wo_backward)
comparacion_aic <- data.frame( Modelo = c("Inicial", "Forward", "Backward"), AIC = c(aic_inicial, aic_forward, aic_backward))
print(comparacion_aic)

hoslem_result_wo_for <- hoslem.test(df_regression$WalkOver,fitted(modelo_wo_forward), g = 10)
print(hoslem_result_wo_for)
hoslem_result_wo_bac<- hoslem.test(df_regression$WalkOver,fitted(modelo_wo_backward), g = 10)
print(hoslem_result_wo_bac)

plot(modelo_wo_forward, wich=1)
plot(modelo_wo_backward, which=1)


#Logistic regression model for Default.
modelo_default_inicial <- glm(Default ~ tourney_level + surface + round + sets + dif_age + dif_rank + year,
                             data = df_regression, family = binomial)
summary(modelo_default_inicial)
modelo_default_forward <- stepAIC(modelo_default_inicial, direction = "forward", trace = FALSE)
summary(modelo_default_forward)
modelo_default_backward <- stepAIC(modelo_default_inicial, direction = "backward", trace = FALSE)
summary(modelo_default_backward)
logistic.display(modelo_default_forward, crude.p.value=TRUE, decimal = 2)
logistic.display(modelo_default_backward, crude.p.value=TRUE, decimal = 2)

aic_inicial <- AIC(modelo_default_inicial)
aic_forward <- AIC(modelo_default_forward)
aic_backward <- AIC(modelo_default_backward)
comparacion_aic <- data.frame( Modelo = c("Inicial", "Forward", "Backward"), AIC = c(aic_inicial, aic_forward, aic_backward))
print(comparacion_aic)

coeficientes <- summary(modelo_default_backward)$coefficients
or <- exp(coeficientes[, "Estimate"])
lower_ci <- exp(coeficientes[, "Estimate"] - 1.96 * coeficientes[, "Std. Error"])
upper_ci <- exp(coeficientes[, "Estimate"] + 1.96 * coeficientes[, "Std. Error"])

forest_data <- data.frame(Variable = rownames(coeficientes), OR = or, LowerCI = lower_ci, UpperCI = upper_ci)
forest_data <- forest_data[-1, ]

ggplot(forest_data, aes(x = Variable, y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  coord_flip() +  
  xlab("Variables") +
  ylab("Odds Ratios (IC 95%)") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") + 
  geom_text(aes(label = sprintf("OR: %.2f\n[%.2f, %.2f]", OR, LowerCI, UpperCI)),
            hjust = -0.1, size = 3, position = position_nudge(y = 0.1))

hoslem_result_de_for <- hoslem.test(df_regression$Default,fitted(modelo_default_forward), g = 10)
print(hoslem_result_de_for)
hoslem_result_de_bac<- hoslem.test(df_regression$Default,fitted(modelo_default_backward), g = 10)
print(hoslem_result_de_bac)


plot(modelo_default_forward, wich=1)
plot(modelo_default_backward, which=1)

