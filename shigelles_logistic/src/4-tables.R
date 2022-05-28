###################
# Tableau article #
###################





#==================================================
#---------------------------------
# packages
#---------------------------------

#a ne lancer qu'une fois
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("tableone")
# install.packages("desctable")
# install.packages("stringr")

# alancer a chaque session
library(dplyr)
library(tidyr)
library(stringr)
# library(tableone)
# library(desctable)

#---------------------------------
# objets
#---------------------------------
source("src/02_data_management.R")

all_ATB <- c("AMXAMP_SIR", "C3G_SIR", "NAL_SIR", "CIP_SIR", "STR_SIR", "AKN_SIR", "GEN_SIR", "SMX_SIR", "TMP_SIR",
             "SXT_SIR", "CHL_SIR", "TET_SIR", "AZM_SIR")

#---------------------------------
#fonctions
#---------------------------------

#Calcule la prevalence de la resistance et son IC95%
get_prev <- function(mydata,myvar){
  # mydata <- dAZM
  # mydata <- dC3G
  # myvar <- "age_cut"
  # myvar <- "age_cut"
  myvec <- mydata[ ,myvar]
  myvec <- as.factor(myvec)
  alllev <- levels(myvec)
  #mylev <- alllev[1]
  .l <- lapply(alllev, function(mylev){
    N <- nrow(mydata[mydata[ ,myvar] == mylev & !is.na(mydata[ ,myvar]), ])
    M <- nrow(mydata[mydata[ ,myvar] == mylev & mydata$status==1 & !is.na(mydata[ ,myvar]) & !is.na(mydata$status), ])
    btest <- binom.test(x = M, n = N)
    p0 <- as.numeric(btest$estimate)
    IC <- as.numeric(btest$conf.int)
    IC <- round(c(p0, IC)*100, 2) # en %
    return(data.frame(variables = myvar, levels = mylev, N=N, n=M, perc_resistance=IC[1], IC_l=IC[2], IC_u=IC[3]))  
  })
  .l <- do.call(rbind, .l)
}

#Calcule l'OR de chaque variable avec IC 95% et pvalue (test Z(binaire ou quanti) ou Likelihood ratio test(quali a plus de 2 classes))
get_univar2 <- function(my_data, my_outcome, my_var){
  #my_data = d2
  #my_var <- "Age_cat"
  #my_var <- factors_notif[2]
  #my_outcome <- "HIV"
  my_data <- as.data.frame(my_data) #sinon probleme pour levels(factor(my_data[ ,my_var]))
  my_data <- na.omit(my_data[ ,c(my_outcome, my_var)])
  print(my_var)
  
  #modele
  my_formula <- as.formula (paste0(my_outcome," ~ ", my_var))
  mod1 <- glm(my_formula, data = my_data, family = "binomial")
  g <- summary(mod1)
  
  #récupérer la pvalue (selon variable catégorielle ou non)
  if(length(levels(factor(my_data[ ,my_var])))>2){
    my_formula <- as.formula (paste0(my_outcome," ~ 1"))
    mod0 <- glm(my_formula, data = my_data, family = "binomial")
    pval <- anova(mod1, mod0, test="LRT")
    pval <- pval$`Pr(>Chi)`[2]
  } else {
    tab <- g$coefficients  
    pval <- tab[grepl(my_var, rownames(tab)), "Pr(>|z|)"]
  }
  
  #récupérer OR et IC
  pval <- round(pval, 3)
  OR <- exp(coef(mod1))
  CI <- exp(confint(mod1))
  ORCI <- round(cbind(OR, CI), 2)
  
  #organiser le tableau
  res <- data.frame(outcome = my_outcome, variable = my_var, levels = NA, ORCI, pval, stringsAsFactors = FALSE)
  res <- res[grepl(my_var, rownames(res)), ]
  #levels <- gsub(my_var, "", rownames(res)) #commenté 20190108
  levels <- str_replace(rownames(res), my_var, "") #ajout 20190108
  levels <- if(any(levels != "")) levels else "1"
  res$levels <- levels
  names(res) <- c("outcome", "variables", "levels", "OR", "IC_l", "IC_u", "pvalue")
  return(res)
}

#calcule la pvalue pour une variable donn?e d'un mod?le multivari?

#Likelihood ratio test(quali a plus de 2 classes) (jai vu que c'?tait pareil que test Z, ? v?rifier)
get_pval_multivar <- function(mod1 = NULL, mod0 = NULL, var_to_test = NULL, all_var = NULL, my_outcome =NULL, my_data =NULL){
  
  if (is.null(mod1) & is.null(mod0)){
    other_var <- all_var[!all_var %in% var_to_test]
    my_data <- na.omit(my_data[ ,c(var_to_test, other_var, my_outcome)])
    #mod1
    my_var_mod1 <- paste(c(var_to_test, other_var), collapse = " + ")
    my_formula <- as.formula (paste0(my_outcome," ~ ", my_var_mod1))
    mod1 <- glm(my_formula, data = my_data, family = "binomial")
    
    #mod0
    if(length(other_var)==0) {
      my_formula <- as.formula (paste0(my_outcome," ~ 1"))
    } else {
      my_var_mod0 <- paste(c(other_var), collapse = " + ")
      my_formula <- as.formula (paste0(my_outcome," ~ ", my_var_mod0))
    }
    mod0 <- glm(my_formula, data = my_data, family = "binomial")
  }
  
  pval <- anova(mod1, mod0, test="LRT")
  pval <- pval$`Pr(>Chi)`[2]
  return(pval)
}

OR_multivar <- function(all_var, my_outcome, my_data){
  myvarm_vec <- paste(all_var, collapse = "+")
  modm <- glm(formula(paste0(my_outcome, " ~ ", myvarm_vec)), family = "binomial", data = my_data)
  OR <- exp(coef(modm))
  CI <- exp(confint(modm))
  ORCI <- round(cbind(OR, CI), 2)
  return(ORCI)
}

my_sub <- function(string, pattern) {str_replace(string, pattern, "")}

#calcul effectif et pourcentage
# var_qual <- function(name_var, mydata, myuse = "no"){
#   #TO TEST
#   #mydata <- CIPR.df
#   #name_var <- "age_cut"
#   #myuse <- "a"
#   
#   var <- mydata[ , name_var]
#   
#   # Calcul de N avec ou sans Na selon myuse
#   freq <- table(var, useNA = myuse)
#   my_levels <- names(freq)
#   my_freq <- as.numeric(freq)
#   
#   # # Calcul des pourcentages
#   # perc <- as.numeric(round(prop.table(freq)*100,2))
#   
#   # Calcul des pourcentages 28022019 : round 1 et pas les missing
#   freqnoNA <- table(var, useNA = "no") #calcul des pourcentages sans les missing
#   perc <- as.numeric(round(prop.table(freqnoNA)*100,1)) #round 1
#   
#   data.frame(var = name_var, level = my_levels, freq = my_freq, percent = perc, stringsAsFactors = FALSE)
# }

var_qual <- function(name_var, mydata, myuse = "no"){
  #TO TEST
  #mydata <- CIPR.df
  #name_var <- "age_cut"
  #myuse <- "a"
  
  var <- mydata[ , name_var]
  
  # Calcul de N avec ou sans Na selon myuse
  freq <- table(var, useNA = myuse)
  my_levels <- names(freq)
  #my_freq <- as.numeric(freq)
  my_freq <- data.frame(freq)
  
  # # Calcul des pourcentages
  # perc <- as.numeric(round(prop.table(freq)*100,2))
  
  # Calcul des pourcentages 28022019 : round 1 et pas les missing
  freqnoNA <- table(var, useNA = "no") #calcul des pourcentages sans les missing
  percnoNA <- prop.table(freqnoNA)
  percnoNA <- data.frame(percnoNA)
  percnoNA$percent <- round(percnoNA$Freq*100, 1)
  percnoNA$Freq <- NULL
  #perc <- as.numeric(round(prop.table(freqnoNA)*100,1)) #round 1
  
  var_des_df <- full_join(my_freq, percnoNA)
  var_des_df$level <- as.character(var_des_df$var)
  var_des_df$var <- as.character(name_var)
  names(var_des_df) <- tolower(names(var_des_df))
  return(var_des_df[ , c("var", "level", "freq", "percent")])
  #data.frame(var = name_var, level = my_levels, freq = my_freq, percent = perc, stringsAsFactors = FALSE)
}


describe_shigelle <- function(my_data, namecol = "freq"){
  
  #TO TEST
  #my_data <- CIPR.df
  
  #tableau avec les variables et s'il faut compter ou non les missing 
  df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "age_sexe","espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                       Nado = c("a", "a", "a", "a","a", "no", "a", "a", "no", "a"), stringsAsFactors = FALSE)
  
 
  #calcul effectif et pourcentage pour chaque variable
  res_vars <- lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    res_var <- var_qual(name_var=name, mydata=my_data, myuse=Nado)
    return(res_var)
  })
  #regrouper variables
  res_vars <- bind_rows(res_vars)
 
  #calcul des effectifs et pourcentage
  myres.df <-  res_vars %>%
    #mise en forme
    #group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>%
    mutate(
      #variable que sur la première ligne
      #var = ifelse(nr !=1, NA, var),
      #N(%)
      #col_freq_p = paste0(freq, " (", percent, ")"),#commenté 20190228
      col_freq_p = ifelse(!is.na(level), paste0(freq, " (", percent, ")"), freq),
      #Level missing au lieu de NA
      level = ifelse(is.na(level), "missing", level),
      #Je ne garde que var level et freq metrop
      freq = NULL, 
      percent = NULL
      #, nr = NULL
    ) %>% data.frame

  # Renomme la colonne de fréquence
  names(myres.df)[names(myres.df) %in% "col_freq_p"] <- namecol
  #names(myres.df)[names(myres.df) %in% "freq"] <- paste0(namecol,"_freq")
  
  return(myres.df)
}

#---------------------------------
# Tableau 1 : population
#---------------------------------

df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "age_sexe","espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                     Nado = c("a", "a", "a", "a", "a", "no", "a", "a", "no", "a"), stringsAsFactors = FALSE)

#DOMTOM
my_data <- d%>% filter(metrop == FALSE) 
all.res <- 
  #calcul effectif et pourcentage pour chaque variable
  lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    var_qual(name, my_data, Nado)
  }) %>% 
  #regrouper variables
  bind_rows() %>% 
  #mise en forme
  group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>% 
  mutate(
    #variable que sur la première ligne
    var = ifelse(nr !=1, NA, var),
    #N(%)
    #freq_dom = paste0(freq, " (", percent, ")"), #commenté 20190228
    freq_dom = ifelse(!is.na(level), paste0(freq, " (", percent, ")"), freq),
    #Level missing au lieu de NA
    level = ifelse(is.na(level), "missing", level),
    #Je ne garde que var level et freq metrop
    freq = NULL, percent = NULL, nr = NULL)

#METROPOLE
my_data <- d%>% filter(metrop == TRUE) 
all.res2 <- 
  #calcul effectif et pourcentage pour chaque variable
  lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    var_qual(name, my_data, Nado)
  }) %>% 
  #regrouper variables
  bind_rows() %>% 
  #mise en forme
  mutate(freq = paste0(freq, " (", percent, ")"))

#Je n'integre que freq metrop au tableau précédent
all.res$freq_metrop <- all.res2$freq 


names_corr <- data.frame(var = df_var$variables, all = c("N", "Age (%)", "Sex (%)", "Age and sex category (%)", "Species (%)", 
                                                         "Flexneri serotypes(%)", "Origin of sample (%)","Travel (%)", 
                                                         "Travel continents (%)", "Year of samples (%)"), stringsAsFactors = FALSE) 
all.res <- names_corr %>% right_join(all.res) 
#all.res$var <- NULL #commente 20190123

#chi2 : association avec l'origine du prélèvement
pvalues <- lapply(c("age_cut", "sexe", "age_sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"), function(var){
  subd <- d[!is.na(d[ ,var]), ]
  #pval <- chisq.test(subd[ ,var], subd$metrop, correct = FALSE)$p.value  
  #pval <- fisher.test(subd[ ,var], subd$metrop)$p.value  
  #return(data.frame(variable = var, pvalue = round(pval,3), stringsAsFactors = FALSE))
  
  #Je choisis le test en fonction des valeurs observées
  chisq <- chisq.test(table(subd[ ,var], subd$metrop), correct = FALSE)
  chisqE <- data.frame(chisq$expected)
  #chisqE[1,1] <- 2
  if(all(chisqE>5)) {
    pval <- chisq$p.value
    test <- "chisq"
  } else {
    if(all(chisqE>3)){
      pval <- chisq.test(table(subd[ ,var], subd$metrop), correct = TRUE)$p.value  
      test <- "chisq corr"
    } else {
      #Si l'espace requis est trop grand, pvalue donne NA
      pval <-  tryCatch(fisher.test(subd[ ,var], subd$metrop, workspace = 1e8)$p.value, error = function(e){NA}) 
      test <- "fisher"
    }
  }
  if (is.numeric(pval)) pval <- round(pval, 3)
 
  return(data.frame(variable = var, pvalue = pval, test, stringsAsFactors = FALSE))
}) %>% bind_rows()

# 2019 01 12
all.res <- left_join(all.res,pvalues, by = c("var"="variable"))
all.res$var <- NULL #20190123
# fin 2019 01 12

write.table(all.res, file = "clipboard", sep = "\t", row.names = FALSE, na = "")


# 2019 01 10
#Calcul de l'age median 
#domtom
myd<-d %>% filter(metrop == FALSE) 
quantile(myd$age, c(0.25,0.5,0.75), na.rm = TRUE)
#metrop 
myd<-d %>% filter(metrop == TRUE) 
quantile(myd$age, c(0.25,0.5,0.75), na.rm = TRUE)
# fin 2019 01 10



#---------------------------------
# Tableau 1 bis : description selon sensible et résistant
#---------------------------------
#20180112
#Description des souches C3G resistantes

#------------------
# création des bdd

#metrop :
d_m <- d %>% filter(metrop == TRUE)

#Je crée les bdd
CIPR.df <- d %>% filter(metrop == TRUE & CIP_SIR == 1)
CIPS.df<-d %>% filter(metrop == TRUE & CIP_SIR == 0)
C3GR.df <- d %>% filter(metrop == TRUE & C3G_SIR == 1)
C3GS.df<-d %>% filter(metrop == TRUE & C3G_SIR == 0)
AZMR.df <- d %>% filter(metrop == TRUE & AZM_SIR == 1)
AZMS.df<-d %>% filter(metrop == TRUE & AZM_SIR == 0)

# # domtom
# d_m <- d %>% filter(metrop == FALSE)
# 
# #Je crée les bdd
# CIPR.df <- d %>% filter(metrop == FALSE & CIP_SIR == 1)
# CIPS.df<-d %>% filter(metrop == FALSE & CIP_SIR == 0)
# C3GR.df <- d %>% filter(metrop == FALSE & C3G_SIR == 1)
# C3GS.df<-d %>% filter(metrop == FALSE & C3G_SIR == 0)
# AZMR.df <- d %>% filter(metrop == FALSE & AZM_SIR == 1)
# AZMS.df<-d %>% filter(metrop == FALSE & AZM_SIR == 0)

#------------------
#description

#Je sors les descriptions des variables pour chaque sous groupe : CIP_R, CIP_S, AZM_R, AZM_S, C3G_R, C3G_S

# 1- Calcul des pvalues comparant R et S pour CIP, AZM et C3G
pvalues_CIP <- lapply(c("age_cut", "sexe", "age_sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"), function(var){
  print(var)
  var_ATB = "CIP_SIR"
  subd <- d_m[!is.na(d_m [ ,var]) & !is.na(d_m[ ,var_ATB]), ]
 subd[ ,var_ATB] <- factor(subd[ ,var_ATB], c("0","1"))
  
  #Je choisis le test en fonction des valeurs observées
  chisq <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = FALSE)
  chisqE <- data.frame(chisq$expected)
  #chisqE[1,1] <- 2
  if(all(chisqE>5)) {
    pval <- chisq$p.value
    test <- "chisq"
  } else {
    if(all(chisqE>3)){
      pval <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = TRUE)$p.value  
      test <- "chisq corr"
    } else {
      #Si l'espace requis est trop grand, pvalue donne NA
      pval <-  tryCatch(fisher.test(subd[ ,var], subd[ ,var_ATB], workspace = 1e8)$p.value, error = function(e){NA}) 
      test <- "fisher"
    }
  }
  if (is.numeric(pval)) pval <- round(pval, 3)
  return(data.frame(var = var, pvalueCIP = pval, testCIP = test, stringsAsFactors = FALSE))
}) %>% bind_rows()

pvalues_AZM <- lapply(c("age_cut", "sexe", "age_sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"), function(var){
  print(var)
  var_ATB = "AZM_SIR"
  subd <- d_m[!is.na(d_m [ ,var]) & !is.na(d_m[ ,var_ATB]), ]
  subd[ ,var_ATB] <- factor(subd[ ,var_ATB], c("0","1"))
  
  #Je choisis le test en fonction des valeurs observées
  chisq <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = FALSE)
  chisqE <- data.frame(chisq$expected)
  #chisqE[1,1] <- 2
  if(all(chisqE>5)) {
    pval <- chisq$p.value
    test <- "chisq"
  } else {
    if(all(chisqE>3)){
      pval <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = TRUE)$p.value  
      test <- "chisq corr"
    } else {
      #Si l'espace requis est trop grand, pvalue donne NA
      pval <-  tryCatch(fisher.test(subd[ ,var], subd[ ,var_ATB], workspace = 1e8)$p.value, error = function(e){NA}) 
      test <- "fisher"
    }
  }
  if (is.numeric(pval)) pval <- round(pval, 3)
  
  return(data.frame(var = var, pvalueAZM = pval, testAZM = test, stringsAsFactors = FALSE))
}) %>% bind_rows()

pvalues_C3G <- lapply(c("age_cut", "sexe", "age_sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"), function(var){
  print(var)
  var_ATB = "C3G_SIR"
  subd <- d_m[!is.na(d_m [ ,var]) & !is.na(d_m[ ,var_ATB]), ]
  subd[ ,var_ATB] <- factor(subd[ ,var_ATB], c("0","1"))
  
  #Je choisis le test en fonction des valeurs observées
  chisq <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = FALSE)
  chisqE <- data.frame(chisq$expected)
  #chisqE[1,1] <- 2
  if(all(chisqE>5)) {
    pval <- chisq$p.value
    test <- "chisq"
  } else {
    if(all(chisqE>3)){
      pval <- chisq.test(table(subd[ ,var], subd[ ,var_ATB]), correct = TRUE)$p.value  
      test <- "chisq corr"
    } else {
      #Si l'espace requis est trop grand, pvalue donne NA
      pval <-  tryCatch(fisher.test(subd[ ,var], subd[ ,var_ATB], workspace = 1e8)$p.value, error = function(e){NA}) 
      test <- "fisher"
    }
  }
  if (is.numeric(pval)) pval <- round(pval, 3)
  return(data.frame(var = var, pvalueC3G = pval, testC3G = test, stringsAsFactors = FALSE))
}) %>% bind_rows()


#2- effectif et pourcentage
CIPR <- describe_shigelle(CIPR.df, "CIP_R")
CIPS <- describe_shigelle(CIPS.df, "CIP_S")
C3GR <- describe_shigelle(C3GR.df, "3GC_R")
C3GS <- describe_shigelle(C3GS.df, "3GC_S")
AZMR <- describe_shigelle(AZMR.df, "AZM_R")
AZMS <- describe_shigelle(AZMS.df, "AZM_S")

#-------------------------
# Mise en forme

#Je crée un template avec les variables et les levels 
all_tab <- list(CIPR, CIPS, C3GR, C3GS, AZMR, AZMS)
nrowtab <- lapply(all_tab, nrow) %>% unlist 
nrowtab_max <- which(nrowtab == max(nrowtab))[1]
template_col <- all_tab[[nrowtab_max]][c("var", "level")]

all_tab <- list(CIPR, CIPS, pvalues_CIP, C3GR, C3GS,pvalues_C3G, AZMR, AZMS, pvalues_AZM)

#Je colle au fur et à mesure
all_res <- template_col
for (tab in all_tab){
  all_res <- left_join(all_res, tab)  
}

#tableau avec les variables et s'il faut compter ou non les missing 
df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "age_sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                     Nado = c("a", "a", "a", "a","a", "no", "a", "a", "no", "a"), stringsAsFactors = FALSE)

names_corr <- data.frame(var = df_var$variables, all = c("N", "Age (%)", "Sex (%)", "Age and sex category (%)", "Species (%)", 
                                                         "Flexneri serotypes(%)", "Origin of sample (%)", "Travel (%)", 
                                                         "Travel continents (%)", "Year of samples (%)"), stringsAsFactors = FALSE) 
all_res <- all_res %>% group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>% 
  mutate(var = ifelse(nr !=1, NA, var), nr = NULL) 
all_res <- names_corr %>% right_join(all_res) %>% mutate(var = NULL)

#ajout 20190228
all_res$pvalueCIP <- ifelse(is.na(all_res$all), NA, all_res$pvalueCIP)
all_res$pvalueC3G <- ifelse(is.na(all_res$all), NA, all_res$pvalueC3G)
all_res$pvalueAZM <- ifelse(is.na(all_res$all), NA, all_res$pvalueAZM)
# all_res$testCIP <- ifelse(is.na(all_res$all), NA, all_res$testCIP)
# all_res$testC3G <- ifelse(is.na(all_res$all), NA, all_res$testC3G)
# all_res$testAZM <- ifelse(is.na(all_res$all), NA, all_res$testAZM)
all_res$testCIP <- NULL
all_res$testC3G <- NULL
all_res$testAZM <- NULL

write.table(all_res, file = "clipboard", sep = "\t", row.names = FALSE, na = "")
# #chi2 : association avec la C3G resistance
# pvalues <- lapply(c("age_cut", "sexe", "espece3", "serotype_flexneri", "voyage_01", "continent", "year_real_cut"), function(var){
#   subd <- d[!is.na(d[ ,var]) & d$metrop == TRUE, ]
#   #pval <- chisq.test(subd[ ,var], subd$C3G_SIR, correct = FALSE)$p.value  
#   pval <- fisher.test(subd[ ,var], subd$C3G_SIR)$p.value  
#   return(data.frame(variable = var, pvalue = round(pval,3), stringsAsFactors = FALSE))
# }) %>% bind_rows()
# 
# all.res <- left_join(all.res_C3G, pvalues, by = c("var"="variable"))
# 
# write.table(all.res, file = "clipboard", sep = "\t", row.names = FALSE, na = "")


#----------------------------------------------------------------------
# Tableau 2 : reg logistique AZM 2014-2016 sonnei et flexneri metropole
#----------------------------------------------------------------------

#tous les serotypes sont NA à partir de 2014 pour flexneri
table(dAZM$serotype, useNA = "a") #que serotype NA
table(dl$serotype_flexneri)

#nombre de resistants
#table(dAZM$status) #202 evenement donc max k = 20 (10evt*20 = 200)
table(dAZM_mod$status)#20180406 : 173 evenement donc max k = 18 (10evt*18 = 180)
table(dAZM_mod$status)#20180122 : 202 evenement donc max k = 20 (10evt*20 = 200) #c'est gracea la categorie missing de voyage

#dim(dAZM)
dim(dAZM_mod)

# 164 patients avaient une valeur manquante (plus maintenant qu'il y a categorie missing pr voyage)
apply(is.na(dAZM[ ,varAZM]), 1, sum) %>% table
#dont 1 patient manquant pour l'age et 163 patients manquants pour le voyage(plus maintenant qu'il y a categorie missing pr voyage) 
apply(is.na(dAZM[ ,varAZM]), 2, sum) 


#-----------------
#prevalence
#.l <- lapply(varAZM, function(myvar) get_prev(dAZM, myvar)) 
#.l <- lapply(varAZM, function(myvar) get_prev(dAZM_mod, myvar)) #20180406
#.l <- lapply(c(varAZM, "voyage_01_NA"), function(myvar) get_prev(dAZM_mod, myvar)) #20190115
.l <- lapply(c(varAZM), function(myvar) get_prev(dAZM_mod, myvar)) #20190122
prevAZM<-do.call(rbind, .l)     
prevAZM <- prevAZM %>% mutate(perc_resistance = paste0(perc_resistance, " (", IC_l, " - ", IC_u, ")"),
                              IC_l = NULL, IC_u = NULL)

#-----------------
# survie univariee
myunivar <- lapply(varAZM, function(myvar)get_univar2(my_data = dAZM_mod, my_outcome = "status", my_var = myvar))
myunivar <- bind_rows(myunivar)
myunivar <- myunivar %>% mutate(OR = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                                IC_l = NULL, IC_u = NULL)

#export de prevalence et de survie univariee
#prevuniv <- left_join(prevAZM, myunivar, by = c("variables", "levels")) #commente 2019 01 15
prevuniv <- full_join(prevAZM, myunivar, by = c("variables", "levels"))#2019 01 15
prevuniv$OR <- ifelse(is.na(prevuniv$OR), 1.00, prevuniv$OR)
prevuniv %>% select(-outcome) %>% write.table(file = "clipboard", sep = "\t", row.names = F)


#verif univarié pour age
mod <- glm(status~age_cut, data = dAZM_mod, family = "binomial")
summary(mod)
exp(confint(mod))
qnorm(0.975)
exp(1.5568 - qnorm(0.975)*0.3193)
#-----------------
#survie multivariée
#les conditions pour le modele sont repectees, j'ai k = 6 pour 203 evt


#automatique
myvar_m <- myunivar %>% 
  filter(pvalue <0.25) %>% 
  pull(variables) %>% unique %>% as.character()
myvarm_vec <- paste(myvar_m, collapse = "+")
mm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dAZM_mod)
modm <- step(mm, direction = "backward")

#manuel

#etape 1
modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dAZM_mod)
OR <- exp(coef(modm))
CI <- exp(confint(modm))
ORCI <- round(cbind(OR, CI), 2); ORCI

all_pval <- unlist(lapply(myvar_m, function(var)get_pval_multivar(var_to_test = var, all_var = myvar_m, my_outcome = "status", my_data = dAZM_mod)))
data.frame(myvar_m, pval = all_pval) 

#etape 2: retrait de espece3
all_univar <- data.frame(myvar_m, pval = all_pval) %>% 
  mutate(max_p = max(pval), 
         to_del = ifelse(pval == max(pval) & pval >= 0.05, TRUE, FALSE)) %>% 
  filter(to_del == FALSE) %>% 
  pull(myvar_m) %>% as.character; all_univar
myvarm_vec <- paste(all_univar, collapse = "+")
modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dAZM_mod)
OR <- exp(coef(modm))
CI <- exp(confint(modm))
ORCI <- round(cbind(OR, CI), 2); ORCI
#les OR n'ont pas beaucoup changé
all_pval <- unlist(lapply(all_univar, function(var) get_pval_multivar(var_to_test = var, all_var = all_univar, my_outcome = "status", my_data = dAZM_mod)))

#30/03 mise en forme 
#aucune autre variable a retirer
all_pval<-data.frame(variables = all_univar, pval = round(all_pval, 3)) 
ORCI <- data.frame(ORCI)
names(ORCI) <- c("OR", "IC_l", "IC_u")
ORCI$variables <- str_extract(row.names(ORCI), paste(prevuniv$variables, collapse = "|"))
ORCI$levels <- mapply(my_sub, row.names(ORCI), ORCI$var)
ORCI <- ORCI %>% mutate(OR_mult = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                        IC_l = NULL, IC_u = NULL, OR = NULL) %>% filter(!is.na(variables)) %>% 
  right_join(all_pval)
#export des resultats complets
prev <- prevAZM %>% left_join(., myunivar, by = c("variables", "levels")) %>% 
  left_join (., ORCI, by = c("variables", "levels")) %>% 
  group_by(variables) %>% 
  mutate(nr = row_number(),
         pvalue = as.numeric(ifelse(nr == 1, lead(pvalue), NA)),
         pval = as.numeric(ifelse(nr == 1, lead(pval), NA)),
         OR_mult = ifelse(nr == 1 & pvalue<0.05, 1, OR_mult),
         OR = ifelse(nr == 1, 1, OR),
         nr = NULL)
prev$pvalue <- ifelse(prev$pvalue == 0, "<0.001", prev$pvalue)
prev$pval <- ifelse(prev$pval == 0, "<0.001", prev$pval)
prev %>% write.table(file = "clipboard", row.names = FALSE, sep = "\t", na = "")
#----------------------------------------------------------------------
# Tableau 3 : reg logistique C3G 2005-2016 sonnei metropole
#----------------------------------------------------------------------

#nombre de resistants (=nb d'evt)
table(dC3G$status)
#92 evenement donc max k = 9 (10evt*9 = 90)
table(dC3G_mod$status)
#79 evenement donc max k = 8 (10evt*8 = 80)

# distribution des valeurs manquantes
apply(is.na(dC3G[ ,varC3G]), 1, sum) %>% table
apply(is.na(dC3G[ ,varC3G]), 2, sum) 

#-----------------
#prevalence
#.l <- lapply(varC3G, function(myvar) get_prev(dC3G, myvar))
.l <- lapply(varC3G, function(myvar) get_prev(dC3G_mod, myvar))#20180406
prevC3G<-do.call(rbind, .l)     
prevC3G <- prevC3G %>% mutate(perc_resistance = paste0(perc_resistance, " (", IC_l, " - ", IC_u, ")"),
                              IC_l = NULL, IC_u = NULL)

#-----------------
# survie univariee
myunivar <- lapply(varC3G, function(myvar)get_univar2(my_data = dC3G_mod, my_outcome = "status", my_var = myvar))
myunivar <- bind_rows(myunivar)
myunivar <- myunivar %>% mutate(OR = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                                IC_l = NULL, IC_u = NULL)

#export de prevalence et de survie univariee
prevuniv <- left_join(prevC3G, myunivar, by = c("variables", "levels"))
prevuniv$OR <- ifelse(is.na(prevuniv$OR), 1.00, prevuniv$OR)
prevuniv %>% select(-outcome) %>% write.table(file = "clipboard", sep = "\t", row.names = F)

#-----------------
#survie multivariée

#automatique
myvar_m <- myunivar %>% 
  filter(pvalue <0.25) %>% 
  pull(variables) %>% unique %>% as.character()
myvarm_vec <- paste(myvar_m, collapse = "+")
mm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
modm <- step(mm, direction = "backward")

#etape 1: modele plein
modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
OR <- exp(coef(modm))
CI <- exp(confint(modm))
ORCI <- round(cbind(OR, CI), 2); ORCI

all_pval <- unlist(lapply(myvar_m, function(var)get_pval_multivar(var_to_test = var, all_var = myvar_m, my_outcome = "status", my_data = dC3G_mod)))
data.frame(myvar_m, pval = all_pval)

#etape 2: retrait de voyage
all_univar <- data.frame(myvar_m, pval = all_pval) %>% 
  mutate(max_p = max(pval), 
         to_del = ifelse(pval == max(pval) & pval >= 0.05, TRUE, FALSE)) %>% 
  filter(to_del == FALSE) %>% 
  pull(myvar_m) %>% as.character; all_univar
myvarm_vec <- paste(all_univar, collapse = "+")
modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
OR <- exp(coef(modm))
CI <- exp(confint(modm))
ORCI <- round(cbind(OR, CI), 2); ORCI
#les OR n'ont pas beaucoup changé
all_pval <- unlist(lapply(all_univar, function(var) get_pval_multivar(var_to_test = var, all_var = all_univar, my_outcome = "status", my_data = dC3G_mod)))

# #les conditions pour le modele sont repectees, j'ai k = 7 pour 92 evt
# #etape 1: modele plein
# all_univar <- unique(myunivar$variables)
# myvarm_vec <- paste(all_univar, collapse = "+")
# modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
# OR <- exp(coef(modm))
# CI <- exp(confint(modm))
# ORCI <- round(cbind(OR, CI), 2); ORCI
# 
# all_pval <- unlist(lapply(all_univar, function(var)get_pval_multivar(var_to_test = var, all_var = all_univar, my_outcome = "status", my_data = dC3G_mod)))
# data.frame(all_univar, pval = all_pval) 
# 
# #etape 2 retrait de age
# all_univar <- data.frame(all_univar, pval = all_pval) %>% arrange(pval)
# all_univar <- c("year_real_cut", "sexe", "voyage_01")
# myvarm_vec <- paste(all_univar, collapse = "+")
# modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
# OR <- exp(coef(modm))
# CI <- exp(confint(modm))
# ORCI <- round(cbind(OR, CI), 2); ORCI
# #les OR n'ont pas beaucoup changé
# all_pval <- unlist(lapply(all_univar, function(var) get_pval_multivar(var_to_test = var, all_var = all_univar, my_outcome = "status", my_data = dC3G_mod)))
# data.frame(all_univar, pval = all_pval) 
# 
# #etape 3 retrait de voyage_01
# all_univar <- data.frame(all_univar, pval = all_pval) %>% arrange(pval)
# all_univar <- c("year_real_cut", "sexe")
# myvarm_vec <- paste(all_univar, collapse = "+")
# modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
# OR <- exp(coef(modm))
# CI <- exp(confint(modm))
# ORCI <- round(cbind(OR, CI), 2); ORCI
# #les OR n'ont pas beaucoup changé
# all_pval <- unlist(lapply(all_univar, function(var) get_pval_multivar(var_to_test = var, all_var = all_univar, my_outcome = "status", my_data = dC3G_mod)))
# data.frame(all_univar, pval = all_pval) 
# 
# #etape 4 retrait de sexe
# all_univar <- data.frame(all_univar, pval = all_pval) %>% arrange(pval)
# all_univar <- c("year_real_cut")
# myvarm_vec <- paste(all_univar, collapse = "+")
# modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dC3G_mod)
# OR <- exp(coef(modm))
# CI <- exp(confint(modm))
# ORCI <- round(cbind(OR, CI), 2); ORCI
# #les OR n'ont pas beaucoup changé
# get_univar2 (dC3G_mod, "status", all_univar)




#30/03 mise en forme 
#aucune autre variable a retirer
all_pval<-data.frame(variables = all_univar, pval = round(all_pval, 3)) 
ORCI <- data.frame(ORCI)
names(ORCI) <- c("OR", "IC_l", "IC_u")
ORCI$variables <- str_extract(row.names(ORCI), paste(prevuniv$variables, collapse = "|"))
ORCI$levels <- mapply(my_sub, row.names(ORCI), ORCI$var)
ORCI <- ORCI %>% mutate(OR_mult = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                        IC_l = NULL, IC_u = NULL, OR = NULL) %>% filter(!is.na(variables)) %>% 
  right_join(all_pval)
#export des resultats complets
prev <- prevC3G %>% left_join(., myunivar, by = c("variables", "levels")) %>% 
  left_join (., ORCI, by = c("variables", "levels")) %>% 
  group_by(variables) %>% 
  mutate(nr = row_number(),
         pvalue = as.numeric(ifelse(nr == 1, lead(pvalue), NA)),
         pval = as.numeric(ifelse(nr == 1, lead(pval), NA)),
         OR_mult = ifelse(nr == 1 & pvalue<0.05, 1, OR_mult),
         OR = ifelse(nr == 1, 1, OR),
         nr = NULL)
prev$pvalue <- ifelse(prev$pvalue == 0, "<0.001", prev$pvalue)
prev$pval <- ifelse(prev$pval == 0, "<0.001", prev$pval)
prev %>% write.table(file = "clipboard", row.names = FALSE, sep = "\t", na = "")# ORCI %>% write.table(file = "clipboard", sep = "\t")
# data.frame(all_univar, pval = all_pval) 

#nombre de missing
table(is.na(dC3G[,"age_cut"]))
table(is.na(dC3G[,"sexe"]))
table(is.na(dC3G[,"voyage_01_m"]))
table(is.na(dC3G[,"year_real_cut"]))



#-----------------------
# autres

#nombre de pansensible
subd <- d %>% filter(metrop == TRUE) %>% select(AKN_SIR, AMXAMP_SIR, C3G_SIR, CAZ1030_SIR, CHL_SIR, CIP_SIR, CRO_SIR, GEN_SIR,
                                                NAL_SIR, SMX_SIR, STR_SIR, SXT_SIR, TET_SIR, TMP_SIR)
table(apply(subd, 1, sum, na.rm = TRUE), useNA = "a") #264 pansensible, 7939 avec au moins une resistance

subd <- d %>% filter(metrop == FALSE) %>% select(AKN_SIR, AMXAMP_SIR, C3G_SIR, CAZ1030_SIR, CHL_SIR, CIP_SIR, CRO_SIR, GEN_SIR,
                                                 NAL_SIR, SMX_SIR, STR_SIR, SXT_SIR, TET_SIR, TMP_SIR)
table(apply(subd, 1, sum, na.rm = TRUE), useNA = "a") #264 pansensible, 7939 avec au moins une resistance


#-------------------------
#-------------------------
#20180108

#Tableau 2 pour CIP resistance

#nombre de resistants (=nb d'evt)
table(dCIP_mod$status)
#667 evenement donc max k = 66 (10evt*66 = 660)

# distribution des valeurs manquantes
apply(is.na(dCIP[ ,varCIP]), 1, sum) %>% table
apply(is.na(dCIP[ ,varCIP]), 2, sum) 

#-----------------
#prevalence
#.l <- lapply(varCIP, function(myvar) get_prev(dCIP, myvar))
.l <- lapply(varCIP, function(myvar) get_prev(dCIP_mod, myvar))#20180406
prevCIP<-do.call(rbind, .l)     
prevCIP <- prevCIP %>% mutate(perc_resistance = paste0(perc_resistance, " (", IC_l, " - ", IC_u, ")"),
                              IC_l = NULL, IC_u = NULL)
#-----------------
# survie univariee
myunivar <- lapply(varCIP, function(myvar)get_univar2(my_data = dCIP_mod, my_outcome = "status", my_var = myvar))
myunivar <- bind_rows(myunivar)
myunivar <- myunivar %>% mutate(OR = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                                IC_l = NULL, IC_u = NULL)

#export de prevalence et de survie univariee
prevuniv <- left_join(prevCIP, myunivar, by = c("variables", "levels"))
prevuniv$OR <- ifelse(is.na(prevuniv$OR), 1.00, prevuniv$OR)
prevuniv %>% select(-outcome) %>% write.table(file = "clipboard", sep = "\t", row.names = F)

#-----------------
#survie multivariée

#automatique
myvar_m <- myunivar %>% 
  filter(pvalue <0.25) %>% 
  pull(variables) %>% unique %>% as.character()
myvarm_vec <- paste(myvar_m, collapse = "+")
mm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dCIP_mod)
modm <- step(mm, direction = "backward")

#etape 1: modele plein
modm <- glm(formula(paste0("status ~ ", myvarm_vec)), family = "binomial", data = dCIP_mod)
OR <- exp(coef(modm))
CI <- exp(confint(modm))
ORCI <- round(cbind(OR, CI), 2); ORCI

all_pval <- unlist(lapply(myvar_m, function(var)get_pval_multivar(var_to_test = var, all_var = myvar_m, my_outcome = "status", my_data = dCIP_mod)))
data.frame(myvar_m, pval = all_pval)

#aucune variable à retirer
all_univar <- data.frame(myvar_m, pval = all_pval) %>% 
  mutate(max_p = max(pval), 
         to_del = ifelse(pval == max(pval) & pval >= 0.05, TRUE, FALSE)) %>% 
  filter(to_del == FALSE) %>% 
  pull(myvar_m) %>% as.character; all_univar

#30/03 mise en forme 
#aucune autre variable a retirer
all_pval<-data.frame(variables = all_univar, pval = round(all_pval, 3)) 
ORCI <- data.frame(ORCI)
names(ORCI) <- c("OR", "IC_l", "IC_u")
ORCI$variables <- str_extract(row.names(ORCI), paste(prevuniv$variables, collapse = "|"))
ORCI$levels <- mapply(my_sub, row.names(ORCI), ORCI$var)
ORCI <- ORCI %>% mutate(OR_mult = paste0(OR, " (", IC_l, " - ", IC_u, ")"),
                        IC_l = NULL, IC_u = NULL, OR = NULL) %>% filter(!is.na(variables)) %>% 
  right_join(all_pval)
#export des resultats complets
prev <- prevCIP %>% left_join(., myunivar, by = c("variables", "levels")) %>% 
  left_join (., ORCI, by = c("variables", "levels")) %>% 
  group_by(variables) %>% 
  mutate(nr = row_number(),
         pvalue = as.numeric(ifelse(nr == 1, lead(pvalue), NA)),
         pval = as.numeric(ifelse(nr == 1, lead(pval), NA)),
         OR_mult = ifelse(nr == 1 & pvalue<0.05, 1, OR_mult),
         OR = ifelse(nr == 1, 1, OR),
         nr = NULL)
prev$pvalue <- ifelse(prev$pvalue == 0, "<0.001", prev$pvalue)
prev$pval <- ifelse(prev$pval == 0, "<0.001", prev$pval)
prev %>% select(-outcome) %>% write.table(file = "clipboard", row.names = FALSE, sep = "\t", na = "")# ORCI %>% write.table(file = "clipboard", sep = "\t")
# data.frame(all_univar, pval = all_pval) 

#-------------------------
#-------------------------
#20180112
#Description des souches C3G resistantes

#domtom
my_data<-d %>% filter(metrop == FALSE & C3G_SIR == 1) 
dim(my_data)
#metrop 
my_data<-d %>% filter(metrop == TRUE & C3G_SIR == 1) 
dim(my_data)

df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                     Nado = c("a", "a", "a", "a", "no", "a", "no", "a"), stringsAsFactors = FALSE)

all.res_C3G <- 
  #calcul effectif et pourcentage pour chaque variable
  lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    var_qual(name, my_data, Nado)
  }) %>% 
  #regrouper variables
  bind_rows() %>% 
  #mise en forme
  group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>% 
  mutate(
    #variable que sur la première ligne
    #var = ifelse(nr !=1, NA, var),
    #N(%)
    freq_metrop_C3G = paste0(freq, " (", percent, ")"),
    #Level missing au lieu de NA
    level = ifelse(is.na(level), "missing", level),
    #Je ne garde que var level et freq metrop
    freq = NULL, percent = NULL
    #, nr = NULL
  )

# #chi2 : association avec la C3G resistance
# pvalues <- lapply(c("age_cut", "sexe", "espece3", "serotype_flexneri", "voyage_01", "continent", "year_real_cut"), function(var){
#   subd <- d[!is.na(d[ ,var]) & d$metrop == TRUE, ]
#   #pval <- chisq.test(subd[ ,var], subd$C3G_SIR, correct = FALSE)$p.value  
#   pval <- fisher.test(subd[ ,var], subd$C3G_SIR)$p.value  
#   return(data.frame(variable = var, pvalue = round(pval,3), stringsAsFactors = FALSE))
# }) %>% bind_rows()
# 
# all.res <- left_join(all.res_C3G, pvalues, by = c("var"="variable"))
# 
# write.table(all.res, file = "clipboard", sep = "\t", row.names = FALSE, na = "")

##################
# C3G R en 2005 et 2006 : 1 C3G resistant par an
my_data %>% filter(year %in% c(2005, 2006))
# C3G R en 2016 : 28 C3G resistants : detail des especes et serotypes
my_data %>% filter(year %in% c(2016)) %>% count(espece)

#origine des souches
my_data %>% mutate(labo = ifelse(labo == "institut", "hospitalise", labo)) %>% 
  count(labo) %>% mutate(n/sum(n))
#notion de voyage
my_data %>% mutate(voyage = 
                     ifelse(voyage %in% c("non precise", "non renseigne"), NA, voyage)) %>% 
  count(voyage) %>% mutate(n/sum(n)) %>% View()


###################################
###################################
# Descripton des souches CIP resistants

#domtom
my_data<-d %>% filter(metrop == FALSE & CIP_SIR == 1) 
dim(my_data)
#metrop 
my_data<-d %>% filter(metrop == TRUE & CIP_SIR == 1) 
dim(my_data)

df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                     Nado = c("a", "a", "a", "a", "no", "a", "no", "a"), stringsAsFactors = FALSE)

all.res_CIP <- 
  #calcul effectif et pourcentage pour chaque variable
  lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    var_qual(name, my_data, Nado)
  }) %>% 
  #regrouper variables
  bind_rows() %>% 
  #mise en forme
  group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>% 
  mutate(
    #variable que sur la première ligne
    #var = ifelse(nr !=1, NA, var),
    #N(%)
    freq_metrop_CIP = paste0(freq, " (", percent, ")"),
    #Level missing au lieu de NA
    level = ifelse(is.na(level), "missing", level),
    #Je ne garde que var level et freq metrop
    freq = NULL, percent = NULL
    #, nr = NULL
  )


###################################
###################################
# Descripton des souches AZM resistants

#domtom
my_data<-d %>% filter(metrop == FALSE & AZM_SIR == 1) 
dim(my_data)
#metrop 
my_data<-d %>% filter(metrop == TRUE & AZM_SIR == 1) 
dim(my_data)

df_var <- data.frame(variables = c("metrop", "age_cut", "sexe", "espece3", "serotype_flexneri", "labo", "voyage_01_m", "continent", "year_real_cut"),
                     Nado = c("a", "a", "a", "a", "no", "a", "no", "a"), stringsAsFactors = FALSE)

all.res_AZM <- 
  #calcul effectif et pourcentage pour chaque variable
  lapply(1:nrow(df_var), function(i){
    name <- df_var$variables[i]
    Nado <- df_var$Nado[i]
    var_qual(name, my_data, Nado)
  }) %>% 
  #regrouper variables
  bind_rows() %>% 
  #mise en forme
  group_by(var) %>% mutate(nr = row_number()) %>% ungroup %>% 
  mutate(
    #variable que sur la première ligne
    #var = ifelse(nr !=1, NA, var),
    #N(%)
    freq_metrop_AZM = paste0(freq, " (", percent, ")"),
    #Level missing au lieu de NA
    level = ifelse(is.na(level), "missing", level),
    #Je ne garde que var level et freq metrop
    freq = NULL, percent = NULL
    #, nr = NULL
  )


all.res <- full_join(all.res_AZM, all.res_C3G, by = c("var", "level")) %>% 
  left_join(all.res_CIP, by = c("var", "level")) %>%
  replace(., is.na(.), "0(0)") %>% 
  #group_by(var) %>% arrange(nr, .by_group = TRUE) %>% ungroup() %>% 
  #arrange(var, level, nr) %>% 
  mutate(var = ifelse(nr !=1, NA, var)) %>% select(-contains("nr"))
write.table(all.res, file = "clipboard", sep = "\t", row.names = FALSE, na = "")

#-------------------------------------
#-------------------------------------
# Tableau prévalence toutes années confondues, par souche


#all_ATB <- c("AKN_SIR", "GEN_SIR", "CHL_SIR", "SXT_SIR", "SMX_SIR", "TMP_SIR", "NAL_SIR", "CIP_SIR", "C3G_SIR", "AMXAMP_SIR", "AZM_SIR", "TET_SIR")

dl %>%
  group_by(metrop, espece3, ATB,status) %>% 
  summarise (M=n()) %>% 
  mutate(N = sum(M), prev = M/N) %>% 
  filter(status == 1)


#2-écrire la fonction sans faire appel à data ni aux arguments
#my_fun(arg1, arg2)
get_prev2 <- function(M,n){
  btest <- binom.test(x = M, n = n)
  p0 <- as.numeric(btest$estimate)
  IC <- as.numeric(btest$conf.int)
  IC <- round(c(p0, IC)*100, 2) # en %
  return(data.frame(n=n, M=M, p0=IC[1], IC_l=IC[2], IC_u=IC[3]))#attention la sortie doit être un data.frame (même si une seule valeur)
}

#3-récupérer les arguments nécessaires
tmp2 <- dl %>% group_by(metrop, espece3, ATB) %>% summarise(n = n(), M = sum(status)) %>% ungroup()
#4-group_by  %>%  do(my_fun(.$arg1, .$arg2))
prev_atb_species <- tmp2%>% group_by(metrop, espece3, ATB) %>% do(get_prev2(.$M, .$n))

prev_atb_species <- prev_atb_species %>% 
  mutate(N = paste0(n, ", ", p0, " (", IC_l, ";", IC_u, ")")) %>% 
  select(metrop, ATB, espece3, N) %>% 
  spread(key = espece3, value = N)

prev_atb_species$ATB <- factor(prev_atb_species$ATB, all_ATB)
prev_atb_species <- prev_atb_species %>% arrange(metrop, ATB) %>% filter(!is.na(ATB)) 

write.table(prev_atb_species, file = "clipboard", sep = "\t", row.names = FALSE, na = "")


#Nombre de souches
d %>% group_by(metrop, espece3) %>% summarise(n = n())


#---------------------------------------------------
# Tableaux complementaires demande Sophie 20190503
#----------------------------------------------------
# Resistance par espece pour les categories age sexe

N_age_sexe<-
bind_rows(
#CIP Sonnei
d %>% filter(metrop == TRUE & espece3 == "sonnei") %>% 
  group_by(age_sexe) %>% count(CIP_SIR) %>% 
  group_by(CIP_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = CIP_SIR) %>% 
  mutate(espece = "sonnei", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("CIP_", ATB)),
#CIP flexneri
d %>% filter(metrop == TRUE & espece3 == "flexneri") %>% 
  group_by(age_sexe) %>% count(CIP_SIR) %>% 
  group_by(CIP_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = CIP_SIR) %>% 
  mutate(espece = "flexneri", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("CIP_", ATB))
,
#AZM sonnei
d %>% filter(metrop == TRUE & espece3 == "sonnei") %>% 
  group_by(age_sexe) %>% count(AZM_SIR) %>% 
  group_by(AZM_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = AZM_SIR) %>% 
  mutate(espece = "sonnei", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("AZM_", ATB)),
#AZM flexneri
d %>% filter(metrop == TRUE & espece3 == "flexneri") %>% 
  group_by(age_sexe) %>% count(AZM_SIR) %>% 
  group_by(AZM_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = AZM_SIR) %>% 
  mutate(espece = "flexneri", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("AZM_", ATB))
,
#C3G sonnei
d %>% filter(metrop == TRUE & espece3 == "sonnei") %>% 
  group_by(age_sexe) %>% count(C3G_SIR) %>% 
  group_by(C3G_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = C3G_SIR) %>% 
  mutate(espece = "sonnei", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("C3G_", ATB)),
#C3G flexneri
d %>% filter(metrop == TRUE & espece3 == "flexneri") %>% 
  group_by(age_sexe) %>% count(C3G_SIR) %>% 
  group_by(C3G_SIR) %>% mutate(n= paste0(n, " (",  round(n/sum(n)*100, 1), ")")) %>% 
  ungroup %>%  spread(key = age_sexe, value = n) %>% rename(ATB = C3G_SIR) %>% 
  mutate(espece = "flexneri", ATB = ifelse (ATB==0, "S","R"), ATB = paste0("C3G_", ATB))
) %>% replace(., is.na(.), "0(0)")  #%>% mutate(var = str_sub(ATB, 1, 3))

d_m <- d %>% filter(metrop == TRUE)

pvalue <- function(mon_espece, var_ATB){
  d_sub <- d_m[d_m$espece3 %in% mon_espece, ]
  subd <- d_sub[!is.na(d_sub$age_sexe) & !is.na(d_sub[ ,var_ATB]), ]
  #Je choisis le test en fonction des valeurs observées
  chisq <- chisq.test(table(subd$age_sexe, subd[ ,var_ATB]), correct = FALSE)
  chisqE <- data.frame(chisq$expected)
  #chisqE[1,1] <- 2
  if(all(chisqE>5)) {
    pval <- chisq$p.value
    test <- "chisq"
  } else {
    if(all(chisqE>3)){
      pval <- chisq.test(table(subd$age_sexe, subd[ ,var_ATB]), correct = TRUE)$p.value  
      test <- "chisq corr"
    } else {
      #Si l'espace requis est trop grand, pvalue donne NA
      pval <-  tryCatch(fisher.test(subd$age_sexe, subd[ ,var_ATB], workspace = 1e8)$p.value, error = function(e){NA}) 
      test <- "fisher"
    }
  }
  if (is.numeric(pval)) {
    pval <- if(pval > 0.001) round(pval, 3) else "<0.001"
  }
  pval <- as.character(pval)
  return(data.frame(ATB = paste0(str_sub(var_ATB,1,3), "_S"), espece=mon_espece, pvalue = pval, test = test, stringsAsFactors = FALSE))
}
pval_age_sexe <- bind_rows(
  pvalue(mon_espece = "flexneri", var_ATB = "C3G_SIR"),
  pvalue(mon_espece = "flexneri", var_ATB = "AZM_SIR"),
  pvalue(mon_espece = "flexneri", var_ATB = "CIP_SIR"),
  pvalue(mon_espece = "sonnei", var_ATB = "C3G_SIR"),
  pvalue(mon_espece = "sonnei", var_ATB = "AZM_SIR"),
  pvalue(mon_espece = "sonnei", var_ATB = "CIP_SIR")
)

all_cat <- left_join(N_age_sexe, pval_age_sexe, by = c("ATB", "espece")) %>% select(ATB, espece, everything())
all_cat %>% write.table(file = "clipboard", row.names = FALSE, sep = "\t", na = "")

  
  
