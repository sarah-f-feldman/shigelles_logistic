###################################
#   Additional stat for Results   #
###################################




#---------------------------------
# packages
#---------------------------------

source("src/00_library.R")
#---------------------------------
# objets
#---------------------------------
source("src/01_data_management.R")

all_ATB <- c("AMXAMP_SIR", "C3G_SIR", "NAL_SIR", "CIP_SIR", "STR_SIR", "AKN_SIR", "GEN_SIR", "SMX_SIR", "TMP_SIR",
             "SXT_SIR", "CHL_SIR", "TET_SIR", "AZM_SIR")

#---------------------------------
# calculs
#---------------------------------

#Nombre de souche
d %>% group_by(metrop) %>% count(n())

#Nombre de patients : pas obtenable en l'etat, il faut reprendre la base non anonymisee, creer un ID pat et merger by num CNR


#Souches entierement sensibles
tmp <- d
tmp$n_R <- rowSums(tmp[, all_ATB], na.rm=T)
tmp %>% group_by(metrop, espece3) %>% filter(n_R == 0) %>% 
  #select(num_CNR, metrop, espece3, AKN_SIR, GEN_SIR, CHL_SIR, SXT_SIR, NAL_SIR, CIP_SIR, C3G_SIR, AZM_SIR, n_R) %>% View
  summarise(n_souches_totalement_S = n()) %>% 
  right_join(., tmp %>% group_by(metrop, espece3) %>% 
               summarise(tot = n()))  %>%
  mutate(n_souches_totalement_S = ifelse(is.na(n_souches_totalement_S), 0, n_souches_totalement_S), 
         perc = paste0(round(n_souches_totalement_S/tot *100, 2), " %"))


#Resistance cipro
dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "CIP_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% 
  filter(status == 1) %>% 
  mutate(perc = sum(n)/tot) %>% 
  filter(year_real %in% c(2005, 2016))

#Resistance 3GC
dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "C3G_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% 
  filter(status == 1) %>% 
  mutate(perc = sum(n)/tot) %>% 
  filter(year_real %in% c(2006, 2016))
#Resistance AZM
dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "AZM_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% 
  filter(status == 1) %>% 
  mutate(perc = sum(n)/tot) 

#------------------------------
# description des manquants retires dans les modeles
# Rappel : pour chaque ATB on part de dCIP, et on fait un na.omit(dCIP[ ,c("status", varCIP)]) pour avoir dCIP_mod

nrow(dCIP)
nrow(dCIP)-nrow(dCIP_mod)
apply(is.na(dCIP[ ,c("status", varCIP)]), 2, sum)
dCIP %>% filter(is.na(age_cut) & !is.na(sexe)) %>% count
dCIP %>% filter(!is.na(age_cut) & is.na(sexe)) %>% count
dCIP %>% filter(is.na(age_cut) & is.na(sexe)) %>% count

nrow(dC3G)
nrow(dC3G)-nrow(dC3G_mod)
apply(is.na(dC3G[ ,c("status", varC3G)]), 2, sum)
dC3G %>% filter(is.na(age_cut) & !is.na(sexe)) %>% count
dC3G %>% filter(!is.na(age_cut) & is.na(sexe)) %>% count
dC3G %>% filter(is.na(age_cut) & is.na(sexe)) %>% count

nrow(d %>% filter(year_real %in% 2014:2016))
nrow(dAZM)
nrow(dAZM)-nrow(dAZM_mod)
apply(is.na(dAZM[ ,c("status", varAZM)]), 2, sum)
dAZM %>% filter(is.na(age_cut) & !is.na(sexe)) %>% count
dAZM %>% filter(!is.na(age_cut) & is.na(sexe)) %>% count
dAZM %>% filter(is.na(age_cut) & is.na(sexe)) %>% count


#----------------------------------
# 20190618 Chi2 de tendance pour CIP

aa <- dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "CIP_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% spread(key = status, value = n)

prop.trend.test(aa$'1', aa$tot, score = 1:nrow(aa))
#prop.trend.test(aa$'1', aa$tot, score = 0:(nrow(aa)-1))#idem

# 20190618 Chi2 de tendance pour AZM

aa <- dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "AZM_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% spread(key = status, value = n)

prop.trend.test(aa$'1', aa$tot, score = 1:nrow(aa))



# 20190618 Chi2 de tendance pour 3GC

aa <- dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "C3G_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% spread(key = status, value = n)

prop.trend.test(aa$'1', aa$tot, score = 1:nrow(aa))


#----------------------------------
#2019 11 08 origine du prelevement

prel.vec <- d %>% filter(metrop == TRUE) %>% pull(origine_prel) 
prel.vec <-  tolower(prel.vec)
prel.vec <- as.factor(prel.vec)

levels(prel.vec) <- c("anal", "Anal/rectal", "Autre", "biopsie", "Biopsie", "Biopsie colique", 
                      "Biopsie rectale", "Crachat", "Inconnue", "Liquide cephalo-rachidien", 
                      "Liquide gastrique", "peritoneal", "Peritoneal", "prelevement humain", 
                      "Prelevement vaginal", "prothese", "pus", "sang", "Sang", "selles", 
                      "Selles", "urethral", "urines", "Urines", "vaginal")      
 
#selles       
prop.table(table(prel.vec, useNA = "a"))

#autres origines gastriques
mysum <- sum(table(prel.vec, useNA = "a")[names(table(prel.vec, useNA = "a")) %in% c("anal", "anal/rectal", "biopsie colique",
                                              "biopsie rectale", "liquide gastrique")]); mysum
mysum/sum(table(prel.vec, useNA = "a"))

#autres origines 
mysum <- sum(table(prel.vec, useNA = "a")[names(table(prel.vec, useNA = "a")) %in% c("biopsie", "crachat",
                                                           "liquide cephalo-rachidien", "peritoneal",
                                                            "prelevement vaginal", 
                                                           "prothese", "pus", "sang", "urethral", "urines", "vaginal")]); mysum
mysum/sum(table(prel.vec, useNA = "a"))

#origine inderterminee
mysum <- sum(table(prel.vec, useNA = "a")[names(table(prel.vec, useNA = "a")) %in% c("prelevement humain","autre", NA)]); mysum
mysum/sum(table(prel.vec, useNA = "a"))       

#----------------------------------
#2019 11 08 serotype
d %>% filter(metrop == TRUE) %>% group_by (espece3, espece) %>%
  summarize(n()) %>% data.frame() %>% 
  write.table(sep ="\t", file = "clipboard", row.names = FALSE)
d %>% filter(metrop == TRUE &  CIP_SIR==1 & year_real == 2016) %>%
  filter(str_detect(espece, "sonnei g")) %>% dim
vec <- d %>% filter(metrop == TRUE & !is.na(serotype)& !is.na(espece3)) %>% mutate(ser = paste(espece3, serotype), 
                                                                  ser = ifelse(str_detect(ser,"6 boyd.. 88"), "flexneri 6 boyd 88", ser )
                                                                  ) %>% pull(ser) %>% tolower %>% table(useNA = "no") #%>% prop.table %>% round(4)
#non il y a des serotype na alors que espece non NA...
vec <- d %>% filter(metrop == TRUE & !is.na(espece)) %>% mutate(espece = tolower(espece),
                                                                ser = str_replace(espece, "shigella",""),
                                                                ser = ifelse(str_detect(ser,"6 boyd.. 88"), "flexneri 6 boyd 88", ser ),
                                                                ser = ifelse (str_detect(espece, "sonnei g"), "sonnei g", ser)) %>% pull(ser) %>% table(useNA = "no") #%>% prop.table %>% round(4)
vec[order(as.numeric(vec), decreasing = TRUE)] %>% data.frame()

#Travel CIP 7 dysenteriae 1

d %>% filter(metrop == TRUE & CIP_SIR =="1") %>% View
d %>% filter(espece3 == "dysenteriae" & serotype == "1") %>% View()
#notion de voyage
#172 souches CIP-R d'inde
d %>% filter(metrop == TRUE & voyage=="inde" & CIP_SIR==1) %>% View

d %>% filter(metrop == TRUE & year_real == 2016 & CIP_SIR==1) %>% View
#voyage 2016 CIP-R
d %>% filter(metrop == TRUE & year_real == 2016 & CIP_SIR==1) %>% pull(voyage) %>% table(useNA = "a") %>% data.frame
114-50-20

#3GC : species
dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "C3G_SIR" & status == "1") %>% group_by(espece3) %>% count
#voyage C3G
d %>% filter(metrop == TRUE & C3G_SIR==1) %>% pull(continent) %>% table
d %>% filter(metrop == TRUE & C3G_SIR==1) %>% pull(voyage) %>% table(useNA = "a") %>% data.frame
105-54-1-11

#AZM
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% dim
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% filter(str_detect(espece, "dysenteriae")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% filter(str_detect(espece, "boydii")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% filter(str_detect(espece, "flexneri")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% filter(str_detect(espece, "sonnei")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & year_real == 2016 & AZM_SIR==1) %>% filter(str_detect(espece, "flexneri")) %>%  pull(espece) %>% table %>% data.frame()
d %>% filter(metrop == TRUE & AZM_SIR==1) %>% pull(voyage) %>% table(useNA = "a") %>% data.frame

#GEN
d %>% filter(metrop == TRUE & GEN_SIR==1) %>% dim
dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "GEN_SIR") %>% 
  select(year_real, espece3, status) %>% 
  group_by(year_real) %>% 
  count(status) %>% 
  mutate(tot = sum(n)) %>% spread(key = status, value = n)
d %>% filter(metrop == TRUE &  GEN_SIR==1) %>% filter(str_detect(espece, "dysenteriae")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & GEN_SIR==1) %>% filter(str_detect(espece, "boydii")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & GEN_SIR==1) %>% filter(str_detect(espece, "flexneri")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & GEN_SIR==1) %>% filter(str_detect(espece, "sonnei")) %>%  pull(espece3) %>% table
d %>% filter(metrop == TRUE & GEN_SIR==1) %>%pull(voyage) %>% table(useNA = "a") %>% data.frame
