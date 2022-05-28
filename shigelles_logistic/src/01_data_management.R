#################
# Data cleaning #
#################

#################
# Data cleaning #
#################

library(lubridate)
library(stringr)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)



dall <- readRDS("data/dall.rds")

#--------------------------------
# Age
dall$age <- as.numeric(ifelse(is.na(dall$Date_prel), dall$Date_rec - dall$DN,  dall$Date_prel - dall$DN))
dall$age <- floor(dall$age/365.25)


#--------------------------------
# serotype and species

dall <- dall %>% mutate(espece = ifelse(is.na(espece), Resultats, espece))
dall <- dall %>%  mutate(espece2 = str_extract(dall$espece, "[Ss]higella"), 
                         espece3 = str_extract(dall$espece,"sonnei|nst|non.s.rotypable|[rR]ough|flexneri|non.agglutinable|dysenteriae|boydii|autoagglutinable|sp"),
                         espece3 = tolower(espece3),
                         espece3 = ifelse(grepl("97-10607", Resultats), "dysenteriae", #au lieu de espece3
                                          ifelse(grepl("1c", Resultats), "flexneri",
                                                  ifelse(grepl("Ewing", Resultats), "boydii",
                                                         espece3)))) %>% 
  mutate(espece3 = recode(espece3, "non agglutinable" = "nst"),
         espece3 = gsub("non s.rotypable", "nst", espece3, ignore.case=T),
         espece3 = as.character(espece3))
  
dall <- dall %>% mutate(info_espece_autre = str_extract(espece3, "nst|rough|autoagglutinable|sp"),
         espece3 = ifelse(espece3 %in% c("nst", "rough", "autoagglutinable", "sp"), NA, espece3))

dall$Resultats <- NULL

table(dall$espece3)
table(dall$info_espece_autre)

#------------------------------
#sexe
dall %>% filter(sexe=="I") 
dall <- dall %>% mutate(sexe = ifelse(sexe == "I", NA, sexe))
dall %>%  count(sexe) #ok


#-------------------------------
#year refait a partir de date prelevement (ou date reception si date prelevement incoherente ou inconnue)

write.table(dall %>% mutate(year_prel = year(Date_prel), year_rec = year(Date_rec), diff_date = year_rec - year_prel) %>% filter(diff_date>1 | diff_date < 0),
            sep = "\t", file = "clipboard", row.names = F)

dall <- dall %>% mutate(year_prel = year(Date_prel), year_rec = year(Date_rec), 
                        year_real = ifelse(!is.na(Date_prel), year_prel, ifelse(!is.na(Date_rec), year_rec, year)))

dall %>% filter(year_real<2004) #plus aucun


#-----------------------------------
#SIR

#---------
#modification of SIR for some AZM
new_AZM <- read.csv2("data/Modif AZM sophie_20171025.csv", stringsAsFactors = F)
new_AZM$num_CNR <- as.character(new_AZM$num_CNR) #car character dans dall
new_AZM$AZM_d <- as.numeric(new_AZM$AZM_d) #car character dans dall
new_AZM$AZM_SIR <- NULL


dall <- left_join(dall, new_AZM, by = "num_CNR") %>% 
     mutate(AZM_d = coalesce(AZM_d.y, AZM_d.x), AZM_d.y = NULL, AZM_d.x = NULL) %>% select(colnames(dall))
 
#----
#verification SIR

# diameter
dall %>% group_by(year, metrop) %>% sample_n(4) %>% select(year, metrop,AMX_d : CAZ_d) 

#pourcecntage de diametre NA par annee
data.frame(unique(dall$year), do.call(rbind, lapply(unique(dall$year), function(year){
  myd <- dall[dall$year == year, grepl("_d", names(dall))]
  myd <- apply(myd, 2, is.na)
  round(colSums(myd)/nrow(myd)*100, 1)
})))

#CAZ 10 et CAZ30 ne sont jamais renseignees en meme temps
dall %>% mutate_at(vars(CAZ10_d, CAZ30_d), function(x)!is.na(x)) %>% mutate(res = CAZ10_d + CAZ30_d) %>% select(res) %>% table
dall %>% mutate_at(vars(CAZ10_SIR, CAZ30_SIR), function(x)!is.na(x)) %>% mutate(res = CAZ10_SIR + CAZ30_SIR) %>% select(res) %>% table
dall %>% mutate_at(vars(CAZ10_d, CAZ30_d, CAZ10_SIR, CAZ30_SIR), function(x)!is.na(x)) %>% mutate(res = CAZ10_d + CAZ30_d + CAZ10_SIR + CAZ30_SIR) %>% select(res) %>% table
#CAZ jamais rempli
dall %>% mutate_at(vars(CAZ_d), function(x)!is.na(x)) %>% select(CAZ_d) %>% table

#qq CTX jamais rempli
dall %>% filter(year_real>=2015) %>% mutate_at(vars(CTX_d, CTX_SIR), function(x)!is.na(x)) %>% mutate(res = CTX_d + CTX_SIR) %>% pull(res) %>%  table

# recherche de valeurs aberrantes dans _d et _SIR
bind_rows(d %>% select(contains("_d")) %>% summarise_all(min, na.rm = T),
          d %>% select(contains("_d")) %>% summarise_all(max, na.rm = T))

#----
#Imputation des SIR NA

#les diametres = 0 sont modifies en NA
dall <- dall %>% gather (key = ATB, value = diam, contains("_d")) %>% 
  mutate(diam = ifelse (diam == 0, NA, diam)) %>% 
  spread(key = ATB, value = diam) %>% 
  select(colnames(dall))
# verif
dall %>% select(num_CNR, year_real, espece3, contains("_d"), metrop) %>%
  gather (key = ATB, value = diam, contains("_d")) %>% filter(diam == 0)


ATB_seuil <- read.xlsx("data/tableauATB_seuil_PourSarah.xlsx", sheetIndex = 1,stringsAsFactors=FALSE)
colnames(ATB_seuil)[3] <- "seuil"
ATB_seuil$seuil <- as.numeric(ATB_seuil$seuil)
#je ne retire pas les seuil NA (sinon les SIR deja imputes a tort restent dans la base)
# Pas d'analyse pour la spectinomycine (SPT) car pas de reference fiable pour les seuils
ATB_seuil[ATB_seuil$ATB == "SPT", "seuil"] <- NA
#je retire les ATB qui ne sont pas dans la base
ATB_seuil <- ATB_seuil[paste0(ATB_seuil$ATB, "_d") %in% colnames(dall), ]

#modif : la reference utilisee est celle de l'annee de reception (approximee par year)et non de l'annee de prelevement year_real
#my_year <- 2015
.l <- lapply(unique(dall$year), function(my_year){
  print(my_year)
  dall_tmp <- dall %>% filter(year == my_year) %>% as.data.frame
  ATB_seuil_tmp <- ATB_seuil %>%  filter(annee == my_year) %>% as.data.frame #seuil selon l'annee de reception (approxime par l'annee du fichier ici)
  
  #Si diametre renseigne je recalcule SIR, sinon je garde le SIR renseigne
  calcul_SIR <- function(my_ATB, my_threshold, data = dall_tmp){
    ATB_d <- data[ , paste0(my_ATB,"_d")]
    ATB_SIR <- data[ , paste0(my_ATB,"_SIR")]
    to_impute <- !is.na(ATB_d)
    
    if(is.na(my_threshold)){
      ATB_SIR <- rep(NA, length(ATB_SIR))
    } else {
      if(any(to_impute)){
        val_diam <- ATB_d[to_impute]
        val_SIR_impute <- ifelse(val_diam < as.numeric(as.character(my_threshold)), "R", "S") #il y avait pb avec my_threshold qui etait en caractere
        ATB_SIR[to_impute] <- val_SIR_impute
      }#else sous entendu : on ne touche pas a ATB_SIR s'il n'y a aucune diametre renseigne pour cette annee et cet ATB      
    }
    return(ATB_SIR)
  }
  antibiogramme.df <- mapply(calcul_SIR, my_ATB = ATB_seuil_tmp$ATB, my_threshold = ATB_seuil_tmp$seuil)
  antibiogramme.df <- if(nrow(dall_tmp)==1) data.frame(t(antibiogramme.df)) else data.frame(antibiogramme.df) 
  colnames(antibiogramme.df) <- paste0(colnames(antibiogramme.df), "_SIR")
  dall_tmp <- bind_cols(dall_tmp %>% select(-one_of(colnames(antibiogramme.df))), antibiogramme.df) #Je retire les colonnes de dall qui ont ete recalculees et je les remplace par celles calculees avec calcul_SIR
  dall_tmp <- dall_tmp %>% select(one_of(colnames(dall))) #je remet les colonnes dans l'ordre
  return(dall_tmp)
})

dall <- do.call(rbind, .l)



#----
#Hommogeneisation des SIR non imputes

#tout en majuscule
dall <- dall %>% mutate_at(vars(contains("_SIR")), toupper)
#I devient R
dall <- dall %>% mutate_at(vars(contains("_SIR")), as.character) %>% 
                 mutate_at(vars(contains("_SIR")), funs(recode(., "I" = "R")))
table(unlist(c(dall %>% select(contains("SIR"))))) 
dall <- dall %>% mutate_at(vars(contains("_SIR")), 
                           funs(recode(., "S" = "S", "R" = "R", "S*" = "S", .default = NA_character_))) 


# test <- data.frame(AMX_d = c(17, 25, NA), index = c(1:3))
# test %>% mutate(AMX_SIR = ifelse(AMX_d < 19, "R", ifelse(AMX_d >=19, "S", NA)))

#--------
#Rassembler AMX et AMP / CAZ10 et CAZ30
#jamais renseigner en meme temps
dall %>% mutate(na1 = !is.na(AMX_SIR), na2 = !is.na(AMP_SIR), natot = na1 + na2) %>% count(natot) #pas de 2
dall %>% mutate(na1 = !is.na(CAZ10_SIR), na2 = !is.na(CAZ30_SIR), natot = na1 + na2) %>% count(natot) #pas de 2
dall <- dall %>% mutate(AMXAMP_SIR = coalesce(AMX_SIR, AMP_SIR), 
                        CAZ1030_SIR = coalesce(CAZ10_SIR, CAZ30_SIR),
                        AMX_SIR = NULL, AMP_SIR = NULL, CAZ10_SIR = NULL, CAZ30_SIR = NULL)



#-------------------------------
#origine prelevement a homogeneiser
#codes inconnus dans epi : C, E, F, I : C = collectivite? colonie de vacances? creche? E = Ecole?
#voyage pas du tout homogene, ça va etre complique a analyser.

#---------------------------------
# recherche et suppression des doublons/patients avec plusieurs prelevement
length(dall$num_CNR)
length(unique(dall$num_CNR))

dall %>% mutate (name_dupl = paste0(NOM, PRENOM)) %>% 
  group_by(name_dupl) %>% mutate(np = n()) %>% 
  #count(name_dupl) %>% 
  filter(np>1 & name_dupl!="NANA") %>% 
  #pull(name_dupl) %>% unique #55 doublons
  arrange(name_dupl) %>% View

#vrais doublons = meme nom de patient, meme serotype et meme antibiogramme.

#--------------------------------
#recoder les voyages (pour les souches resistantes a C3G, CIP, AZM)
dall %>% filter(CAZ10_SIR == "R" | CAZ30_SIR == "R" | CTX_SIR =="R" | CRO_SIR =="R" | CIP_SIR =="R" | AZM_SIR =="R") %>% pull(voyage) %>% table
dall %>% filter(CAZ10_SIR == "R" | CAZ30_SIR == "R" | CTX_SIR =="R" | CRO_SIR =="R" | CIP_SIR =="R" | AZM_SIR =="R") %>% pull(voyage_continent) %>% table

dall <- dall %>% mutate(voyage = tolower(voyage),
                  voyage = ifelse(grepl("inde.*", voyage), "inde", voyage),
                  voyage = gsub("pas de voyage", "aucun", voyage),
                  voyage = gsub("e", "e", voyage),
                  voyage_continent = ifelse(voyage %in% c("allemagne", "espagne", "grece", "portugal", "suisse", "iles canaries", "italie"), "europe", voyage),
                  voyage_continent = ifelse(voyage %in% c("asie", "thailande", "sri lanka", "cambodge", "chine", "inde", "myanmar", "nepal", "nepal,ind", "mongolie", "viet nam"), "asie", voyage_continent),
                  voyage_continent = ifelse(voyage %in% c("tunisie", "maroc", "algerie"), "maghreb", voyage_continent))

#---------------------------------
#recoder ville ou hospitalise

#ville ou hospitalise
#NB : ce n'est pas la colonne epi (epi remplie que si TIAC)
#SI L.A.M ou LABM ou CM: c'est ville
#Si Hôp, hosp, IHU, Hop, CHU, CH, GHU, GH, CHRU, CHR, CHI, CHG
vec_lab <- dall$labo
vec_lab[grepl("^CH", vec_lab)] <- "hospitalise"
vec_lab[grepl("^GH", vec_lab)] <- "hospitalise"
vec_lab[grepl("GCS", vec_lab)] <- "hospitalise"
vec_lab[grepl("^H", vec_lab)] <- "hospitalise"
vec_lab[grepl("^C.H.", vec_lab)] <- "hospitalise"
vec_lab[grepl("^ho", vec_lab, ignore.case = T)] <- "hospitalise"
vec_lab[grepl("h.sp", vec_lab, ignore.case = T)] <- "hospitalise"
vec_lab[grepl("h.pit", vec_lab, ignore.case = T)] <- "hospitalise"
vec_lab[grepl("^hô", vec_lab, ignore.case = T)] <- "hospitalise"
vec_lab[grepl("facult.", vec_lab, ignore.case = T)] <- "hospitalise"

vec_lab[grepl("institu", vec_lab, ignore.case = T)] <- "institut"
vec_lab[grepl("IP", vec_lab, ignore.case = F)] <- "institut"

vec_lab <- ifelse(!vec_lab %in% c("hospitalise", "institut") & !is.na(vec_lab), "ville", vec_lab)
table(vec_lab)

dall$labo <- vec_lab
#---------------------------
#Je supprime la colonne SPT_SIR car fait bugger les analyses (car entierement NA)
dall$SPT_SIR <- NULL



#saveRDS(dall, "data/dall_dmfini.rds")
#saveRDS(dall, "data/dall_dmfini20171017.rds") #Resultats corriges, dm voyage et labo, rassemblement AMX et AMP, CAZ10 et CAZ30, debug de la fonction SIR (le seuil n'etait pas numerique...) 
#saveRDS(dall, "data/dall_dmfini20171018.rds") #97-10607 a partir de Resultats et non de espece3. Retrait de la colonne Resultats (redondant avec espece)
#saveRDS(dall, "data/dall_dmfini20171023.rds") #diametre = 0 remplace par NA
#saveRDS(dall, "data/dall_dmfini20171027.rds") #correction de quelques diametre et fonction calcul_SIR refaite : prend year au lieu de year_real et met SIR NA si seuil NA, CAZ10, 30 AMX et AMP supprime (on ne garde que les fusions CAZ1030 et AMXAMP) 
saveRDS(dall, "data/dall_dmfini20171128.rds") #correction de de diametre CAZ30 et CTX pour 7 patients
#write.csv2(dall, file = "data/dall20170929.csv")
# write.csv2(dall, file = "data/dall20171017.csv") #Resultats corriges
# write.csv2(dall, file = "../Pour Sophie 20171018/tableau20171018.csv") #Resultats corriges
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

#---------------------------
#data management suite

d <- dall
#d <- readRDS("data/dall_dmfini20171027.rds")
#d <- readRDS("data/dall_dmfini20171128.rds")

# Sophie 20171018 : Virer les 2004
d <- d %>% filter(year_real > 2004)
#je mets les age negatif NA
d <- d %>% mutate(age = ifelse(age<0, NA, age), DN = as_date(ifelse(age<0, NA, DN)))

#mail sophie 29/11 : je supprime info AZM pour date de reception avant avril 2014
d[d$Date_rec < as_date("2014-04-01"), c("AZM_d", "AZM_SIR")] <- c(NA, NA)


#implementation de C3G_R 
d <- d %>% mutate_at(vars(contains("SIR")), funs(recode(., "S" = F , "R" = T))) 
#si au moins un TRUE resultat est TRUE. Si tout FALSE ou tout NA resultat FALSE
d$C3G_SIR <- apply(d[ ,c("CAZ1030_SIR", "CTX_SIR", "CRO_SIR")], 1, any, na.rm=T)
#Je veux transformer les C3G_SIR en NA si tout est NA 
d$C3G_SIR1noNA <- apply(apply(d[ ,c("CAZ1030_SIR", "CTX_SIR", "CRO_SIR")], 2, function(x)!is.na(x)), 1, any)#si que des NA alors any est FALSE
d$C3G_SIR[!d$C3G_SIR1noNA] <- NA #si que des NA alors C3G_SIR est NA 
d$C3G_SIR1noNA <- NULL

# S et R en 0,1
d <- d %>% mutate_at(vars(contains("SIR")), funs(as.numeric)) 
#creation colonne inde
d$inde <- ifelse(d$voyage == "inde" & !is.na(d$voyage), 1, 0)

#saveRDS(d, "data/dall_dmfini20171122.rds") 
#saveRDS(d, "data/dall_dmfini20171129.rds") # pas bon
saveRDS(d, "data/dall_dmfini20171130.rds") 


#---------------------------
#data management suite

#d <- readRDS("data/dall_dmfini20171130.rds")

# departement_bis : departement, CP patient ou CP labo 
d$departement_bis <- d$departement

d <- d %>% 
  mutate(CP_patlab = coalesce(CP_patient, CP_labo),
         CP_patlab = str_sub(CP_patlab, 1, 3)) 

d <- d %>% mutate(departement_bis = ifelse(metrop == FALSE & is.na(departement), CP_patlab, departement))
# filter(metrop == FALSE & is.na(departement)) %>% pull(departement_bis)
# d$departement_bis[d$metrop == FALSE & is.na(d$departement)] <- departement_bis
#je modifie le 96 (96 appartient au protectorat tunisien qui n'existe plus depuis 1956): je prend CP du labo
d %>% filter(departement_bis == 96)
d<-d %>% mutate(departement_bis = ifelse(departement_bis == 96, CP_patlab, departement_bis))
table(d$departement_bis[d$metrop == FALSE], useNA = "a")
#Je modifie 978 en 974 : apres verification ce ne sont que des labos a saint paul dont le CP est 97460 (et le cedex est 97863)
d <- d %>% mutate(departement_bis = ifelse(departement_bis == 978, 974, departement_bis))
table(d$departement_bis[d$metrop == FALSE], useNA = "a")
#si CP patient etait renseigne, alros je privilegie CP patient
d <- d %>% mutate(departement_bis = ifelse(!is.na(CP_patient), CP_patlab, departement_bis))
table(d$departement_bis[d$metrop == FALSE], useNA = "a")

# saveRDS(d, "data/bdd_shigelles_20180101.rds")
# 
# #---------------------------
# #data management suite
# d <- readRDS("data/bdd_shigelles_20180101.rds")

#----------

Voyage
#Correction a la main des voyages
#table(d$voyage) %>% data.frame %>% write.csv2(., "data/voyages_continent_toclean.csv", row.names = FALSE)
voy <- read.csv2("data/voyages_continent_cleaned.csv", stringsAsFactors = FALSE)
#merger voy et d
d <- left_join(d %>% select(-voyage_continent), voy %>% select(-Freq), by = "voyage")


#variable voyage et continent
#je rajoute cat?gorie missing
d$continent[d$continent %in% "none"] <- NA
#je reorganise les levels de continent
d$continent <- as.factor(d$continent)
d$continent <- factor(d$continent, levels = c("africa", "asia", "europe", "north america", 
                                              "south america", "oceania", "unknown"))
d$voyage_01 <- factor(d$voyage_01)#30/03

# voyage01_NA avec non precise et non renseigne qui deviennent NA 
d$voyage_01_NA <- d$voyage_01
d$voyage_01_NA[d$voyage %in% c("non precise", "non renseigne")] <- NA
# voyage01_m avec non precise, non renseigne et NA dans une classe missing
d$voyage_01_m <- as.character(d$voyage_01_NA)
d$voyage_01_m[is.na(d$voyage_01_m)] <- "missing_cat"
d$voyage_01_m <- as.factor(d$voyage_01_m)

d$inde <- as.character(d$voyage_01)
d$inde <- ifelse(d$voyage %in% "inde", "3", d$inde)
d$inde <- factor(d$inde, labels = c("pas de voyage", "autres voyages", "inde"))

d$inde_ref <- d$inde
d$inde_ref <- relevel(d$inde_ref, ref = "inde")

d$inde_NA <- as.character(d$voyage_01_NA)
d$inde_NA <- ifelse(d$voyage %in% "inde", "3", d$inde_NA)
d$inde_NA <- factor(d$inde_NA, labels = c("pas de voyage", "autres voyages", "inde"))

d$inde_m <- as.character(d$voyage_01_m)
d$inde_m <- ifelse(d$voyage %in% "inde", "3", d$inde_m)
d$inde_m <- factor(d$inde_m, labels = c("pas de voyage", "autres voyages", "inde", "missing_cat"))

d$asie <- as.character(d$voyage_01)
d$asie <- ifelse(d$continent %in% "asia", "3", d$asie)
d$asie <- factor(d$asie, labels = c("pas de voyage", "autres voyages", "asie"))

d$asie_NA <- as.character(d$voyage_01_NA)
d$asie_NA <- ifelse(d$continent %in% "asia", "3", d$asie_NA)
d$asie_NA <- factor(d$asie_NA, labels = c("pas de voyage", "autres voyages", "asie"))

d$asie_m <- as.character(d$voyage_01_m)
d$asie_m <- ifelse(d$continent %in% "asia", "3", d$asie_m)
d$asie_m <- factor(d$asie_m, labels = c("pas de voyage", "autres voyages", "asie", "missing_cat"))


#----------
#decoupage des annees
d$year_real_cut <- NA
d$year_real_cut[d$year_real %in% 2005:2008] <- "2005_2008"
d$year_real_cut[d$year_real %in% 2009:2012] <- "2009_2012"
d$year_real_cut[d$year_real %in% 2013:2016] <- "2013_2016"
#table(d$year_real_cut)

#-----------
#decoupage des ages
d$age_cut <- NA 
d$age_cut[d$age < 15] <- "< 15"
d$age_cut[d$age %in% 15:49] <- "15-49"
d$age_cut[d$age > 49] <- "50 +"
table(d$age_cut, useNA = "a")

#------------
# variable age sexe

d$age_sexe <-  NA
d$age_sexe[d$age<15 & !is.na(d$age) & d$sexe %in% "F" & !is.na(d$sexe)] <- "F_14less"
d$age_sexe[d$age>=15 & d$age<50 & !is.na(d$age) & d$sexe %in% "F" & !is.na(d$sexe)] <- "F_15to49"
d$age_sexe[d$age>=50 & !is.na(d$age) & d$sexe %in% "F" & !is.na(d$sexe)] <- "F_50plus"
d$age_sexe[d$age<15 & !is.na(d$age) & d$sexe %in% "M" & !is.na(d$sexe)] <- "M_14less"
d$age_sexe[d$age>=15 & d$age<50 & !is.na(d$age) & d$sexe %in% "M" & !is.na(d$sexe)] <- "M_15to49"
d$age_sexe[d$age>=50 & !is.na(d$age) & d$sexe %in% "M" & !is.na(d$sexe)] <- "M_50plus"

#---------------
# serotype flexneri

d$serotype_flexneri <- d$serotype
d$serotype_flexneri[d$espece3 != "flexneri" | 
                      (! d$serotype %in% c("1b","2a", "3a") & !is.na(d$serotype)) | 
                      is.na(d$serotype)] <- NA
d$serotype_flexneri[d$espece3 == "flexneri" & ! d$serotype %in% c("1b","2a", "3a") & !is.na(d$serotype)] <- "others" 
d$serotype_flexneri[d$espece3 == "flexneri" & is.na(d$serotype)] <- "missing"
d$serotype_flexneri <- factor(d$serotype_flexneri)


# certaines info flexneri etait dans espece et serotype n'etait pas missing
d %>% filter(serotype_flexneri %in% c("missing")) %>% select(espece, espece3, serotype_flexneri, serotype)
d %>% filter(serotype_flexneri %in% c("missing")) %>% 
  select(espece, espece3, serotype_flexneri, serotype) %>% 
  pull(espece) %>% table
table(d$serotype_flexneri)
d[d$espece3 %in% "flexneri" & str_detect(d$espece, "1b") & d$serotype_flexneri == "missing", "serotype_flexneri"] <- "1b"
d[d$espece3 %in% "flexneri" & str_detect(d$espece, "2a") & d$serotype_flexneri == "missing", "serotype_flexneri"] <- "2a"
d[d$espece3 %in% "flexneri" & str_detect(d$espece, "3a") & d$serotype_flexneri == "missing", "serotype_flexneri"] <- "3a"
d[d$espece3 %in% "flexneri" & d$serotype_flexneri %in% "missing", "serotype_flexneri"] <- "others"
table(d$serotype_flexneri)
d$serotype_flexneri <- factor(d$serotype_flexneri) #refaire pour supprimer missing qui vaut maintenant 0


###################################
# Long formqt
#base de donn?e en format long
#dl <- d %>% select(num_CNR, year_real, year_real_cut, espece3, contains("SIR"), metrop, age_cut, sexe, voyage_01, serotype_flexneri, serotype) %>%
#20180322
dl <- d %>% select(year_real, year_real_cut, espece3, contains("SIR"), metrop, age, age_cut, sexe, voyage_01, 
                   inde, #20190108
                   inde_NA, inde_m, inde_ref, #20190116
                   asie, asie_NA, asie_m, #20190116
                   voyage_01_NA, voyage_01_m, #20190115
                   serotype_flexneri, serotype) %>%
  gather (key = ATB, value = status, contains("SIR"))
dl <- dl %>% filter(!is.na(status))


#####################################
# dataframes pour Table 2 de l'article : reg logistique AZM 2014-2016 sonnei et flexneri metropole
dAZM <- dl %>% filter(metrop == TRUE & ATB =="AZM_SIR" & espece3 %in% c("sonnei", "flexneri"))
dAZM$shigelle <- ifelse(dAZM$espece3 == "sonnei", "sonnei", paste0("flexneri_",dAZM$serotype_flexneri))
dAZM$year_real <- as.factor(dAZM$year_real)
#varAZM <- c("age_cut", "sexe", "voyage_01", "espece3") #commente 20190115
varAZM <- c("age_cut", "sexe", "voyage_01_m", "espece3", "year_real") #20190122
#varAZM <- c("age_cut", "sexe", "voyage_01", "voyage_01_m", "voyage_01_NA", "espece3") # 20190115
dAZM_mod <- na.omit(dAZM[ ,c("status", varAZM)])

########################################
# dataframes pour Table 3 de l'article : reg logistique C3G 2005-2016 sonnei metropole
dC3G <- dl %>% filter(metrop == TRUE & ATB =="C3G_SIR" & espece3 %in% c("sonnei"))
#varC3G <- c("age_cut", "sexe", "voyage_01", "year_real_cut") #commente 20190115
varC3G <- c("age_cut", "sexe", "voyage_01_m", "year_real_cut") #20190122
#varC3G <- c("age_cut", "sexe", "voyage_01", "voyage_01_m", "voyage_01_NA", "year_real_cut") # 20190115
dC3G_mod <- na.omit(dC3G[ ,c("status", varC3G)])

#######################################
# dataframe pour reg logistique CIP 2005-2016 sonnei flexneri metropole
dCIP <- dl %>% filter(metrop == TRUE & ATB =="CIP_SIR" & espece3 %in% c("sonnei", "flexneri"))
#varCIP <- c("age_cut", "sexe", "inde", "voyage_01", "year_real_cut", "espece3") #commente 20190116
#varCIP <- c("age_cut", "sexe", "inde", "inde_ref", "asie", "voyage_01", "year_real_cut", "espece3") #20190116
varCIP <- c("age_cut", "sexe", "inde_m", "year_real_cut", "espece3") #20190122
#varCIP <- c("age_cut", "sexe", "inde", "voyage_01", "voyage_01_m", "voyage_01_NA", "year_real_cut", "espece3") #20190115
dCIP_mod <- na.omit(dCIP[ ,c("status", varCIP)])


# 2018 05 01 essai en changeant la reference pour les annees :
dC3G_mod$year_real_cut <- as.factor(dC3G_mod$year_real_cut)
#dC3G_mod$year_real_cut <- relevel(dC3G_mod$year_real_cut, ref = "2009_2012")
#dC3G_mod$year_real_cut <- relevel(dC3G_mod$year_real_cut, ref = "2013_2016")
dC3G_mod$year_real_cut <- relevel(dC3G_mod$year_real_cut, ref = "2005_2008") #25/02/2019
# fin 2018 05 01

# 2018 05 01
dCIP_mod$year_real_cut <- as.factor(dCIP_mod$year_real_cut)
#dCIP_mod$year_real_cut <- relevel(dCIP_mod$year_real_cut, ref = "2013_2016")
dCIP_mod$year_real_cut <- relevel(dCIP_mod$year_real_cut, ref = "2005_2008") #25/02/2019
# fin # 2018 05 01
