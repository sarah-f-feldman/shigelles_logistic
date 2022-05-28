###################
#     Figure 3    #
###################


source("src/00_library.R")
source("src/01_data_management.R")

#2019 01 17
#liste des ATB d'interet dans l'ordre
ATB_ordered <- c("AMX/AMP_SIR", "3GC_SIR", "NAL_SIR", "CIP_SIR", "STR_SIR", "AKN_SIR", "GEN_SIR", "SMX_SIR", "TMP_SIR",
                 "SXT_SIR", "CHL_SIR", "TET_SIR", "AZM_SIR")
#fin 2019 01 17

#----------------------------------
#En 2005-2006, metrop, toute especes confondues
# nombre de souche resistantes sur nombre de souches testees

s2 <- dl %>% filter(!is.na(status) & year_real %in% 2005:2006 & metrop == TRUE) %>%
  group_by(ATB, espece3, status) %>% summarise(n = n()) %>%  group_by(ATB, espece3) %>% 
  mutate(is_status1 = sum(status),
         row_to_keep = ifelse(is_status1 == 0, 1, status)
  ) %>% 
  mutate(sum_n = sum(n),
         n = ifelse(is_status1 == 0, 0, n),
         perc_R = (n/sum_n)+0.00001) %>% 
  filter(row_to_keep == "1" & !is.na(espece3)) %>% arrange (ATB)


#Je renomme 2 niveaux
s2$ATB[s2$ATB == "C3G_SIR"] <- "3GC_SIR"
s2$ATB[s2$ATB == "AMXAMP_SIR"] <- "AMX/AMP_SIR"

#Je filtre sur ATB d'interet
s2 <- s2 %>% filter(ATB %in% ATB_ordered)
#Je met les niveaux dans l'ordre de ATB d'interet
s2$ATB <- factor(s2$ATB, c(ATB_ordered))


p0 <- data.frame(x=s2$ATB,
                strain = s2$espece3,
                n = s2$n,
                sum_n = s2$sum_n,
                value = round(s2$perc_R*100,1),
                value_bin = s2$perc_R*100 < 50,
                value_cut = cut(s2$perc_R*100, breaks = c(0, 25, 50, 75, 100))) %>% na.omit


p <- p0 %>%
  ggplot(aes(x=x, y = strain, size = value, fill = value_cut,
             #width =width
             )) +
  geom_point(shape = 21, alpha = 1) + 
  theme_minimal() +
  scale_fill_manual(values = c("#9ad0f3", "#0072B2","darkorange1", "brown3"), label = c("0-25","25-50","50-75","75-100"))+
  scale_size(range = c(1, 15)) +
  geom_text(aes(label =  paste0(n,"/",sum_n)),  vjust = -2, size = 4, colour = "black") +
  scale_x_discrete(breaks=ATB_ordered,
                   labels=gsub("_SIR", "", ATB_ordered),
                   limits = ATB_ordered,
                   position = "top", drop = FALSE) +
  
  theme(axis.text.x = element_text(angle=45, vjust = -3), axis.title = element_blank())
M1 <- p +  guides(size = "none") +
  guides(fill = guide_legend("Resistance, %", override.aes = list(size = c(1,5,10,15)))) +
  labs(title = "Shigella antibioresistance, 2005-2006, Metropole") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(t=+20)),
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        axis.ticks.length = unit(.85, "cm"))


M1
#----------------------------------
#En 2015-2016, metrop, toute especes confondues
# nombre de souche resistantes sur nombre de souches testees

s2 <- dl %>% filter(!is.na(status) & year_real %in% 2015:2016 & metrop == TRUE) %>%
  group_by(ATB, espece3, status) %>% summarise(n = n()) %>%  group_by(ATB, espece3) %>% 
  mutate(is_status1 = sum(status),
         row_to_keep = ifelse(is_status1 == 0, 1, status)
  ) %>% 
  mutate(sum_n = sum(n),
         n = ifelse(is_status1 == 0, 0, n),
         perc_R = (n/sum_n)+0.00001) %>% 
  filter(row_to_keep == "1" & !is.na(espece3)) %>% arrange (ATB)

#Je renomme 2 niveaux
s2$ATB[s2$ATB == "C3G_SIR"] <- "3GC_SIR"
s2$ATB[s2$ATB == "AMXAMP_SIR"] <- "AMX/AMP_SIR"

#Je filtre sur ATB d'interet
s2 <- s2 %>% filter(ATB %in% ATB_ordered)
#Je met les niveaux dans l'ordre de ATB d'interet
s2$ATB <- factor(s2$ATB, c(ATB_ordered))


#data.frame(x = gsub("_SIR", "", s2$ATB),
p <- data.frame(x = s2$ATB,
                strain = s2$espece3,
                n = s2$n,
                sum_n = s2$sum_n,
                value = round(s2$perc_R*100,1),
                value_bin = s2$perc_R*100 < 50,
                value_cut = cut(s2$perc_R*100, breaks = c(0, 25, 50, 75, 100))) %>% na.omit %>%
  ggplot(aes(x=x, y = strain, size = value, fill = value_cut)) +
  geom_point(shape = 21, alpha = 1) + theme_minimal() +
  scale_fill_manual(values = c("#9ad0f3", "#0072B2","darkorange1", "brown3"), label = c("0-25","25-50","50-75","75-100"))+
  scale_size(range = c(1, 15)) +
  geom_text(aes(label =  paste0(n,"/",sum_n)),  vjust = -2, size = 4, colour = "black") +
  scale_x_discrete(breaks=unique(s2$ATB),
                   labels=gsub("_SIR", "", unique(s2$ATB)),
                   position = "top") +
  theme(axis.text.x = element_text(angle=45, vjust = -3), axis.title = element_blank())
M2 <- p +  guides(size = "none") +
  guides(fill = guide_legend("Resistance, %", override.aes = list(size = c(1,5,10,15)))) +
  labs(title = "Shigella antibioresistance, 2015-2016, Metropole") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(t=+20)),
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        axis.ticks.length = unit(.85, "cm"))


#----------------------------------
#En 2005-2006, DOMTOM, toute especes confondues
# nombre de souche resistantes sur nombre de souches testees

s2 <- dl %>% filter(!is.na(status) & year_real %in% 2005:2006 & metrop == FALSE) %>%
  group_by(ATB, espece3, status) %>% summarise(n = n()) %>%  group_by(ATB, espece3) %>% 
  mutate(is_status1 = sum(status),
         row_to_keep = ifelse(is_status1 == 0, 1, status)
  ) %>% 
  mutate(sum_n = sum(n),
         n = ifelse(is_status1 == 0, 0, n),
         perc_R = (n/sum_n)+0.00001) %>% 
  filter(row_to_keep == "1" & !is.na(espece3)) %>% arrange (ATB)

#Je renomme 2 niveaux
s2$ATB[s2$ATB == "C3G_SIR"] <- "3GC_SIR"
s2$ATB[s2$ATB == "AMXAMP_SIR"] <- "AMX/AMP_SIR"

#Je filtre sur ATB d'interet
s2 <- s2 %>% filter(ATB %in% ATB_ordered)
#Je met les niveaux dans l'ordre de ATB d'interet
s2$ATB <- factor(s2$ATB, c(ATB_ordered))


#data.frame(x = gsub("_SIR", "", s2$ATB),
p <- data.frame(x = s2$ATB,
                strain = s2$espece3,
                n = s2$n,
                sum_n = s2$sum_n,
                value = round(s2$perc_R*100,1),
                value_bin = s2$perc_R*100 < 50,
                value_cut = cut(s2$perc_R*100, breaks = c(0, 25, 50, 75, 100))) %>% na.omit %>%
  ggplot(aes(x=x, y = strain, size = value, fill = value_cut)) +
  geom_point(shape = 21, alpha = 1) + theme_minimal() +
  #scale_fill_brewer(palette = 3, direction = 2)+
  scale_fill_manual(values = c("#9ad0f3", "#0072B2","darkorange1", "brown3"), label = c("0-25","25-50","50-75","75-100"))+
  scale_size(range = c(1, 15)) +
  geom_text(aes(label =  paste0(n,"/",sum_n)),  vjust = -2, size = 4, colour = "black") +
  scale_x_discrete(breaks=ATB_ordered,
                   labels=gsub("_SIR", "", ATB_ordered),
                   limits = ATB_ordered,
                   position = "top", drop = FALSE) +
  theme(axis.text.x = element_text(angle=45, vjust = -3), axis.title = element_blank())
D1 <- p +  guides(size = "none") +
  guides(fill = guide_legend("Resistance, %", override.aes = list(size = c(1,5,10)))) +
  labs(title = "Shigella antibioresistance, 2005-2006, Dom Tom") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(t=+20)),
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        axis.ticks.length = unit(.85, "cm"))

#----------------------------------
#En 2015-2016, DOMTOM, toute especes confondues
# nombre de souche resistantes sur nombre de souches testees

s2 <- dl %>% filter(!is.na(status) & year_real %in% 2015:2016 & metrop == FALSE) %>%
  group_by(ATB, espece3, status) %>% summarise(n = n()) %>%  group_by(ATB, espece3) %>% 
  mutate(is_status1 = sum(status),
         row_to_keep = ifelse(is_status1 == 0, 1, status)
  ) %>% 
  mutate(sum_n = sum(n),
         n = ifelse(is_status1 == 0, 0, n),
         perc_R = (n/sum_n)+0.00001) %>% 
  filter(row_to_keep == "1" & !is.na(espece3)) %>% arrange (ATB)

#2019 01 17
#Je renomme 2 niveaux
s2$ATB[s2$ATB == "C3G_SIR"] <- "3GC_SIR"
s2$ATB[s2$ATB == "AMXAMP_SIR"] <- "AMX/AMP_SIR"

#Je filtre sur ATB d'interet
s2 <- s2 %>% filter(ATB %in% ATB_ordered)
#Je met les niveaux dans l'ordre de ATB d'interet
s2$ATB <- factor(s2$ATB, c(ATB_ordered))
#fin 2019 01 17


p <- data.frame(x = s2$ATB,
                strain = s2$espece3,
                n = s2$n,
                sum_n = s2$sum_n,
                value = round(s2$perc_R*100,1),
                value_bin = s2$perc_R*100 < 50,
                value_cut = cut(s2$perc_R*100, breaks = c(0, 25, 50, 75, 100))) %>% na.omit %>%
  ggplot(aes(x=x, y = strain, size = value, fill = value_cut)) +
  geom_point(shape = 21, alpha = 1) + theme_minimal() +
  scale_fill_manual(values = c("#9ad0f3", "#0072B2","darkorange1", "brown3"), label = c("0-25","25-50","50-75","75-100"))+
  scale_size(range = c(1, 15)) +
  geom_text(aes(label =  paste0(n,"/",sum_n)),  vjust = -2, size = 4, colour = "black") +
  scale_x_discrete(breaks=ATB_ordered,
                   labels=gsub("_SIR", "", ATB_ordered),
                   limits = ATB_ordered,
                   position = "top", drop = FALSE) +
  theme(axis.text.x = element_text(angle=45),
        axis.title = element_blank())
D2 <- p +  guides(size = "none") +
  guides(fill = guide_legend("Resistance, %", override.aes = list(size = c(1,5,10,15)))) +
  labs(title = "Shigella antibioresistance, 2015-2016, Dom Tom") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(t=+30)),
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 15),
        axis.ticks.length = unit(.85, "cm"))


#------------------------------

M1  
M2 
D1 + D2 + plot_layout(ncol = 1)
