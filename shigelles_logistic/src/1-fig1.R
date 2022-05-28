###############
#    Fig 1    #
###############



source("src/00_library.R")
source("src/01_data_management.R")





#==============================
# Ciprofloxacin resistance

dlCIP <- dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "CIP_SIR") %>% 
  select(year_real, espece3, status) 

# Reproduction Sophie : 
dlCIP <- dlCIP  %>% 
  group_by(year_real, espece3) %>% count(status) %>% 
  group_by(year_real) %>% mutate(tot = sum(n)) %>% 
  group_by(year_real, status) %>% mutate(n_r_s = sum(n)) %>% 
  # rapport du nb d'espece R une année donneé sur nb total de souches (R et S) cette année là
  mutate(sp_perc_resist = n/tot) %>%
  # verification : la somme des pourcentages de resistance vaut la resistance totale
  filter(status == 1) %>% mutate(yr_perc_resist = n_r_s/tot, cumsum(sp_perc_resist))
dlCIP

g1 <- ggplot(dlCIP, aes(x = factor(year_real), y = sp_perc_resist, fill = espece3))+
  geom_bar(stat = "identity") + 
  labs(y="Prevalence", fill = "CIP") +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05),
                     minor_breaks = seq(0, 0.25, by = 0.01), 
                     labels = scales::percent_format(accuracy = 1), 
                     lim = c(0,0.25)) 
#==============================
# resistance to ceftazidime, cefotaxime or ceftriaxone


dl3GC <-  dl %>% filter(!is.na(status) & metrop == T & year_real %in% c(2005:2016) & ATB == "C3G_SIR") %>% 
  select(year_real, espece3, status)

# Reproduction Sophie : 
dl3GC <- dl3GC  %>% 
  group_by(year_real, espece3) %>% count(status) %>% 
  group_by(year_real) %>% mutate(tot = sum(n)) %>% 
  group_by(year_real, status) %>% mutate(n_r_s = sum(n)) %>% 
  # rapport du nb d'espece R une année donneé sur nb total de souches (R et S) cette année là
  mutate(sp_perc_resist = n/tot) %>%
  # verification : la somme des pourcentages de resistance vaut la resistance totale
  filter(status == 1) %>% mutate(yr_perc_resist = n_r_s/tot, cumsum(sp_perc_resist))
dl3GC

g2 <- ggplot(dl3GC, aes(x = factor(year_real), y = sp_perc_resist, fill = espece3))+
  geom_bar(stat = "identity") + 
  labs(y="Prevalence", fill = "3GC") +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05),
                     minor_breaks = seq(0, 0.25, by = 0.01), 
                     labels = scales::percent_format(accuracy = 1), 
                     lim = c(0,0.25)) 

#==============================
# Azythromycin

dlAZM <- dl %>% 
  filter(!is.na(status) & metrop == T & year_real %in% c(2014,2015,2016) & ATB == "AZM_SIR") %>% 
  select(year_real, espece3, status)

dlAZM$year_real <- factor(dlAZM$year_real, c(2005:2016))



# Reproduction Sophie : 
dlAZM <- dlAZM  %>% 
  group_by(year_real, espece3) %>% count(status) %>% 
  group_by(year_real) %>% mutate(tot = sum(n)) %>% 
  group_by(year_real, status) %>% mutate(n_r_s = sum(n)) %>% 
  # rapport du nb d'espece R une année donneé sur nb total de souches (R et S) cette année là
  mutate(sp_perc_resist = n/tot) %>%
  # verification : la somme des pourcentages de resistance vaut la resistance totale
  filter(status == 1) %>% mutate(yr_perc_resist = n_r_s/tot, cumsum(sp_perc_resist))
dlAZM

g3 <- ggplot(dlAZM, aes(x = year_real, y = sp_perc_resist, fill = espece3))+
  geom_bar(stat = "identity") + 
  labs(y="Prevalence", fill = "AZM") +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05),
                     minor_breaks = seq(0, 0.25, by = 0.01), 
                     labels = scales::percent_format(accuracy = 1), 
                     lim = c(0,0.25)) +
  scale_x_discrete(breaks=as.character(c(2005:2016)),
                   labels=as.character(c(2005:2016)),
                   limits = as.character(c(2005:2016)),
                   drop = FALSE) 






g1+g2+g3 + plot_layout(ncol = 1)

