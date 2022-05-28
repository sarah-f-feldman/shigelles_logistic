##################
# Figure 2 (CMI) #
##################

source("src/00_library.R")
source("src/01_data_management.R")

#--------------------
# 3GC CAZ CRO

CMI <- read.csv2("data/CMI_C3G20190508.csv", stringsAsFactors = FALSE)
Rtested <- nrow(CMI)
totR <- nrow(d[!is.na(d$C3G_SIR) & d$C3G_SIR == 1, ])
CMI$MIC_CAZ_mgperL[str_detect(CMI$MIC_CAZ_mgperL, "256")] <- "256"
CMI$MIC_CRO_mgperL[str_detect(CMI$MIC_CRO_mgperL, "256")] <- "256"
CMI$MIC_CAZ_mgperL <- as.factor(as.numeric(CMI$MIC_CAZ_mgperL))
CMI$MIC_CRO_mgperL <- as.factor(as.numeric(CMI$MIC_CRO_mgperL))
levels(CMI$MIC_CRO_mgperL)[levels(CMI$MIC_CRO_mgperL) == "256"] <- ">256"
levels(CMI$MIC_CAZ_mgperL)[levels(CMI$MIC_CAZ_mgperL) == "256"] <- ">256"

CMI <- CMI %>% group_by(MIC_CRO_mgperL, MIC_CAZ_mgperL) %>% mutate(N_same_combination = n())


#Jitter et CAZ en fonction de CRO
ggplot(CMI, aes(MIC_CRO_mgperL, MIC_CAZ_mgperL)) +
  geom_jitter(width = 0.1) +
  #labs(x = "MIC of CRO(mg/L)",
  labs(x = expression(MIC[CRO]*" "*(mg/L)),
       y =  expression(MIC[CAZ]*" "*(mg/L)),
       title = paste0("MIC of C3G-R isolates; \nN tested = ",Rtested,"; N 3GC-R = ", totR))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(b=10)),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)
  )


ggplot(CMI, aes(y=MIC_CRO_mgperL, x=MIC_CAZ_mgperL, size = N_same_combination)) +
  geom_point() +
  labs(y = "MIC of CRO(mg/L)", 
       x = "MIC of CAZ(mg/L)",
       title = paste0("MIC of CAZ and CRO\nN tested = ",Rtested,"; N = ???"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(t=+20)),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15)
  )



ggplot(CMI, aes(y=MIC_CRO_mgperL, x=MIC_CAZ_mgperL)) +
  geom_jitter(width = 0.1) +
  labs(y= "MIC of CRO(mg/L)", 
       x= "MIC of CAZ(mg/L)",
       title = paste0("MIC of CAZ and CRO\nN tested = ",Rtested,"; N = ???"))+
theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(b=10)),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 15, margin = margin(r = 10)),
      axis.title.x = element_text(size = 15, margin = margin(t = 10))
) 
  
ggplot(CMI, aes(MIC_CRO_mgperL, MIC_CAZ_mgperL)) +
  geom_raster(aes(fill=N_same_combination))+ 
  scale_fill_gradient(low = "light cyan2",high = "steel blue4")+
  theme(panel.background = element_blank())

CMI <- read.csv2("data/CMI_C3G20190508.csv", stringsAsFactors = FALSE)
CMI$MIC_CAZ_mgperL[str_detect(CMI$MIC_CAZ_mgperL, "256")] <- "256"
CMI$MIC_CRO_mgperL[str_detect(CMI$MIC_CRO_mgperL, "256")] <- "256"
CMI$MIC_CAZ_mgperL <- as.numeric(CMI$MIC_CAZ_mgperL)
CMI$MIC_CRO_mgperL <- as.numeric(CMI$MIC_CRO_mgperL)
CMI$MIC_CAZ_mgperL[CMI$MIC_CAZ_mgperL>=32] <- 32
CMI$MIC_CRO_mgperL[CMI$MIC_CRO_mgperL>=32] <- 32

CMI$MIC_CAZ_mgperL <- as.factor(as.numeric(CMI$MIC_CAZ_mgperL))
CMI$MIC_CRO_mgperL <- as.factor(as.numeric(CMI$MIC_CRO_mgperL))
CMI <- CMI %>% group_by(MIC_CRO_mgperL, MIC_CAZ_mgperL) %>% mutate(N_same_combination = n())

#--------------------
# CIP
CMI <- read.csv2("data/CMI_CIP20190508.csv", stringsAsFactors = FALSE)

Rtested <- sum(CMI$N_isolates)
totR <- sum(d$CIP_SIR, na.rm = TRUE)
CMI$MIC_mgperL <- ifelse (CMI$MIC_mgperL %in% "supeq32", "32", CMI$MIC_mgperL)
CMI$MIC_mgperL <- factor(CMI$MIC_mgperL,levels = c("1", "2", "4", "8", "16", "32")) #sup ou egal à 32
ggplot(CMI, aes(x = as.factor(MIC_mgperL), y = N_isolates, group = 1))+
geom_point() + geom_line() +
labs(x = "MIC (mg/L)", 
     y = "Number of isolates",
     title = paste0("MIC of CIP-R isolates\nN tested = ",Rtested,"; N CIP-R = ",totR))+
theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(b=10)),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 15, margin = margin(r = 10)),
      axis.title.x = element_text(size = 15, margin = margin(t = 10))
) 
  

#---------------------------
#AZM
CMI <- read.csv2("data/CMI_AZM20190508.csv")

CMI <- CMI %>% count(CMI_AZM) %>% rename(N_isolates=n, MIC_mgperL = CMI_AZM)
Rtested <- sum(CMI$N_isolates)
totR <- sum(d$AZM_SIR, na.rm = TRUE)

CMI$MIC_mgperL <- factor(CMI$MIC_mgperL,levels = c("24", "32", "48", "64", "96", "128", "192", ">256"))
ggplot(CMI, aes(x = as.factor(MIC_mgperL), y = N_isolates, group = 1))+
  geom_point() + geom_line() +
  labs(x = "MIC (mg/L)", 
       y = "Number of isolates",
       title = paste0("MIC of AZM-R isolates\nN tested = ",Rtested,"; N AZM-R = ",totR))+
theme(plot.title = element_text(hjust = 0.5, size = 20, margin=margin(b=10)),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 15, margin = margin(r = 10)),
      axis.title.x = element_text(size = 15, margin = margin(t = 10))
) 

     
