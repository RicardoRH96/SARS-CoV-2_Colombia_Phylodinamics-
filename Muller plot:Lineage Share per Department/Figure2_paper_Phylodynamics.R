### Script modified on June 13 2022 by Cinthy 


library(ggplot2)
library(data.table)
library(purrr)
library(ggpubr)
library(dplyr)
library(maps)
library(viridis)
library(ggmap)
library(raster)
library(sp)
library(ggsci)
library(readxl)
library(lubridate)
library(plotly)
library(htmlwidgets)
library(extrafont)
library(tidyverse)


setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Database/rawdataset/Colombia")


colombia <- as.data.frame(fread("Colombia_allseq_130622_CJ.tsv"))
#To know the total of lineages
#length(unique(colombia$lineage))

#colombia$date <- as.numeric(colombia$date)
colombia$date <- as.Date(colombia$date)


colombia$date <- format(colombia$date, "%Y-%m")

colombia <- colombia %>% distinct(strain, .keep_all = TRUE)

totalgenomes <- as.data.frame(table(colombia$division))
#cundinamarca = 261+ 3217
totalgenomes[totalgenomes == 261] <- 3478

totalgenomes <- totalgenomes %>% filter(Var1 != "Bogota")
totalgenomes <- totalgenomes %>% filter(Var1 != "Na")
id <- c(1,12,23,
        27,28,29,
        30,31,32,
        2,3,4,
        5,6,7,
        10,8,9,
        11,13,14,
        15,16,17,
        18,19,20,
        21,22,24,
        25,26)
totalgenomes <- cbind(totalgenomes,id) #gives me a number of lenghts error because of differing number of rows causing by mispellings

#homogenize department names

totalgenomes$Var1 <- gsub("AtlÃ¡ntico", "Atlantico", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Bogota D.C.", "Bogota", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("BogotÃ¡", "Bogota", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("BolÃ­var", "Bolivar", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("BoyacÃ¡", "Boyaca", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("CÃ³rdoba", "Cordoba", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("CALDAS", "Caldas", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("CaquetÃ¡", "Caqueta", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("CAUCA", "Cauca", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Distrito Capital", "Bogota", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("GuanÃ­a", "Guainia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Guiania", "Guainia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("NariÃ±o", "Narino", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Norte De Santander", "Norte de Santander", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("QuindÃ­o", "Quindio", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("San AndrÃ©s y Providencia", "San Andres y Providencia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("San Andres", "San Andres y Providencia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("San Andres islas", "San Andres y Providencia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("San Andres Y Providencia", "San Andres y Providencia", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Valle del cauca", "Valle del Cauca", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("Valle Del Cauca", "Valle del Cauca", totalgenomes$Var1)
totalgenomes$Var1 <- gsub("VaupÃ©s", "Vaupes", totalgenomes$Var1)










names(totalgenomes) <- c("states","Genomes","id")

##############################################################
#plot all lineages vs VOC
ggplot(colombia, aes(lineage_b, fill=VOC)) +
  geom_bar() +
  coord_flip()+
  theme_minimal()

#############################################
#plot VOC through time 

lin_share <- as.data.frame(table(colombia$date,colombia$VOC))
names(lin_share) <- c("date","VOC","Count")


data <- lin_share %>%
  group_by(date, VOC) %>%
  summarise(n = sum(Count)) %>%
  mutate(Percentage = n/sum(n))


data_num_date <- data %>%
  mutate_at(vars(date), ym)
  
class(data_num_date$date)
time <- ggplot(data = data_num_date, mapping = aes(x=date, y=Percentage, fill = VOCI)) +
  stat_smooth(se=FALSE, geom="area",
              method = 'loess',
              span = 0.1,aes(fill=VOC), 
              position = "stack")+
  #geom_area(alpha=0.6, size=0.01)+
  #geom_histogram(position = "fill") +
  #theme_minimal()+
  #scale_fill_brewer(palette="Spectral")+
  #scale_color_npg()+
  #scale_fill_npg()+
  #theme(text = element_text(size=20, family = "Times New Roman"))+
  xlab("Time") +
  ylab ("Prevalence (%)")+ 
 # theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b-%Y")+
  #wave1 Aug 2020 214/365
  geom_vline(xintercept =as.numeric(data_num_date$date[66]), colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave2 Jan 2021  25/365
  geom_vline(xintercept =as.numeric(data_num_date$date[120]), colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave3 Jun 2021 152/365
  geom_vline(xintercept =as.numeric(data_num_date$date[170]), colour = "grey", lty=1, size=3, alpha=0.5)+
  #wave4 Jan 2022
  geom_vline(xintercept =as.numeric(data_num_date$date[243]), colour = "grey", lty=1, size=3, alpha=0.5)+
  #lockdown 25 march 2020
  geom_vline(xintercept =as.numeric(data_num_date$date[12]), colour = "black", lty=2)+
  #end of lockdown 18 june 2020
  geom_vline(xintercept =as.numeric(data_num_date$date[40]), colour = "black", lty=2)+
  #economic reactivation 10 Agust 2020
  geom_vline(xintercept =as.numeric(data_num_date$date[56]), colour = "black", lty=2)+
  #Reopening of domestic and international flights september 2020
  geom_vline(xintercept =as.numeric(data_num_date$date[77]), colour = "black", lty=2)+
  #vaccination phase 1 Febraury 2021
  geom_vline(xintercept =as.numeric(data_num_date$date[125]), colour = "black", lty=2)+
  #vaccination phase 2 June 2021
  geom_vline(xintercept =as.numeric(data_num_date$date[243]), colour = "black", lty=2)+
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#64CC80",
                                    "B.1.1.348" = "#80A680",
                                    "B.1.420" = "magenta",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
  scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#64CC80",
                                  "B.1.1.348" = "#80A680",
                                  "B.1.420" = "magenta",
                                  "Alpha" = "#DC0000FF",
                                  "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                  "Omicron" = "#8491B4FF"))+
  theme_set(theme_classic(base_size=12))+
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), plot.title = element_text(face="bold", size=18))
#theme(legend.position = "none") 


###################################################3
#plot number of lineages per state
states <- ggplot(colombia, aes(division, fill = VOCI)) +
  geom_bar() +
  coord_flip()+
  theme_minimal()+
  #scale_fill_brewer(palette="Spectral")+
  scale_colour_manual("",values = c("B.1" = "#91D1C2FF",
                                    "B.1.111" = "#64CC80",
                                    "B.1.1.348" = "#80A680",
                                    "B.1.420" = "magenta",
                                    "Alpha" = "#DC0000FF",
                                    "Mu" = "#F39B7FFF",
                                    "Lambda" = "#3C5488FF",
                                    "Gamma" = "#00A087FF",
                                    "Delta" = "#4DBBD5B2",
                                    "Omicron" = "#8491B4FF"))+
  scale_fill_manual("",values = c("B.1" = "#91D1C2FF",
                                  "B.1.111" = "#64CC80",
                                  "B.1.1.348" = "#80A680",
                                  "B.1.420" = "magenta",
                                  "Alpha" = "#DC0000FF",
                                  "Mu" = "#F39B7FFF",
                                  "Lambda" = "#3C5488FF",
                                  "Gamma" = "#00A087FF",
                                  "Delta" = "#4DBBD5B2",
                                  "Omicron" = "#8491B4FF"))+
  theme(text = element_text(size=20))+
  xlab("Colombian States") +
  ylab ("No. genomes") 


############ map of colombian's states with number of available genomes


co <- getData("GADM", country = "CO", level = 1, download = TRUE)
col_depto <- fortify(co)  # make compatible to ggplot2

locat = as.vector(bbox(co))

ncmap = get_map(location = locat, source = "stamen", maptype = "toner", zoom = 6)
# ggmap(ncmap) not nice

class(col_depto$id)
totalgenomes$id <- as.character(totalgenomes$id)
col_depto <- left_join(col_depto, totalgenomes, by = "id")

col_depto2 <- col_depto %>% filter(id == "2")

mapbase <- ggplot(col_depto, aes(long, lat, group = group)) + 
  geom_polygon(aes(fill = Genomes), color = "blue") + 
  coord_equal() + 
 # geom_text(aes(label = id)) 
  scale_fill_gradient(high = "black", low = "white", guide = "colorbar") +
  geom_path(color = "grey")+
  theme(text = element_text(size=20))+
  labs(title = "Genomes per Colombian state")
  

### to see per state, you could use (id == "32") in the fill 

p_all<-ggarrange(time, states, mapbase, ncol = 3, nrow = 1)

setwd("C:/Users/Nikla/Dropbox/SARS-CoV-2_Colombia_Phylodinamics-/Figures")
ggsave("Figure2.jpg",plot=p_all,dpi=500,width=30,height=10)




##############333
table(colombia$division,colombia$VOC)

lin_share <- as.data.frame(table(colombia$division,colombia$lineage_c))
names(lin_share) <- c("state","lineage","Count")

colombia$date <- format(colombia$date, "%Y")


data <- lin_share %>%
  group_by(state, lineage) %>%
  summarise(n = sum(Count)) %>%
  mutate(Percentage = n/sum(n))


data_num_date <- data %>%
  mutate_at(vars(date), ym)


#####################################3

library(dplyr)
arrests <- USArrests 
arrests$region <- tolower(rownames(USArrests))
head(arrests)

# Retrieve the states map data and merge with crime data
states_map <- map_data("state")
arrests_map <- left_join(states_map, arrests, by = "region")

# Create the map
ggplot(arrests_map, aes(long, lat, group = group))+
  geom_polygon(aes(fill = Assault), color = "white")+
  scale_fill_viridis_c(option = "C")












##########################################################
world_map <- map_data("world")
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")

colombia <- map_data("world", region = "Colombia")

colombia.lab.data <- colombia %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

ggplot(colombia, aes(x = long, y = lat)) +
  geom_polygon(aes( group = group, fill = region))+
  geom_text(aes(label = region), data = colombia.lab.data,  size = 3, hjust = 0.5)+
  scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")



############################################


