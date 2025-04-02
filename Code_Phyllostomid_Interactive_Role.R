##%######################################################%##
#                                                          #
####                Niche Structure and                 ####
####         Interactive role of Phyllostomidae         ####
#                                                          #
##%######################################################%##

pacman::p_load("here", "readr", "tidyverse", "dplyr", "spThin",                 #PHASE I
       "bipartite", "circlize", "corrr", "ggcorrplot",                          #PHASE II
       "devtools", "sf", "terra", "raster", "rgl",                              #PHASE III
       "ellipsenm",  "ntbox", "stringr", "rgdal", "ggplot2", "stars"            #PHASE IV
       )

### To run this code for other species (view the "Selected_Sp" object, line 114) 
### and change the lines with side notes for the desired species. 

here()
setwd("~Data")                                                                  # Your path


# Bat-Plant Interaction database (Taxonomy checked)
OG_NTW <- read_csv("OG_NTW.csv")

######################################### PHASE I: Network extraction (1-9) ----
# 1) Rename and reorganize data
OG_NTW <- OG_NTW %>% 
  rename(Bat = `Bat specie`,Plant = `Plant specie`) %>% 
  dplyr::select(!c(Interaction,6:8, 11:18)) %>% 
  filter(Interaction_type == "Frugivory") %>% 
  dplyr::select(!c(Interaction_type))
Sites <- as.factor(OG_NTW$ID_sites)

# 2) Filter only networks with a minimum of 6 nodes
Optimal_NTWs <- OG_NTW %>% 
  group_by(ID_sites) %>% 
  filter(n() > 6)

# 3)  Identify networks formed by less of 6 species of bats and plants
# 6 bats
Unvalid_NTWs_Bats <- Optimal_NTWs %>%
  group_by(ID_sites) %>%
  summarize(n_sp_bat = n_distinct(Bat)) %>% 
  filter(n_sp_bat < 6) %>% 
  dplyr::select(ID_sites) 
# 6 plants
Unvalid_NTWs_Plants <- Optimal_NTWs %>%
  group_by(ID_sites) %>%
  summarize(n_sp_plant = n_distinct(Plant)) %>% 
  filter(n_sp_plant < 6) %>% 
  dplyr::select(ID_sites)
# bind into a single df
Unvalid_NTWs <- bind_rows(Unvalid_NTWs_Bats, Unvalid_NTWs_Plants)
Unvalid_NTWs <- distinct(Unvalid_NTWs)

# 4) Selecting only useful networks (6 interactors and 6 especies per network)
Valid_NTWs <- anti_join(Optimal_NTWs, Unvalid_NTWs, by="ID_sites")

# 5) Filtering networks by geographic threshold 
set.seed(111)
Geo_Unique <- unique(Valid_NTWs$ID_sites)
Geo_ID <- rep("A", length(Geo_Unique))
Geo_Coords <- Valid_NTWs %>% 
  ungroup() %>% 
  dplyr::select(Latitude, Longitude) %>% 
  distinct()
Geo_NTWs <- bind_cols(Geo_ID, Geo_Coords)
Geo_NTWs <- rename(Geo_NTWs, Sp = `...1`)

# Aplying the filter
Geo_Filter <- thin(
  Geo_NTWs,
  lat.col = "Latitude",
  long.col = "Longitude",
  spec.col = "Sp",
  thin.par = 5,
  reps = 1,
  locs.thinned.list.return = T,
  write.files = F,
  verbose = F
)
summaryThin(Geo_Filter, T)
FinalCoords <- as.data.frame(Geo_Filter[[1]])
FinalLon <- FinalCoords %>% 
  dplyr::select(Longitude)

# 6) Selecting the final Networks
Final_NTWs <- merge(Valid_NTWs, FinalLon, by="Longitude")
Final_NTWs <- Final_NTWs %>% 
  relocate(Longitude, .after = Latitude) %>% 
  as_tibble
NTWs_List <- Final_NTWs %>% 
  distinct(ID_sites)

# Extracting all networks into individual DF and exporting to csv's (Optional) 
#for(i in unique(Final_NTWs$ID_sites)) {
#  nam <- paste("Network", i, sep = "_")
#  assign(nam, Final_NTWs[Final_NTWs$ID_sites==i,])
#}
#by(Final_NTWs, Final_NTWs$ID_sites, FUN=function(i) write_csv(i, paste0(i$ID_sites[1], ".csv")))

# 7) Obtaining candidate species for the study (participants in 6 or more networks)
Selected_Sp <- Final_NTWs %>%
  group_by(ID_sites, Bat) %>%         #Cols used for grouping
# Species found in more than two times per network
  mutate(Interactions_Per_Newtwork = length(Bat),
         withinGroup = ifelse(Interactions_Per_Newtwork > 1,T,F)) %>%  
# Grouping by species
  group_by(Bat) %>%
# Number of networks every species is interacting in
  mutate(Total_NTWs = ifelse(n_distinct(ID_sites) == 1, 0, n_distinct(ID_sites)),
         Multiple_NTWs            = ifelse(Total_NTWs  > 0,T,F)) %>% 
# Merging results
  dplyr::select(ID_sites,Bat, withinGroup, 
                Interactions_Per_Newtwork, Multiple_NTWs, Total_NTWs) %>% 
  ungroup()

# 8) List of candidate species
Selected_Sp <- Selected_Sp %>% 
  filter(Total_NTWs > 5) %>% 
  distinct(Bat)

######################################### PHASE II: Network Analysis (10-20) ----
# 10) Extracting networks for a specific species                                # Spp change
NTWs_by_Spp <- Final_NTWs %>% 
  group_by(ID_sites) %>% 
  mutate(col_3 = as.integer("Artibeus fimbriatus" %in% Bat)) %>%                  
  rename(Artibeus_fimbriatus = col_3) %>%                       
  filter(Artibeus_fimbriatus == "1") %>% 
  dplyr::select(!Artibeus_fimbriatus)
Redes_Murci_porSp <- NTWs_by_Spp %>% distinct(ID_sites)
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
List_by_Sp <- split_tibble(NTWs_by_Spp, 'ID_sites') %>% 
  lapply("[",1:3) # Column selection

# 11) Matrix transformation
Matrix_List <- list()
for (i in 1:length(List_by_Sp)) {
  Matrix <- data.frame(frame2webs(List_by_Sp[[i]],
                                  varnames = c("Bat","Plant","ID_sites"), 
                                  type.out="list"))
  Matrix_List[[i]] <- Matrix
}

# 12) Bipartite network analysis: "species level" 
Matrix_Stats <- list()
for (i in 1:length(Matrix_List)) {
  Stats <- specieslevel(Matrix_List[[i]], level = "lower")
  Matrix_Stats[[i]] <- Stats
}

# 13) Metrics selection: N. degree, closeness, betwenness
Stats_Selected <- list()
for (i in 1:length(Matrix_Stats)) {
  Selected <- Matrix_Stats[[i]] %>% 
    dplyr::select(c(normalised.degree, 
             closeness, betweenness))
  Stats_Selected[[i]] <- Selected
}                        

# 14) Matrix transposition
Species_Full <- list()
for (i in 1:length(Stats_Selected)) {
  Full <- data.frame(t(Stats_Selected[[i]]))
 Species_Full[[i]] <- Full
}

# 15) Metrics for a specific species                                            # Spp change
Species_Selected <- list()
for (i in 1:length(Species_Full)) {
  Selected <- Species_Full[[i]] %>% 
    dplyr::select(Artibeus.fimbriatus)                                        
  Species_Selected[[i]] <- Selected
}

# 16) Returning matrix to original position
Centralidad_Murci <- list()
for (i in 1:length(Species_Selected)) {
  Centralidad <- data.frame(t(Species_Selected[[i]]))
  Centralidad_Murci[[i]] <- Centralidad
}

# Extract networks from DF to individual DF and exporting to csv (Optional)
#for(i in unique(Final_NTWs$ID_sites)) {
#  nam <- paste("Network", i, sep = "_")
#  assign(nam, Final_NTWs[Final_NTWs$ID_sites==i,])
#}

# 17) Merging the information of all the Networks in a DF
Final_Stats <- list_rbind(Centralidad_Murci)

# 18) PCA to reduce the number of dimensions
#Corr matrix
corr_matrix <- cor(Final_Stats)

#PCA
pca <- prcomp(Final_Stats, scale = T)
summary(pca) # Proporción de la varianza explicada (PC1 - 43%)
ResultsPCA <- data.frame(pca$x) #Contribución por elemento
PC1 <- dplyr::select(ResultsPCA, -PC2, -PC3)

# 19) Final Centrality value for the final model
Centrality <- PC1 + abs(min(PC1)) + 0.00001

# 20 NTWs as spatial points                                                     # Spp change
Points_by_NTWS <- split_tibble(NTWs_by_Spp, 'ID_sites') %>% 
  lapply("[",1:5) 
Points_NTWs <- list_rbind(Points_by_NTWS)
Points_NTWs <- Points_NTWs %>% distinct(ID_sites, .keep_all = T) %>% 
  dplyr::select(Longitude,Latitude) %>% 
  mutate(Sp = "Artibeus_fimbriatus", .before = Longitude) %>% ungroup()

########### PHASE III: Climatic vars and Species Occurrence records (21-22) ----
# 21) Spatial Data
# Worldclim vars
vars <- raster::stack(list.files("~Data",      # Your path
                                 "Artibeus_fimbriatus.+.tif"))                  # Spp change

# Distribution polygon
Shp <- st_read("Artibeus_fimbriatus_Shp.shp")                                   # Spp change

# 22) Occ: Test & Train data (Random Data splitting was done only one time) 
Points_GBIF <-  read_csv("Artibeus_fimbriatus_Occ.csv")                         # Spp change
#Points_TT <- ellipsenm::split_data(Points_GBIF, 
#                                   method = "random", 
#                                   longitude = "lon", latitude = "Lat", 
#                                   raster_layer = NULL, train_proportion = 0.75, 
#                                   save = F)
#write.csv(Points_TT[[2]], "Artibeus_fimbriatus_Occ_Train.csv")
#write.csv(Points_TT[[3]], "Artibeus_fimbriatus_Occ_Test.csv")

M_train <- read_csv("Artibeus_fimbriatus_Occ_Train.csv")
M_test <- read_csv("Artibeus_fimbriatus_Occ_Test.csv")

############################# PHASE IV: Ecological Niche Modelling (23-28) -----
# 23) Climate information extraction
# Train
M_EnvTrain <- raster::extract(vars,M_train[,c("Lon","Lat")],df=TRUE)
M_EnvTrain <- M_EnvTrain[,-1]
M_EnvTrain <- na.omit(M_EnvTrain)
head(M_EnvTrain)

# Test
M_EnvTest <- raster::extract(vars,M_test[,c("Lon","Lat")],df=TRUE)
M_EnvTest <- M_EnvTest[,-1]
M_EnvTest <- na.omit(M_EnvTest)
head(M_EnvTest)

# 24) Non correlated variables
Env_matrix <- cor(M_EnvTrain)
ggcorrplot(Env_matrix)
env_varsL <- ntbox::correlation_finder(cor(M_EnvTrain,method = "spearman"),
                                       threshold = 0.8,
                                       verbose = F)
env_vars <- env_varsL$descriptors

# 25) Model Calibration
nvarstest <- c(3,4,5) # Number of possible variables
level <- 0.95 # Proportion of train occ used to adjust the volume
env_bg <- ntbox::sample_envbg(vars,1000) # Background points for the partial-ROC
omr_criteria <- 0.05 # Omision rate
proc <- TRUE # Partial-ROC value returned
set.seed(111)

# 26) Model competition
e_selct <- ntbox::ellipsoid_selection(env_train = M_EnvTrain,
                                      env_test = M_EnvTest,
                                      env_vars = env_vars,
                                      level = level,
                                      nvarstest = nvarstest,
                                      env_bg = env_bg,
                                      omr_criteria= omr_criteria,
                                      proc = proc)
# Highest auc and omr
e_selct <- e_selct %>% 
  mutate(Final_Rank = 
           ((rank_by_omr_train_test + rank_omr_aucratio)*rank_by_omr_train_test))
e_selct  <- e_selct %>% arrange(Final_Rank)  

# 27) Final Model
bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")[[1]]

# E espace
best_mod <- ntbox::cov_center(M_EnvTrain[,bestvarcomb],
                              mve = T,
                              level = 0.99,
                              vars = 1:length(bestvarcomb))
# G Espace
Model_Final <- ntbox::ellipsoidfit(vars[[bestvarcomb]],
                             centroid = best_mod$centroid,
                             covar = best_mod$covariance,
                             level = 0.99,size = 3)
# NOTE: DO NOT CLOSE THE RGL PLOT UNTIL SAVING (LINE 351)

Suitability_Rast <- Model_Final$suitRaster

# 28) Extracting Suitability values for each network
Suitability_Data <- raster::extract(Model_Final$suitRaster,
                                    Points_NTWs[,c("Longitude","Latitude")],
                                    df=TRUE)
Suitability_NTWs <- as_tibble(Suitability_Data$suitability)

####################################### PHASE V: Statistical Models (29-33) ----
# 29) Merging suitability and centrality in a single DF
FinalFrame <- as_tibble(bind_cols(Suitability_NTWs, Centrality))
FinalFrame <- FinalFrame %>% 
  rename(Suitability = value, Interactive_Role = PC1) 
FinalFrame <- na.omit(FinalFrame)
attach(FinalFrame)
plot(Suitability, Interactive_Role)
#par(mfrow = c(2,2), mar = c(4,4,4,4))

# 30) LM                                                                        # Spp change
LM_Artibeus_fimbriatus <- lm(Interactive_Role ~ Suitability)
summary(LM_Artibeus_fimbriatus)
# plot(LM_Artibeus_fimbriatus)
LM_Murci <- ggplot(data = FinalFrame, aes(x = Suitability, y = Interactive_Role)) +
  geom_smooth(method = "lm") +
  geom_point(size = 5) +
  theme_classic()
LM_Murci_Plot <- LM_Murci + ggtitle("LM: Artibeus_fimbriatus")

# 31) GLM                                                                       # Spp change
GLM_Artibeus_fimbriatus <- MASS::glm.nb(Interactive_Role ~ Suitability)
summary(GLM_Artibeus_fimbriatus)
#plot(GLM_Artibeus_fimbriatus)
#R-squared
deviance <- summary(GLM_Artibeus_fimbriatus)$deviance
null_deviance <- summary(GLM_Artibeus_fimbriatus)$null.deviance
rsquared <- 1 - (deviance / null_deviance)
rsquared

# 32) Quadratic Regression model                                                # Spp change
QM_Artibeus_fimbriatus <- lm(Interactive_Role ~ Suitability + I(Suitability^2))
summary(QM_Artibeus_fimbriatus)
# plot(QM_Artibeus_fimbriatus)
QM_Murci <- ggplot(data = FinalFrame, aes(x = Suitability, y = Interactive_Role)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_point(size = 5) + 
  theme_classic()
QM_Murci_Plot <- QM_Murci + ggtitle("QM: Artibeus_fimbriatus")

# 33) Power regression model                                                    # Spp change
PM_Artibeus_fimbriatus <- lm(log(Interactive_Role) ~ log(Suitability))
summary(PM_Artibeus_fimbriatus)
#plot(PM_Artibeus_fimbriatus)
PM_Murci <- ggplot(data = FinalFrame, aes(x = Suitability, y = Interactive_Role)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_point(size = 5) + 
  theme_classic()
PM_Murci_Plot <- PM_Murci + ggtitle("PM: Artibeus_fimbriatus")

######################### Phase VI: Exporting Figures and Results (34 - 40) ----
# NOTE: 34, 35 and 40 need to be exported only once.

# 34) List of selected networks
write.csv(Final_NTWs, "Outputs/NTWs_Final_List.csv")

# 35) Fig 1: Networks in America 
# America Shapefile
Shp_America <- st_read("Shp_America.shp")
America <- ggplot() +
 geom_sf(data=Shp_America,
         mapping = aes(fill = Species),
         show.legend = FALSE,
         colour = alpha("black", 0.8),
         fill = NA) +
 labs(x="Longitude",
      y="Latitude") +
   theme_classic()
ggsave("Outputs/Fig_America.pdf", America, width = 8, height = 6, units = "in")

# Network Points
Map_NTWs <- Final_NTWs %>% dplyr::select(Latitude, Longitude) %>% 
  distinct(Latitude, .keep_all = T)

# America with Networks
Fig1 <- ggplot() + 
 geom_sf(data=Shp_America, 
         show.legend = FALSE, 
         colour = alpha("black", 0.5), 
         fill = "lightgray") +
 geom_point(Map_NTWs, mapping = aes(x= Longitude, y=Latitude), color = alpha("blue",0.5), 
            fill= "cornflowerblue",  shape=21, size = 4) +
 labs(x="Longitude",
      y="Latitude") +
 theme(plot.title = element_text(hjust = 2)) + theme_classic()
ggsave("Outputs/Fig1.pdf", Fig1, width = 8, height = 6, units = "in")

# 36) Elipsoid. NOTE: Move the 3D plot for better visualization                 # Spp change
rgl.snapshot("Outputs/Artibeus_fimbriatus_Ellipsoid.png", fmt = "png", top = TRUE)

# 37) Fig 2 (Suitability Map with network)                                      # Spp change
# Suitability Raster
Suitability_Rast
writeRaster(Suitability_Rast, "Outputs/Artibeus_fimbriatus_Rast.tiff")

# NTW points
Points_NTWs_Full <- bind_cols(Points_NTWs, Interactive_Role = FinalFrame$Interactive_Role)

# Final Suitability Map
Suitability_Rast_Stars <- st_as_stars(Suitability_Rast)
Final_Map <- ggplot() + 
  geom_stars(data = Suitability_Rast_Stars) + 
  coord_equal() +
  scale_fill_gradient(low = "white", high = "black", na.value = "transparent")+
  geom_sf(data=Shp_America, 
          mapping = aes(fill = Species), 
          show.legend = FALSE, 
          colour = alpha("black", 0.7), 
          fill = NA) +
  geom_point(Points_NTWs_Full,
             mapping = aes(x= Longitude, y=Latitude, color = Interactive_Role), 
             shape=19, size = 5) +
  scale_colour_gradient(low = "orange", high = "red") +
  labs(x="Longitude",
       y="Latitude",
       title=expression(paste("Mapa de idoneidad y redes de ", italic("Artibeus fimbriatus")))) +
  theme(plot.title = element_text(hjust = 0.5)) + theme_classic()
ggsave("Outputs/Artibeus_fimbriatus_SuitMap.pdf", Final_Map, width = 8, height = 6, units = "in")

# 38) Results Interactive role suitability
write.csv(FinalFrame, "Outputs/Artibeus_fimbriatus_Results.csv")

# 39) Statistic models regression line
LM_Murci_Plot
ggsave("Outputs/Artibeus_fimbriatus_LM.pdf", LM_Murci_Plot, width = 8, height = 6, units = "in")
PM_Murci_Plot
ggsave("Outputs/Artibeus_fimbriatus_PM.pdf", PM_Murci_Plot, width = 8, height = 6, units = "in")
QM_Murci_Plot
ggsave("Outputs/Artibeus_fimbriatus_QM.pdf", QM_Murci_Plot, width = 8, height = 6, units = "in")

# 40) Networks plots (Chord Diagram)
circos.par(start.degree = 90, clock.wise = T)
pdf("Outputs/Plot_All_NTWs.pdf", width = 8, height = 8)
for(i in unique(Final_NTWs$ID_sites)) {
  rowI <- which(Final_NTWs$ID_sites == i)
  Network <- table(Final_NTWs[rowI, c(1,2)])
  Full_NTWs <- data.matrix(Network)
  chordDiagram(Full_NTWs,
               grid.col = "lightblue", col = "grey",
               annotationTrack = c("name","grid","axes"), big.gap = 30)
  title(i)
}
dev.off()
circos.clear()

##%######################################################%##
#                                                          #
####                      THE END :D                    ####
#                                                          #
##%######################################################%##