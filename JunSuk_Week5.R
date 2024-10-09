library(geodata)
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(usdm)
library(ENMeval)
library(dismo)
library(predicts)

# 1. Download Worldclim variables
worldclim <- worldclim_global(var="bio", "res"="5", path=tempdir())


# 2. Select rasters and crop for Austalia
Aus_shp <- read_sf("C:\\Users\\iamth\\Desktop\\UMD\\Coursework\\2024 Autumn\\Spatial ecology\\Week 2\\HW\\oz_nz_aea.shp") %>%
  filter(NAME=="Australia") %>%
  st_transform(crs=4326)

Ausclim <- worldclim[[c(2:7,10,11,15,18,19)]] %>%
  crop(Aus_shp, mask=T) %>%
  rename_with(~c("bio2","bio3","bio4","bio5","bio6","bio7","bio10","bio11","bio15","bio18","bio19"))
  

# 3. Download records for Austral grass tree and clean data (remove duplicates)
grasstree <- sp_occurrence(genus = "Xanthorrhoea",
                           species = "australis",
                           download=T,
                           geo=T,
                           removeZeros = T) %>%
  filter(coordinateUncertaintyInMeters<=10000, year>1990) %>%
  distinct(lat,lon, .keep_all=T) %>%
  st_as_sf(coords=c('lon','lat'), crs=4326) %>%
  st_intersection(Aus_shp)

west <- which.min(data.frame(st_coordinates(grasstree))$X)
grasstree <- grasstree[-west,]

dups <- cellFromXY(Ausclim[[2]], st_coordinates(grasstree)) %>%
  duplicated()

grasstree <- grasstree[!dups,]


###### Q1. In general terms, how would you expect the resolution of a raster to influence the number of spatial duplicates?
### The lower the resolution, the more the number of spatial duplicates since the grid would be larger.


# 4. Extract bioclimatic variables from raster stack
grassclim <- terra::extract(Ausclim, grasstree, cell=T, xy=T, ID=F) %>%
  distinct(cell, .keep_all = TRUE) %>%
  mutate(type="occurence")


# 5. Generate 10000 background points and extract bioclimatic variables from raster stack. Combine.
set.seed(900)
background <- spatSample(Ausclim, 10000, method='random', na.rm=T, cell=T, xy=T) %>%
  mutate(type="background")

occbg <- rbind(grassclim, background)


# 6. Removed highly correlated variables
corr <- vifstep(select(occbg, -c(x,y,type)), th=10)@excluded
occbg <- select(occbg, -c(corr))


# 7. Divide data into 80% training and 20% testing through kfold
occ <- occbg %>%
  filter(type=="occurence") %>%
  select(c(x,y))

bg <- occbg %>%
  filter(type=="background") %>%
  select(c(x,y))

set.seed(900)
kfold <- get.randomkfold(filter(occbg, type=="occurence"),
                         filter(occbg, type=="background"),
                         5)

kfold <- as.numeric(unlist(kfold))
train <- subset(occbg, kfold>1) %>%
  mutate(type="train")
test <- subset(occbg, kfold==1) %>%
  mutate(type="test")


# 8. Map training, testing, and background data on one of the bioclimatic rasters
sdm_data <- rbind(filter(occbg, type=="background"),
                  train,
                  test)

sdm_data.sf <- sdm_data %>%
  select(c(x,y,type)) %>%
  st_as_sf(coords=c('x','y'), crs=4326)
  

ggplot() +
  geom_spatraster(data=Ausclim[[1]], aes(fill=bio2)) +
  scale_fill_gradient(low='white', high='green') +
  geom_sf(data=sdm_data.sf, aes(color=type), size=1) +
  scale_color_manual(values=c("background"="#FFFF00", "train"="#4876FF", "test"="#EE2C2C")) +
  theme_minimal()


# 9. Fit and predict Mahalanobis model
Ausclim <- Ausclim[[c("bio2","bio3","bio6","bio15","bio18","bio19")]]

mm <- mahal(stack(Ausclim),
            sdm_data %>%
              filter(type=="train") %>%
              dplyr::select(c(x,y)) %>%
              as.matrix())

pm <- raster::predict(stack(Ausclim),
                      mm,
                      progress='text')


# 10. Plot Mahalanobis prediction
mm.prob <- app(rast(pm), function(x, k=nlyr(Ausclim)){
  x <- 1-x
  x <- x^2
  p_value <- pchisq(x, df = k, lower.tail = FALSE)
  return(p_value)
})

ggplot() + 
  geom_spatraster(data = mm.prob, aes(fill = lyr.1)) +
  scale_fill_gradient(low='white', high='green', na.value = "transparent", 
                       name = "Prob") + 
  geom_sf(data = filter(sdm_data.sf, type=="train")) +  
  geom_sf(data = Aus_shp, fill=NA, color='black') +
  theme_minimal()


# 11. Use training data and bioclimatic rasters to fit and predict MaxEnt model. Assess variable importance and plot variable response curves. Map prediction.
filePath <- "C:\\Users\\iamth\\Desktop\\UMD\\Coursework\\2024 Autumn\\Spatial ecology\\Week 6\\maxent"

mx <- MaxEnt(x=Ausclim,
            p=sdm_data %>%
              filter(type=="train") %>%
              dplyr::select(c(x,y)),
            a=sdm_data %>%
              filter(type=="background") %>%
              dplyr::select(c(x,y)),
            path=filePath,
            args=c("jackknife","responsecurves"))

mx_pred <- raster::predict(mx, Ausclim)

ggplot() + 
  geom_spatraster(data = mx_pred, aes(fill = maxent)) +
  scale_fill_gradient(low='white', high='green', na.value = "transparent", 
                      name = "Prob") + 
  geom_sf(data = filter(sdm_data.sf, type=="train")) +  
  geom_sf(data = Aus_shp, fill=NA, color='black') +
  theme_minimal()


##### Q2. Which are the two most important variables and the least important variable associated with Austral grass tree distribution?
### A2. bio19 and bio 18 are the two most important variables. bio3 is the least important variable.

##### Q3. What can you say about bioclimatic controls on the distribution of Austral grass trees?
### A3. Seeing that the prediction of distribution is most affected by precipitation (bio19, bio18), areas with lower precipitation at both the warmest and coldest quarters will have the highest prediction accuracy.
### Although the contributions of temperature variables (bio2, bio3, bio6) are nominal, areas with higher diurnal range or isothermality will have lower prediction accuracy.

##### Q4. Comparing the predicted distributions of the two SDMs, how are they similar or different? Where do models over or under predict? What might account for these model errors?
### A4. The maxent model predictions are less dispersed and more concentrated in areas where the training data exist. However, they do both predict the possibility of distribution in the southwestern part of Australia, which could be areas they overpredict.
### Including the absent areas in the modeling process may have caused this.


# 12. Evaluate mahal and maxent models using test data
mm.prob.test <- terra::extract(mm.prob,
                               sdm_data %>%
                                 filter(type=="test") %>%
                                 dplyr::select(c(x,y)))$lyr.1
mm.prob.background <- terra::extract(mm.prob,
                               sdm_data %>%
                                 filter(type=="background") %>%
                                 dplyr::select(c(x,y)))$lyr.1

mx_pred.test <- terra::extract(mx_pred,
                               sdm_data %>%
                                 filter(type=="test") %>%
                                 dplyr::select(c(x,y)))$maxent
mx_pred.background <- terra::extract(mx_pred,
                               sdm_data %>%
                                 filter(type=="background") %>%
                                 dplyr::select(c(x,y)))$maxent

mm.eval <- pa_evaluate(p=mm.prob.test,
                       a=mm.prob.background)
mx.eval <- pa_evaluate(p=mx_pred.test,
                       a=mx_pred.background)

mm.eval
mx.eval

##### Q5. Which model performed best? How might you handle this seeming contradiction between the difference in the spatial predictions but similarity in AUC?
### A5. The AUC for the maxent model is higher, but very close to the AUC in the Mahalanobis model. I would say that the maxent model performed better because it was less spatially dispersed.
### I would go with the model that predicts less area outside the native range of the grass tree, since it wouldn't be necessary to conserve areas that the grass trees are currently not distributed.

##### Q6. How might you improve SDMs for Austral grass trees?
### A6. It seems the bioclimatic variables used for the prediction are not sufficient. I would include more bioclimatic variables and add other variables (e.g., land use change, population, distribution of competing species) that may affect the distribution of the grass trees. 