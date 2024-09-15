library(tidyverse)
library(terra)
library(tidyterra)
library(ncf)
library(gstat)
library(sf)

# 1. Generate regular spatial samples
Cwren <- rast("C:\\Users\\iamth\\Desktop\\UMD\\Coursework\\2024 Autumn\\Spatial ecology\\Week 3\\HW\\carolinaWren.tif")

Cwren_s <- spatSample(Cwren, size=500, method='regular', na.rm=T, as.df=T, xy=T,
                      ext=c(-1000000, 1000000, 3000000, 5000000))

# 2. Map sample over original map
Cwren_s <- Cwren_s %>% filter(carolinaWren>0)

plot(Cwren)
points(Cwren_s$x, Cwren_s$y, col='grey30', pch=20)

# 3. Create and plot correlogram for regular samples

# Calculate maximum distance between points
Cwren_s_sf <- st_as_sf(Cwren_s, coords=c("x","y")) %>%
  st_set_crs(st_crs(Cwren))

dist <- st_distance(Cwren_s_sf)
which(dist==max(dist), arr.ind=T)
max_dist <- dist[345,1]

# Plot correlograms with different increments
c <- correlog(x=Cwren_s$x,
              y=Cwren_s$y,
              z=Cwren_s$carolinaWren,
              increment=max_dist/15,
              resamp=1000)

c2 <- correlog(x=Cwren_s$x,
              y=Cwren_s$y,
              z=Cwren_s$carolinaWren,
              increment=max_dist/150,
              resamp=1000)

par(mfrow=c(1,2))
plot(c, main='Higher increment')
plot(c2, main='Lower increment')

### How does the correlogram change as the increment argument is changed?
# With lower increment, the correlogram captures the correlation between points within closer distance compared to the higher increment.
# On the contrary, the higher increment shows a 'smoothed' relationship between points with less fluctuations by distance

### What is your interpretation of the resulting plot? 
# The sample has higher values in points around the center, with decreasing values as the points are placed farther from the center

### Does the correlogram look different if you sample random locations instead of using a regularly spaced grid of samples?
Cwren_s2 <- spatSample(Cwren, size=500, method='random', na.rm=T, as.df=T, xy=T,
                      ext=c(-1000000, 1000000, 3000000, 5000000))
Cwren_s2 <- Cwren_s2 %>% filter(carolinaWren>0)

Cwren_s2_sf <- st_as_sf(Cwren_s2, coords=c("x","y")) %>%
  st_set_crs(st_crs(Cwren))

dist2 <- st_distance(Cwren_s2_sf)
which(dist2==max(dist2), arr.ind=T)
max_dist2 <- dist2[471,380]

c3 <- correlog(x=Cwren_s2$x,
               y=Cwren_s2$y,
               z=Cwren_s2$carolinaWren,
               increment=max_dist2/15,
               resamp=1000)

plot(c, main='Regular sample')
plot(c3, main='Random sample')
# Apparently, no. Although the intensity of correlation differs, the spatial patterns of both samples are almost identical.

# 4. Variogram
par(mfrow=c(1,1))
v <- variogram(log(carolinaWren)~1, Cwren_s2_sf, alpha=c(0,45,90,135))
plot(v)

### Is the abundance pattern isotropic or anisotropic (use variograms to support your answer)?
# Isotropy seems a reasonable assumption as the strength and pattern of the spatial correlation are broadly the same in four directions

### If you were to design a study to measure abundance of the Carolina wren across its geographic range, is
### there a distance at which you should space the sample sites to minimize spatial autocorrelation in the observations?
# The sill of the variogram in all 4 directions are close to 800000, so I would space the sample sites by 800km.