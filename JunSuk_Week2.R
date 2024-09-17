library(tidyverse)
library(geodata)
library(terra)
library(tidyterra)
library(sf)
library(stars)
library(ggpubr)

# 1. Download Worldclim
worldclim <- worldclim_global(var="bio", "res"="5", path=tempdir())

# 2. Raster stack with four variables and crop to Australia
Aus_shp <- read_sf('https://github.com/fitzLab-AL/spatialEcology2024_classShare/blob/main/codingAssignments/ca-2/oz_nz_aea.shp') %>%
  filter(NAME=="Australia") %>%
  st_transform(crs=4326)

#####code above did work for me, not sure why

Ausclim <- worldclim[[c(10,11,18,19)]] %>%
  crop(Aus_shp, mask=T)

# 3. Download records for Xanthorrhoea australis
grasstree <- sp_occurrence(genus="Xanthorrhoea", species="australis",
                           ext=ext(Ausclim), geo=TRUE, path=tempdir()) %>%
  st_as_sf(coords=c('lon','lat'), crs=4326) %>%
  st_intersection(Aus_shp) %>%
  select(acceptedScientificName, institutionCode, year)

st_write(grasstree, file="grasstree in Australia.shp")

# 4. Map Xanthorrhoea australis

#####Nice plot, but in the data are in a 
##### strange / wrong projection. See assignment.
#####Also, there is an erroneous record in the 
#####center of Australia.

grasstree_bio10 <- ggplot() +
  geom_spatraster(data=project(Ausclim[[1]], "EPSG:5070")) +
  scale_fill_viridis_b() +
  geom_sf(data=st_transform(grasstree, crs=5070), aes(color=year)) +
  scale_color_gradientn(colours=c("pink","red","purple")) +
  geom_sf(data=st_transform(Aus_shp, crs=5070), fill=NA, color='black')

ggsave(grasstree_bio10, file="Grasstree by year on bio10.tif",
       width=24, height = 36, units = "cm")

# 5. Extract variables and compare
grasstree_clim <- extract(Ausclim, grasstree, ID=FALSE, cells=FALSE, xy=T) %>%
  transmute(bio10_grasstree=wc2.1_5m_bio_10,
            bio11_grasstree=wc2.1_5m_bio_11,
            bio18_grasstree=wc2.1_5m_bio_18,
            bio19_grasstree=wc2.1_5m_bio_19,
            x=x, y=y)

#####Check the assignment - we needed just a sample
##### of the data, not the entire dataset for Australia.

Ausclim_df <- as.data.frame(Ausclim, xy=T) %>%
  transmute(bio10_aus=wc2.1_5m_bio_10,
            bio11_aus=wc2.1_5m_bio_11,
            bio18_aus=wc2.1_5m_bio_18,
            bio19_aus=wc2.1_5m_bio_19,
            x=x, y=y)

Aus_grasstree_clim <- Ausclim_df %>%
  left_join(grasstree_clim, by=c("x"="x", "y"="y")) %>%
  select(-c(x,y)) %>%
  gather(key="Variable", value="value", 1:8)

bio10_plot <- Aus_grasstree_clim %>%
  filter(Variable %in% c("bio10_aus", "bio10_grasstree")) %>%
  ggplot(aes(x=Variable, y=value, fill=Variable)) +
  geom_jitter(color='grey50', size=0.4, alpha=0.8) +
  geom_violin() +
  scale_fill_discrete() +
  theme_light() +
  theme(legend.position = 'none',
        plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        title = element_blank())

bio11_plot <- Aus_grasstree_clim %>%
  filter(Variable %in% c("bio11_aus", "bio11_grasstree")) %>%
  ggplot(aes(x=Variable, y=value, fill=Variable)) +
  geom_jitter(color='grey50', size=0.4, alpha=0.8) +
  geom_violin() +
  scale_fill_discrete() +
  theme_light() +
  theme(legend.position = 'none',
        plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        title = element_blank())

bio18_plot <- Aus_grasstree_clim %>%
  filter(Variable %in% c("bio18_aus", "bio18_grasstree")) %>%
  ggplot(aes(x=Variable, y=value, fill=Variable)) +
  geom_jitter(color='grey50', size=0.4, alpha=0.8) +
  geom_violin() +
  scale_fill_discrete() +
  theme_light() +
  theme(legend.position = 'none',
        plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        title = element_blank())

bio19_plot <- Aus_grasstree_clim %>%
  filter(Variable %in% c("bio19_aus", "bio19_grasstree")) %>%
  ggplot(aes(x=Variable, y=value, fill=Variable)) +
  geom_jitter(color='grey50', size=0.4, alpha=0.8) +
  geom_violin() +
  scale_fill_discrete() +
  theme_light() +
  theme(legend.position = 'none',
        plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        title = element_blank())

ggarrange(bio10_plot, bio11_plot, bio18_plot, bio19_plot,
          labels=c("Mean warmest temp", "Mean coldest temp", "Warmest prec", "Coldest prec"),
          ncol = 2, nrow = 2)

## Compared to the entire Australia, X.aus prefers colder temperatures throughout the year and higher precipitations during cold seasons.

# 6. Create raster of the number of Xanthorrhoea australis
v <- grasstree %>%
  mutate(value=1) %>%
  vect()
r <- rast(v, res=res(Ausclim))

grasstree_rast <- rasterize(v, r, "value", fun="sum", background=0)