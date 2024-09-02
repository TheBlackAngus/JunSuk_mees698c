library(tidyverse)
library(tidygeocoder)
library(daymetr)
library(lubridate)

# 1. Create dataframe with three locations
places <- data.frame(
 Name = c("High", "Mid", "Low"),
 Address = c("Fairbanks, Alaska", "Frostburg, Maryland", "Galveston, Texas")
 )

# 2. Geocode coordinates into dataframe
places <- places %>%
  geocode(address = Address, method = "arcgis")

# 3. Download daymet data
get_daymet <- function(i){
  temp_lat <- places[i, ] %>% pull(lat)
  temp_lon <- places[i, ] %>% pull(long)
  temp_name <- places[i, ] %>% pull(Name)
  temp_address <- places[i, ] %>% pull(Address)
  temp_daymet <- download_daymet(
    lat = temp_lat,
    lon = temp_lon,
    start=1980,
    end=2020
  ) %>%
    .$data %>%
    mutate(Name = temp_name,
           Address = temp_address) %>%
    group_by(year, yday) %>%
    transmute(Name = Name,
              Address = Address,
              year= year,
              yday = yday,
              date = as.Date(paste(year, yday, sep="-"), "%Y-%j"),
              temp = (tmax..deg.c.+tmin..deg.c.)/2)
  return(temp_daymet)
}

places_daymet <- lapply(1:nrow(places), get_daymet) %>%
  bind_rows()

# 4. Aggregate for summer season each year
places_daymet_summer <- places_daymet %>%
  filter(yday %in% c(170:260)) %>%
  group_by(Name, year) %>%
  transmute(Name = Name,
            Address = Address,
            year = year,
            summmer_temp = mean(temp)) %>%
  distinct()

# 5. Plot summer temperature
plot <- places_daymet_summer %>%
  ggplot(aes(x=year, y=summmer_temp, colour=Address)) +
  scale_x_continuous(n.breaks=8) +
  scale_y_continuous(n.breaks=5) +
  geom_line(linewidth=2) +
  geom_smooth(linetype='dashed') +
  labs(title='Average summer temperature of 3 locations by year',
       x='Year',
       y='Temperature (C)') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.key.size = unit(2,'cm'),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        title = element_text(size=24))

plot

ggsave(plot, file="Summer temperatures by year.png",
         width=36, height = 24, units = "cm")
