
# load library tidyverse

library(tidyverse)
library(ggplot2)
library(scatterpie)

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, header = F)

Continent <- read.delim(file = "Continent.Coord.xls", row.names = 1, as.is = T)

MAF <- read.delim(file = "MAF.100G.xls", row.names = 1, as.is = T)


Continent$Risk <- colMeans(MAF[1:30, ])[rownames(Continent)]

Continent$noRisk <- 1 - Continent$Risk

Continent$Region <- rownames(Continent)

# create data for world coordinates using 
# map_data() function
world_coordinates <- map_data("world")
  
# create world map using ggplot() function
pp <- ggplot() +
  
# geom_map() function takes world coordinates 
# as input to plot world map
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region),
    color = "white", fill = Colors["Yellow", ], size = 0.2
  ) +

  geom_scatterpie(aes(x=long, y=lat, group=Region), data=Continent,
                           cols=c("Risk", "noRisk"), pie_scale = 1.5) + 

  scale_fill_manual(values = Colors[c("Red", "Blue"), ]) +

  coord_equal() +

  theme_void() +

  theme(legend.title = element_blank())

ggsave("WorldMap.haplo.pdf", pp, width = 7, height = 4)


Continent$Risk <- colMeans(MAF["rs893403", ])[rownames(Continent)]

Continent$noRisk <- 1 - Continent$Risk

Continent$Region <- rownames(Continent)

# create data for world coordinates using 
# map_data() function
world_coordinates <- map_data("world")

# create world map using ggplot() function
pp <- ggplot() +

# geom_map() function takes world coordinates 
# as input to plot world map
  geom_map(
    data = world_coordinates, map = world_coordinates,
    aes(long, lat, map_id = region),
    color = "white", fill = Colors["Yellow", ], size = 0.2
  ) +

  geom_scatterpie(aes(x=long, y=lat, group=Region), data=Continent,
                           cols=c("Risk", "noRisk"), pie_scale = 1.5) +

  scale_fill_manual(values = Colors[c("Red", "Blue"), ]) +

  coord_equal() +

  theme_void() +

  theme(legend.title = element_blank())

ggsave("WorldMap.NEJ.pdf", pp, width = 7, height = 4)

