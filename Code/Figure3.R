library(terra)
library(tidyterra)
library(ggplot2)
library(rnaturalearth)
library(sf) 
library(ggspatial)
library(rcolors)
library(ggpubr)

#global map projection of Robinson 
crs_robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

#Get land polygons at coarse scale
global_land <- ne_countries(scale = 110, returnclass = "sf")

#Then transform projection to Robinson (same as your raster)
crs_robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
global_land <- st_transform(global_land, crs_robin)
global_land <- st_make_valid(global_land)

#(Optionally) merge all polygons into one for simpler plotting
global_land <- st_union(global_land)
global_land_vect <- vect(global_land)

#Define parallels and meridians
lat_breaks_dashed <- c(-30, 0, 30, 60)
lat_breaks_solid  <- c(-61.5, 88)   # top and bottom edges

lon_breaks_dashed <- c(-120, -60, 0, 60, 120)
lon_breaks_solid  <- c(-180, 180) # left and right edges

#Labels
lat_labels <- c("30°S", "0°", "30°N", "60°N")
lon_labels <- c("120°W", "60°W", "0°", "60°E", "120°E")

lat_lab_df <- data.frame(lon = -180, lat = lat_breaks_dashed, lab = lat_labels)
lon_lab_df <- data.frame(lon = lon_breaks_dashed, lat = -61.5, lab = lon_labels)  # bottom edge in degrees

#projection map__________________________________________________________________________________
rfn_top <- "CUE_Top_mean_prediction.tif"
r <- rast(rfn_top)
r <- terra::project(r, crs_robin)

Fig3a <- ggplot()+
  # horizontal dashed latitude lines
  annotation_spatial_hline(
    intercept = lat_breaks_dashed, crs = 4326,
    linetype = "dashed", linewidth = 0.3, colour = "grey40"
  ) +
  
  # vertical dashed longitude lines
  annotation_spatial_vline(
    intercept = lon_breaks_dashed, crs = 4326,
    linetype = "dashed", linewidth = 0.3, colour = "grey40"
  ) +
  
  # solid boundary lines
  annotation_spatial_vline(
    intercept = lon_breaks_solid, crs = 4326,
    linetype = "solid", linewidth = 0.4, colour = "black"
  ) +
  annotation_spatial_hline(
    intercept = lat_breaks_solid, crs = 4326,
    linetype = "solid", linewidth = 0.4, colour = "black"
  ) +

  # latitude labels
  geom_spatial_text(
    data = lat_lab_df,
    aes(x = lon, y = lat, label = lab),
    crs = 4326, size = 3, hjust=1
  ) +

  # longitude labels (bottom edge)
  geom_spatial_text(
    data = lon_lab_df,
    aes(x = lon, y = lat, label = lab),
    crs = 4326, size = 3,vjust=-1
  ) +
  geom_spatraster(data=r, maxcell = 1e7)+
  geom_spatvector(data = global_land_vect, fill=NA, size=0.2)+
  coord_sf(crs = crs_robin)+
  scale_y_continuous(expand = c(0,0), limits = c(-65*10^5, 86*10^5))+
  scale_x_continuous(expand = c(0,0), limits = c(-180*10^5, 180*10^5))+
  labs(title = "")+
  scale_fill_stepsn(
    colors = get_color("MPL_YlGn",n=6),
    na.value = NA,
    name="CUE",
    limits = c(0,0.6),
    breaks = seq(0,0.6,0.1),
    guide = guide_colorbar(direction = "horizontal",
                           nrow=1,
                           title.position="left",
                           barwidth=20
    ))+
  theme_void()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8,family = "Arial", color="black" ),
        legend.title = element_text(size=10, family = "Arial", color="black"),
        plot.background = element_blank()
)
Fig3a

rfn_sub <- "CUE_Sub_mean_prediction.tif"
r <- rast(rfn_sub)
r <- terra::project(r, crs_robin)

Fig3c <- ggplot()+
  # horizontal dashed latitude lines
  annotation_spatial_hline(
    intercept = lat_breaks_dashed, crs = 4326,
    linetype = "dashed", linewidth = 0.3, colour = "grey40"
  ) +
  
  # vertical dashed longitude lines
  annotation_spatial_vline(
    intercept = lon_breaks_dashed, crs = 4326,
    linetype = "dashed", linewidth = 0.3, colour = "grey40"
  ) +
  
  # solid boundary lines
  annotation_spatial_vline(
    intercept = lon_breaks_solid, crs = 4326,
    linetype = "solid", linewidth = 0.4, colour = "black"
  ) +
  annotation_spatial_hline(
    intercept = lat_breaks_solid, crs = 4326,
    linetype = "solid", linewidth = 0.4, colour = "black"
  ) +
  # latitude labels
  geom_spatial_text(
    data = lat_lab_df,
    aes(x = lon, y = lat, label = lab),
    crs = 4326, size = 3, hjust=1
  ) +
  
  # longitude labels (bottom edge)
  geom_spatial_text(
    data = lon_lab_df,
    aes(x = lon, y = lat, label = lab),
    crs = 4326, size = 3,vjust=-1
  ) +
  geom_spatraster(data=r, maxcell = 1e7)+
  geom_spatvector(data = global_land_vect, fill=NA, size=0.2)+
  coord_sf(crs = crs_robin)+
  scale_y_continuous(expand = c(0,0), limits = c(-65*10^5, 86*10^5))+
  scale_x_continuous(expand = c(0,0), limits = c(-180*10^5, 180*10^5))+
  labs(title = "")+
  scale_fill_stepsn(
    colors = get_color("MPL_YlGn",,n=6),
    na.value = NA,
    name="CUE",
    limits = c(0,0.6),
    breaks = seq(0,0.6,0.1),
    guide = guide_colorbar(direction = "horizontal",
                           nrow=1,
                           title.position="left",
                           barwidth=20
    ))+
  theme_void()+
  theme(legend.position = "bottom",
        plot.background = element_blank())
Fig3c

df_top <- read.csv("CUE_Top_mean_prediction.csv")

top_summary <- df_top %>%
  group_by(lat = round(y, 0.5)) %>%
  summarise(mean_cue = mean(mean_pred, na.rm=T),
            mean_upper = mean(upper_95, na.rm=T),
            mean_lower = mean(lower_95, na.rm=T))
head(top_summary)

Fig3b <- ggplot(top_summary, aes(x=lat, y=mean_cue))+
  geom_ribbon(aes(ymin = mean_lower, ymax = mean_upper), fill="gray", alpha=0.2)+
  geom_smooth(se=F,color="black",linetype="dashed", size=0.5)+
  geom_line(color="darkgreen", size=0.5, alpha=0.6)+
  theme_bw()+
  scale_x_continuous(expand = c(0,0), breaks = c(-30,0,30,60), labels = c("30˚S","0","30˚N","60˚N"))+
  scale_y_continuous(expand = c(0,0),limits = c(0.2,0.5), labels = seq(0.2,0.5,0.1))+
  coord_flip()+
  labs(x="",y="CUE")+
  theme(axis.ticks.length.y = unit(0.1, "cm"),
        plot.background = element_blank())
Fig3b  

df_sub <- read.csv("CUE_Sub_mean_prediction.csv")

sub_summary <- df_sub %>%
  group_by(lat = round(y, 0.5)) %>%
  summarise(mean_cue = mean(mean_pred, na.rm=T),
            mean_upper = mean(upper_95, na.rm=T),
            mean_lower = mean(lower_95, na.rm=T))
head(sub_summary)

Fig3d <- ggplot(sub_summary, aes(x=lat, y=mean_cue))+
  geom_ribbon(aes(ymin = mean_lower, ymax = mean_upper), fill="gray", alpha=0.2)+
  geom_smooth(se=F,color="black",linetype="dashed", size=0.5)+
  geom_line(color="forestgreen", size=0.5)+
  theme_bw()+
  scale_x_continuous(expand = c(0,0), breaks = c(-30,0,30,60), labels = c("30˚S","0","30˚N","60˚N"))+
  scale_y_continuous(expand = c(0,0),limits = c(0.2,0.5), labels = seq(0.2,0.5,0.1))+
  coord_flip()+
  labs(x="",y="CUE")+
  theme(axis.ticks.length.y = unit(0.1, "cm"),
        plot.background = element_blank())
Fig3d 

library(patchwork)
Fig3 <- Fig3a + Fig3b + Fig3c + Fig3d + plot_layout(ncol = 2, nrow = 2)
Fig3

ggsave("figures/figures/Fig3.tiff", width = 12, height = 8, dpi = 300)
#plot uncertainty__________________________________________________________________________________________
u_top <- "CUE_Top_uncertainty.tif"
u <- rast(u_top)
u <- terra::project(u, crs_robin)

FigS6a <- ggplot()+
  geom_spatraster(data=u, maxcell = 1e7)+
  geom_spatvector(data = global_land_vect, fill=NA, size=0.2)+
  coord_sf(crs = crs_robin)+
  scale_y_continuous(expand = c(0,0), limits = c(-65*10^5, 86*10^5))+
  scale_x_continuous(expand = c(0,0), limits = c(-180*10^5, 180*10^5))+
  labs(title = "Uncertainty of CUE in Topsoil")+
  scale_fill_stepsn(
    colors = get_color("OrRd",n=10),
    na.value = NA,
    name="95% CI",
    limits = c(0,0.5),
    breaks = seq(0,0.5,0.1),
    guide = guide_colorbar(direction = "horizontal",
                           nrow=1,
                           title.position="top",
                           barwidth=20
    ))+
  theme_minimal()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 10)
  )
FigS6a

u_sub <- "CUE_Sub_uncertainty.tif"
u <- rast(u_sub)
u <- terra::project(u, crs_robin)

FigS6b <- ggplot()+
  geom_spatraster(data=u, maxcell = 1e7)+
  geom_spatvector(data = global_land_vect, fill=NA, size=0.2)+
  coord_sf(crs = crs_robin)+
  scale_y_continuous(expand = c(0,0), limits = c(-65*10^5, 86*10^5))+
  scale_x_continuous(expand = c(0,0), limits = c(-180*10^5, 180*10^5))+
  labs(title = "Uncertainty of CUE in Subsoil")+
  scale_fill_stepsn(
    colors = get_color("OrRd",n=10),
    na.value = NA,
    name="95% CI",
    limits = c(0,0.5),
    breaks = seq(0,0.5,0.1),
    guide = guide_colorbar(direction = "horizontal",
                           nrow=1,
                           title.position="top",
                           barwidth=20
    ))+
  theme_minimal()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 10))

FigS6b

FigS6 <- ggarrange(FigS6a, FigS6b, common.legend = T, legend = "bottom",
               ncol = 1, 
               labels = c("a)","b)"), 
               font.label = list(size = 12, family = "ARL", face = "plain"))
FigS6
ggsave("figures/figures/FigS6.tiff", width = 6, height = 6, dpi = 300)


















              