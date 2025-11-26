library(tidyverse)
library(raster)
library(sf)
library(ggpubr)
library(sp)
library(classInt)
library(cowplot)

colmat <- function(nquantiles = 3, upperleft = "#0096EB", upperright = "#820050", 
                   bottomleft= "#BEBEBE", bottomright = "#FFE60F",
                   xlab = "x label", ylab = "y label", plotLeg = TRUE,
                   saveLeg = TRUE) {
  require(classInt)
  my.data <- seq(0, 1, .01)
  my.class <- classInt::classIntervals(my.data,
                                       n = nquantiles,
                                       style = "quantile"
  )
  my.pal.1 <- findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)
  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i, ] <- findColours(my.class, my.col)
  }
  col.matrix.plot <- col.matrix %>%
    as.data.frame(.) %>%
    mutate("Y" = row_number()) %>%
    mutate_at(.tbl = ., .vars = vars(starts_with("V")), .funs = list(as.character)) %>% 
    gather(data = ., key = X, value = HEXCode, na.rm = FALSE, -Y) %>%
    mutate("X" = as.integer(sub("V", "", .$X))) %>%
    distinct(as.factor(HEXCode), .keep_all = TRUE) %>%
    dplyr::select(-c(4)) %>%
    mutate("X" = rep(seq(from = 1, to = nquantiles, by = 1), each = nquantiles),
           "Y" = rep(seq(from = 1, to = nquantiles, by = 1), times = nquantiles)) %>%
    mutate("UID" = row_number())
  if (plotLeg) {
    p <- ggplot(col.matrix.plot, aes(X, Y, fill = HEXCode)) +
      geom_raster() +
      scale_fill_identity() +
      coord_equal(expand = FALSE) +
      theme_void() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 12, colour = "black",hjust = 0.5, 
                                      vjust = 1),
            axis.title.y = element_text(angle = 90, hjust = 0.5)) +
      xlab(bquote(.(xlab) ~  symbol("\256"))) +
      ylab(bquote(.(ylab) ~  symbol("\256")))
    print(p)
    assign(
      x = "BivLegend",
      value = p,
      pos = .GlobalEnv
    )
  }
  if (saveLeg) {
    ggsave(filename = "bivLegend.pdf", plot = p, device = "pdf",
           path = "./", width = 4, height = 4, units = "in",
           dpi = 300)
  }
  seqs <- seq(0, 100, (100 / nquantiles))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
}

bivariate.map <- function(rasterx, rastery, colormatrix = col.matrix,
                          nquantiles = 3, export.colour.matrix = TRUE,
                          outname = paste0("colMatrix_rasValues", names(rasterx))) {
  quanmean <- getValues(rasterx)
  temp <- data.frame(quanmean, quantile = rep(NA, length(quanmean)))
  brks <- with(temp, quantile(temp,
                              na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))
  ))
  r1 <- within(temp, quantile <- cut(quanmean,
                                     breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE
  ))
  quantr <- data.frame(r1[, 2])
  quanvar <- getValues(rastery)
  temp <- data.frame(quanvar, quantile = rep(NA, length(quanvar)))
  brks <- with(temp, quantile(temp,
                              na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))
  ))
  r2 <- within(temp, quantile <- cut(quanvar,
                                     breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE
  ))
  quantr2 <- data.frame(r2[, 2])
  as.numeric.factor <- function(x) {
    as.numeric(levels(x))[x]
  }
  col.matrix2 <- colormatrix
  cn <- unique(colormatrix)
  for (i in 1:length(col.matrix2)) {
    ifelse(is.na(col.matrix2[i]),
           col.matrix2[i] <- 1, col.matrix2[i] <- which(
             col.matrix2[i] == cn
           )[1]
    )
  }
  if (export.colour.matrix) {
    exportCols <- as.data.frame(cbind(
      as.vector(col.matrix2), as.vector(colormatrix),
      t(col2rgb(as.vector(colormatrix)))
    ))
    colnames(exportCols)[1:2] <- c("rasValue", "HEX")
    assign(
      x = outname,
      value = exportCols,
      pos = .GlobalEnv
    )
  }
  cols <- numeric(length(quantr[, 1]))
  for (i in 1:length(quantr[, 1])) {
    a <- as.numeric.factor(quantr[i, 1])
    b <- as.numeric.factor(quantr2[i, 1])
    cols[i] <- as.numeric(col.matrix2[b, a])
  }
  r <- rasterx
  r[1:length(r)] <- cols
  return(r)
}
nBreaks <- 6

col.matrix <- colmat(nquantiles = nBreaks, xlab = "", ylab = "",
                     upperleft = "#FFD166",    # soft yellow
                     upperright = "#EF476F",   # warm pink/red
                     bottomleft = "#06D6A0",   # mint green
                     bottomright = "#118AB2",  # soft blue
                     # upperleft = "#FFDDD2", upperright = "#FFB4A2",
                     # bottomleft = "#CDEAC0", bottomright = "#A2D2FF",
                     saveLeg = F, plotLeg = T) #using export to export then PPT combine

Top<-raster("CUE_Top_mean_prediction.tif")
Sub<-raster("CUE_Sub_mean_prediction.tif")
rh_resampled <-resample(x=Top,y=Sub)

bivmap <- bivariate.map(rasterx = rh_resampled, Sub,
                        export.colour.matrix = TRUE, outname = "bivMapCols",
                        colormatrix = col.matrix, nquantiles = nBreaks)

bivMapDF <- as.data.frame(bivmap, xy = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::rename("BivValue" = 3) %>%
  gather(key = Variable, value = bivVal, na.rm = FALSE, BivValue)

Fig4a <- ggplot(bivMapDF, aes(x = x, y = y)) +
  geom_raster(aes(fill = bivVal)) +
  scale_fill_gradientn(colours = col.matrix, na.value = "transparent") + 
  theme_bw() +
  coord_quickmap(expand = FALSE) +
  scale_x_continuous(breaks = seq(-120,120, 60), labels = c("120˚W","60˚W","0˚","60˚E","120˚E"))+
  scale_y_continuous(breaks = c(-30,0,30,60), labels = c("30˚S","0˚","30˚N","60˚N"))+
  labs(x=paste('Longitude (\u00b0)'),
       y=paste('Latitude (\u00b0)'))+
  theme(legend.position = "none",
        plot.background = element_blank(),
        # text = element_blank()
  )
Fig4a
##==============================================================================================================
df_top <- read.csv("CUE_Top_mean_prediction.csv")

# load Köppen–Geiger raster (download needed, e.g. from Beck et al. 2018)
kg_raster <- raster("koppen_geiger_0p1.tif")

# turn your data into sf
top_data_sf <- st_as_sf(df_top, coords = c("x","y"), crs = 4326)

# extract climate zone for each point
df_top$climate <- raster::extract(kg_raster, top_data_sf)

legend_raw <- read.delim("legend.txt", header = F, stringsAsFactors = F)

legend_df <- legend_raw %>%
  slice(-1, -2) %>%
  filter(str_detect(V1, "^\\s*\\d+:")) %>%
  mutate(
    code  = as.integer(str_extract(V1, "^\\s*\\d+")),
    class = str_extract(V1, "^\\s*\\d+:\\s*(\\w{1,3})") %>% str_replace("^\\s*\\d+:\\s*", ""),
    rgb   = str_extract(V1, "\\[[0-9 ]+\\]"),
    desc  = V1 %>%
      # remove leading code+class
      str_remove("^\\s*\\d+:\\s*\\w{1,3}\\s*") %>%
      # remove trailing RGB
      str_remove("\\s*\\[[0-9 ]+\\]\\s*$") %>%
      str_trim()
  ) %>%
  dplyr::select(code, class, desc)

legend_df

df_top <- df_top %>%
  left_join(legend_df, by=c("climate"="code"))

top_df <- df_top %>%
  mutate(
    climate_group = str_extract(desc, "^[A-Za-z]+")
  ) %>%
  filter(!is.na(climate_group) & climate_group != "") %>%
  mutate(climate_group = factor(climate_group, levels = c("Tropical","Temperate","Arid","Cold","Polar")))

top_df %>%
  group_by(climate_group) %>%
  summarise(CUE=mean(mean_pred), lower_95 = mean(lower_95), upper_95 = mean(upper_95))

df_sub <- read.csv("CUE_Sub_mean_prediction.csv")

sub_data_sf <- st_as_sf(df_sub, coords = c("x","y"), crs = 4326)

# extract climate zone for each point
df_sub$climate <- raster::extract(kg_raster, sub_data_sf)

df_sub <- df_sub %>%
  left_join(legend_df, by=c("climate"="code"))

sub_df <- df_sub %>%
  mutate(
    climate_group = str_extract(desc, "^[A-Za-z]+")
  ) %>%
  filter(!is.na(climate_group) & climate_group != "") %>%
  mutate(climate_group = factor(climate_group, levels = c("Tropical","Temperate","Arid","Cold","Polar")))

sub_df %>%
  group_by(climate_group) %>%
  summarise(CUE=mean(mean_pred), lower_95 = mean(lower_95), upper_95 = mean(upper_95))

combined_df <- bind_rows(
  top_df %>% mutate(dataset = "Top"),
  sub_df %>% mutate(dataset = "Sub")
)

combined_data <- combined_df %>%
  mutate(Location = case_when(
    y <=30 & y >=-30 ~ "30S-30N",
    y > 30 ~ "> 30N",
    y < -30 ~ "> 30S"
  ))

combined_with_overview <- combined_data %>%
  bind_rows(
    combined_data %>%
      group_by(dataset) %>%
      mutate(Location = "Overview")   # overwrite Location with "Overview"
  ) %>%
  mutate(dataset=factor(dataset, levels = c("Top","Sub"), labels = c("Topsoil","Subsoil")))

combined_with_overview$Location <- factor(
  combined_with_overview$Location,
  levels = c("Overview", "30S-30N", "> 30N", "> 30S")
)

ci_summary <- combined_with_overview %>%
  group_by(Location, dataset) %>%
  summarise(
    mean_pred = mean(mean_pred, na.rm = TRUE),
    lower_95  = mean(lower_95, na.rm = TRUE),
    upper_95  = mean(upper_95, na.rm = TRUE),
    .groups = "drop"
  )

Fig4b <- ggplot(combined_with_overview, aes(
  x = Location,y = mean_pred,fill = dataset)) +
  geom_errorbar(
    data = ci_summary,
    aes(x=Location,
        ymin = lower_95, ymax = upper_95),
    width = 0.2, position = position_dodge(width = 0.6)) +
  # Boxplots
  geom_boxplot(
    width = 0.4, outlier.shape = NA,
    position = position_dodge(width = 0.6),
    coef=0, # whiskers go to min/max but visually ignored
    size=0.5 # thinner outline so CI stands out
  ) +
  geom_vline(xintercept = 1.5, linetype="dashed")+
  theme_bw() +
  labs(x = "", y = "CUE") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.6)) +
  scale_x_discrete(breaks=c("Overview", "30S-30N", "> 30N", "> 30S"),
    labels=c(
    "Overview" = "Global",
    "30S-30N" = "30˚S-30˚N",
    "> 30N" = "> 30˚N",
    "> 30S" = "> 30˚S"  ))+
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  theme(legend.position = c(0.1,0.9),
        legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(0.8,"line"),
        plot.background = element_blank(),
        )
Fig4b
##==============================================================================================================
# load Köppen–Geiger raster (download needed, e.g. from Beck et al. 2018)
bio_raster <- raster("landcover_0.25deg.tif")

# turn your data into sf
top_data_sf <- st_as_sf(df_top, coords = c("x","y"), crs = 4326)

# extract climate zone for each point
df_top$biome <- raster::extract(bio_raster, top_data_sf)
df_top <- df_top %>%
  mutate(biome=round(biome,0)) %>%
  mutate(biome=factor(biome))

legend_raw <- read.delim("LandCoverLabel.txt", header = T, stringsAsFactors = F)

df_top <- df_top %>%
  left_join(legend_raw, by=c("biome"="Value"))

df_top %>%
  group_by(Label) %>%
  summarise(CUE=mean(mean_pred), lower_95 = mean(lower_95), upper_95 = mean(upper_95))

sub_data_sf <- st_as_sf(df_sub, coords = c("x","y"), crs = 4326)

# extract climate zone for each point
df_sub$biome <- raster::extract(bio_raster, sub_data_sf)
df_sub <- df_sub %>%
  mutate(biome=round(biome,0)) %>%
  mutate(biome=factor(biome))

df_sub <- df_sub %>%
  left_join(legend_raw, by=c("biome"="Value"))

combined_df <- bind_rows(
  df_top %>% mutate(dataset = "Topsoil"),
  df_sub %>% mutate(dataset = "Subsoil")
)
combined_df$dataset <- factor(combined_df$dataset, levels = c("Topsoil","Subsoil"))

combined_df <- combined_df %>%
  filter(Label == "Evergreen Needleleaf forest" | Label == "Evergreen Broadleaf forest" |
           Label == "Deciduous Needleleaf forest" | Label == "Deciduous Broadleaf forest"|
           Label == "Mixed forest" | Label == "Woody savannas" | Label == "Savannas"| Label=="Grasslands"|
           Label == "Closed shrublands" | Label == "Open shrublands")

combined_df <- combined_df %>%
  mutate(Type = case_when(
    Label == "Evergreen Needleleaf forest" | Label == "Evergreen Broadleaf forest" |
      Label == "Deciduous Needleleaf forest" | Label == "Deciduous Broadleaf forest"|
      Label == "Mixed forest"  ~ "Forest",
    Label == "Woody savannas" | Label == "Savannas" ~ "Savannas",
    Label == "Closed shrublands" | Label == "Open shrublands" ~ "shrublands",
    Label=="Grasslands" ~ "Grasslands"
  ))

ci_summary <- combined_df %>%
  group_by(Type, dataset) %>%
  summarise(
    mean_pred = mean(mean_pred, na.rm = TRUE),
    lower_95  = mean(lower_95, na.rm = TRUE),
    upper_95  = mean(upper_95, na.rm = TRUE),
    .groups = "drop"
  )

Fig4c <- ggplot(combined_df, aes(
  x = Type,y = mean_pred,fill = dataset)) +
  geom_errorbar(
    data = ci_summary,
    aes(x=Type,
        ymin = lower_95, ymax = upper_95),
    width = 0.2, position = position_dodge(width = 0.6)) +
  # Boxplots
  geom_boxplot(
    width = 0.4, outlier.shape = NA,
    position = position_dodge(width = 0.6),
    coef=0, # whiskers go to min/max but visually ignored
    size=0.5 # thinner outline so CI stands out
  ) +
  theme_bw() +
  labs(x = "", y = "CUE") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.6)) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  theme(legend.position = c(0.1,0.9),
        legend.title = element_blank(), 
        # legend.text = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.8,"line"),
        plot.background = element_blank(),
        # axis.text = element_blank(),
        # axis.title.y = element_blank()
        )
Fig4c

below_row <- ggarrange(Fig4b,Fig4c,
          ncol = 2, 
          labels = c("b)","c)"),
          font.label = list(size = 10, family = "ARL", face = "plain")
          )
below_row

Fig4 <- ggarrange(Fig4a,below_row,
          ncol = 1, 
          labels = c("a)",""),
          font.label = list(size = 10, family = "ARL", face = "plain")
)
Fig4
ggsave("figures/figures/Fig4.tif", height = 8, width = 8, dpi = 300)

combined_df$Label <- factor(combined_df$Label, levels = c(
  "Evergreen Needleleaf forest","Evergreen Broadleaf forest","Deciduous Needleleaf forest","Deciduous Broadleaf forest","Mixed forest",
  "Grasslands","Woody savannas", "Savannas","Closed shrublands","Open shrublands"))

ci_summary2 <- combined_df %>%
  group_by(Label, dataset) %>%
  summarise(
    mean_pred = mean(mean_pred, na.rm = TRUE),
    lower_95  = mean(lower_95, na.rm = TRUE),
    upper_95  = mean(upper_95, na.rm = TRUE),
    .groups = "drop"
  )

FigS3 <- ggplot(combined_df, aes(
  x = Label,y = mean_pred,fill = dataset)) +
  geom_errorbar(
    data = ci_summary2,
    aes(x=Label,
        ymin = lower_95, ymax = upper_95),
    width = 0.2, position = position_dodge(width = 0.6)) +
  # Boxplots
  geom_boxplot(
    width = 0.4, outlier.shape = NA,
    position = position_dodge(width = 0.6),
    coef=0, # whiskers go to min/max but visually ignored
    size=0.5 # thinner outline so CI stands out
  ) +
  # Mean points + CI bars (from summarised data)
  geom_point(
    data = ci_summary2,shape=21,
    aes(x = Label,
        y = mean_pred, fill = dataset),position = position_dodge(width = 0.6),size = 2) +
  theme_bw() +
  labs(x = "", y = "CUE") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.6)) +
  scale_x_discrete(labels=c(
    "Evergreen\nNeedleleaf\nforest",
    "Evergreen\nBroadleaf\nforest",
    "Deciduous\nNeedleleaf\nforest",
    "Deciduous\nBroadleaf\nforest",
    "Mixed\nforest", "Grasslands",
    "Woody\nsavannas", "Savannas",
    "Closed\nshrublands","Open\nshrublands"
  ))+
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  theme(legend.position = c(0.1,0.9),
        legend.title = element_blank(), 
        legend.background = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size=8, family="Arial",color="black"),
        axis.text = element_text(size=8, family="Arial",color="black"),
        axis.title.y = element_text(size=10, family="Arial",color="black"))
FigS3
ggsave("figures/figures/FigS3.tif", height = 5, width = 7, dpi = 300)



