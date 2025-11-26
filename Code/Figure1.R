library(tidyverse)
library(ggforce)
library(rstatix)
library(ggpubr)
library(ggpmisc)
library(plotbiomes)
library(extrafont)

#import data
data <- read.csv("DATA.csv")

#statistical to test CUE 
model <- lmerTest::lmer(CUE ~ depth_class + (1|Paper.ID), data=data)
summary(model)
anova(model)
em <- emmeans::emmeans(model,  ~ depth_class)
data_letter <- multcomp::cld(em, Letters= letters, adjudt="tukey")

sample_counts <- data %>%
  group_by(depth_class) %>%
  summarise(n=n()) %>%
  mutate(x=c(1,2,3))

FigS1a <- data %>%
  ggplot(aes(x=depth_class, y=CUE)) +
  geom_sina(aes(color=depth_class),alpha=0.3)+
  geom_boxplot(width=0.2,outlier.shape = NA,aes(color=depth_class))+
  geom_text(data = sample_counts, 
            aes(x=x+0.25, y=0.05,label = paste0("(", n,")"),color=depth_class),
            inherit.aes = FALSE, size = 3) +
  geom_text(data = data_letter, aes(x=depth_class, label=.group, y=0.61))+
  scale_fill_manual(values = c('#2A83C7','#E88B33',"#8da0cb"))+
  scale_color_manual(values = c('#2A83C7','#E88B33',"#8da0cb"))+
  coord_cartesian(ylim = c(0,0.7))+
  theme_bw()+
  annotate("text",
           x = 2,
           y = 0.7,
          label = "italic(F) == 11.67 * \"; \"* italic(P) * \" < \" * 0.001",
          size = 4,
           parse=T) +
  labs(x="Soil Depth (cm)", y="CUE")+
  theme(legend.position = "none",
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title = element_text(color="black",size=12, family="Arial")
  )

FigS1a

data_TD <- data %>%
  mutate(depth_type=case_when(
    depth_class == "0-30" ~ "Topsoil",
    depth_class %in% c("30-60","60-100") ~ "Subsoil"
  )) %>%
  mutate(depth_type=factor(depth_type, levels = c("Topsoil","Subsoil"), labels = c("Topsoil","Subsoil")))

sample_counts_TD <- data_TD %>%
  group_by(depth_type) %>%
  summarise(n=n()) %>%
  mutate(x=c(1,2))

#statistical to test CUE 
model <- lmerTest::lmer(CUE ~ depth_type + (1|Paper.ID), data=data_TD)
summary(model)
anova(model)
em <- emmeans::emmeans(model,  ~ depth_type)
data_letter <- multcomp::cld(em, Letters= letters, adjudt="tukey")

Fig1b <- data_TD %>%
  ggplot(aes(x=depth_type, y=CUE)) +
  geom_sina(aes(color=depth_type),alpha=0.3)+
  geom_boxplot(width=0.2,outlier.shape = NA,aes(color=depth_type))+
  geom_text(data = sample_counts_TD, 
            aes(x=x+0.25, y=0.05,label = paste0("(", n,")"),color=depth_type),
            inherit.aes = FALSE, size = 3) +
  geom_text(data = data_letter, aes(x=depth_type, label=.group, y=0.61))+
  scale_fill_manual(values = c('#2A83C7','#E88B33'))+
  scale_color_manual(values = c('#2A83C7','#E88B33'))+
  coord_cartesian(ylim = c(0,0.7))+
  theme_bw()+
  annotate("text",
           x = 1.5,
           y = 0.7,
           label = "italic(F) == 22.21 * \"; \"* italic(P) * \" < \" * 0.001",
           size = 4,
           parse=T) +
  labs(x="Soil layer", y="CUE")+
  theme(legend.position = "none",
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title = element_text(color="black",size=12, family="Arial")
  )

Fig1b

data_TD_climate <- data_TD %>%
  mutate(Climatezone = case_when(
           Climatezone %in% c("1","2") ~ "Tropical",
           Climatezone %in% c("4","5","6","7") ~ "Arid",
           Climatezone %in% c("11","12","14","15") ~ "Temperate",
           Climatezone %in% c("21","22","23") ~ "Cold",
           Climatezone == "29" ~ "Cold",
           TRUE ~ NA_character_
           ))

data_TD_climate <- data_TD_climate %>%
  mutate(Climatezone=factor(Climatezone, levels = c("Tropical","Temperate","Arid","Cold")))

data_TD_climate %>%
  group_by(Climatezone,depth_type) %>%
  summarise(n=n())

stat.test <- data_TD_climate %>%
  group_by(Climatezone) %>%
  wilcox_test(CUE ~ depth_type) %>%
  adjust_pvalue(method="bonferroni") %>%
  add_significance("p.adj")

Fig1d <- ggplot(data_TD_climate, aes(x=Climatezone, y=CUE)) +
  theme_bw()+
  geom_boxplot(aes(color=depth_type), position = position_dodge(0.9), width=0.3)+
  stat_compare_means(aes(group = depth_type), label = "p.signif")+
  geom_sina(aes(color=depth_type),alpha=0.3)+
  scale_fill_manual(values = c('#2A83C7','#E88B33'))+
  scale_color_manual(values = c('#2A83C7','#E88B33'))+
  coord_cartesian(ylim = c(0,0.7))+
  labs(x="Climate zone", y="CUE")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title = element_text(color="black",size=12, family="Arial")
  )
Fig1d

data$Latitude<-as.numeric(data$Latitude)

data_LAT <- data %>%
  mutate(Latitude = abs(Latitude)) %>%
  mutate(depth_class=factor(depth_class, levels = c("0-30", "30-60","60-100")))

FigS1b <- data_LAT %>%
  ggplot(aes(x=Latitude, y=CUE, fill = depth_class,color = depth_class))+
  theme_bw()+
  geom_point(alpha=0.3, size=3,shape=21,stroke=0.3,color="white") +
  stat_poly_line(formula = y ~ poly(x,2), se=T)+
  stat_poly_eq(formula = y ~ poly(x,2),
               aes(label=paste(..eq.label..,
                               ..adj.rr.label..,
                               ..p.value.label..,
                               sep = "~~~")), parse = T, vstep = 0.065,
               label.x.npc = 'left', label.y.npc = 'top', size=3) +
  scale_fill_manual(values = c('#2A83C7','#E88B33',"#8da0cb"))+
  scale_color_manual(values = c('#2A83C7','#E88B33',"#8da0cb"))+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by=20))+
  scale_y_continuous(limits = c(0,0.7), breaks = seq(0,0.7, by=0.2))+
  labs(x=paste('Absolute latitude (\u00b0)'), y="CUE")+
  theme(legend.position = "none",
        strip.text = element_text(color="black",size=12,  family="Arial"),
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title=element_text(color="black",size=12, family="Arial"))
FigS1b   

data_TD_2 <- data_LAT %>%
  mutate(depth_type=case_when(
    depth_class == "0-30" ~ "Topsoil",
    depth_class %in% c("30-60","60-100") ~ "Subsoil"
  )) %>%
  mutate(depth_type=factor(depth_type, levels = c("Topsoil","Subsoil")))

Fig1c <- data_TD_2 %>%
  ggplot(aes(x=Latitude, y=CUE, fill = depth_type,color = depth_type))+
  theme_bw()+
  geom_point(alpha=0.3, size=3,shape=21,stroke=0.3,color="white") +
  stat_poly_line(formula = y ~ poly(x,2), se=T)+
  stat_poly_eq(formula = y ~ poly(x,2),
               aes(label=paste(..eq.label..,
                               ..adj.rr.label..,
                               ..p.value.label..,
                               sep = "~~~")), parse = T, vstep = 0.065,
               label.x.npc = 'left', label.y.npc = 'top', size=3) +
  scale_fill_manual(values = c('#2A83C7','#E88B33'))+
  scale_color_manual(values = c('#2A83C7','#E88B33'))+
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60, by=20))+
  scale_y_continuous(limits = c(0,0.7), breaks = seq(0,0.7, by=0.2))+
  labs(x=paste('Absolute latitude (\u00b0)'), y="CUE")+
  theme(legend.position = "none",
        strip.text = element_text(color="black",size=12,  family="Arial"),
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title=element_text(color="black",size=12, family="Arial"))
Fig1c  

data <- data %>%
  mutate(
  Latitude = as.numeric(Latitude),
  Longitude = as.numeric(Longitude),
  CUE = as.numeric(CUE),
)

data <- data %>%
  mutate(depth_type=case_when(
    depth_class == "0-30" ~ "Topsoil",
    depth_class %in% c("30-60","60-100") ~ "Subsoil"
  )) %>%
  mutate(depth_type=factor(depth_type, levels = c("Topsoil","Subsoil")))

col <- c("#2A83C7","#E88B33")

Fig1a <- whittaker_base_plot()+
  theme_bw()+
  geom_point(data = data, aes(x=MAT, y=MAP/10, color = depth_type), size=3, stroke=0.5, shape=19)+
  labs(y="Precipitation (cm)", color="Soil layer")+
  scale_y_continuous(limits = c(0,460))+
  scale_color_manual(values = col)+
  theme(panel.background = element_blank(),
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.title = element_text(family = "Arial", size=10,color="black"),
        legend.text = element_text(family = "Arial", size=8,color="black"),
        axis.text = element_text(family = "Arial", size=10,color="black"))

Fig1a

below_row <- ggarrange(Fig1b, Fig1c,
                       ncol = 2, 
                       labels = c("b)", "c)"), 
                       font.label = list(size = 14, family = "ARL", face = "plain")
                       )

Fig.1.1 <- ggarrange(Fig1a, below_row, 
                  ncol = 1, 
                  labels = "a)", 
                  font.label = list(size = 14, family = "ARL", face = "plain")
                  )

Fig.1 <- ggarrange(Fig.1.1, Fig1d, 
                   ncol = 1, 
                   labels = c("","d)"), 
                   heights = c(1, 0.5),  # You can adjust relative heights
                   font.label = list(size = 14, family = "ARL", face = "plain")
)

print(Fig.1)

ggsave("figures/figures/Fig1.tiff", width = 8, height = 10, dpi = 300)

FigS1<- ggarrange(FigS1a, FigS1b,
                       ncol = 2, 
                       labels = c("a)", "b)"), 
                       font.label = list(size = 12, family = "ARL", face = "plain"))

print(FigS1)

ggsave("figures/figures/FigS1.tiff", width = 7.5, height = 4, dpi = 300)

##here plot world map for supplementary=========================================================================================

world_map <- map_data('world')

map1<-ggplot()+
  geom_polygon(data = world_map, aes(x=long, y=lat, group = group), color=NA,
               fill='#d9d9d9')+
  theme_bw()+
  geom_point(data = data,size=3,
             aes(y = Latitude, x = Longitude, color= CUE),stroke=0.1,alpha=0.3)+
  scale_color_gradient(low="#fef0d9", high = "#b30000")+
  scale_x_continuous(breaks = seq(-180,180,30))+
  scale_y_continuous(breaks = seq(-90,90,30))+
  coord_cartesian(xlim = c(-155,175), ylim = c(-55,80))+
  labs(color="CUE",
       x=paste('Longitude (\u00b0)'),
       y=paste('Latitude (\u00b0)'))+
  theme(legend.position = c(0.1,0.3),
        legend.background = element_blank(),
        legend.title = element_text(margin = margin(b=1), color="black",size=10, family="Arial"),
        legend.text = element_text(color="black",size=10, family="Arial"),
        panel.grid = element_blank(),
        axis.text = element_text(color="black",size=12, family="Arial"),
        axis.title = element_text(color="black",size=12, family="Arial")
  )
map1
ggsave("figures/figures/FigS5.tiff", width = 7.5, height = 4, dpi = 300)









