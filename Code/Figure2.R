library(tidyverse)
library(rfPermute)
library(ggplot2)
library(ggpubr)
library(plspm)
library(scales)
library(RColorBrewer)

#import data
data <- read.csv("Env_factors.csv", header = T)
data <- data %>% select(-X)

#check missing values
na_columns <- names(data)[colSums(is.na(data)) > 0]
print(na_columns)

data$depth_class <- as.factor(data$depth_class)

num_cols <- setdiff(names(data),c("depth_class"))
data[num_cols] <- lapply(data[num_cols],as.numeric)

data <- data %>%
  rename(pH = Mean_pH, Sand = Sand.content, Clay = Clay.content, PS = PS_Bio15, TS = TS_Bio4, 
         Silt = Silt.content, BD = Bulk.Density, AGB = AbovegroundBiomass, BGB = BelowgroundBiomass)

##random forest overall soil profile_______________________________________________________________________________
set.seed(629)
overall_data <- data %>%
  mutate(depth_type=case_when(
    depth_class == "0-30" ~ "Topsoil",
    depth_class %in% c("30-60","60-100") ~ "Subsoil"
  )) %>%
  mutate(depth_type=factor(depth_type, levels = c("Topsoil","Subsoil"), labels = c("Topsoil","Subsoil")))

rfP <- rfPermute(CUE ~ depth_type + MAP + MAT + PS + TS + AI + Elevation +
                   Slope + BD + Moisture + pH + CEC + Sand + Clay + Silt + RootDepth +
                   Bedrock + AGB + BGB + LAI + Shannon_EVI + GPP,
                 data=overall_data,ntree=1000, mtry=10, nrep = 9, num.cores = 15,
importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE, rsq=TRUE)

summary(rfP)

mse_overall <- data.frame(
  Trees = 1:rfP$rf$ntree,
  MSE = rfP$rf$mse
)
mse_overall <- mse_overall %>%
  mutate(Layer="Overall")

importance_scale <- data.frame(importance(rfP,scale = T),check.names = F)
importance_scale

importance_scale <- importance_scale[order(importance_scale$`%IncMSE`,decreasing = T),]
importance_scale

importance_scale$name <- rownames(importance_scale)
importance_scale$name <- factor(importance_scale$name, levels = importance_scale$name)

importance_scale$sig_label <- with(importance_scale,
                                       ifelse(`%IncMSE.pval` < 0.001, '***',
                                              ifelse(`%IncMSE.pval` < 0.01, '**',
                                                     ifelse(`%IncMSE.pval` < 0.05, '*', ''))))

num_sig <- sum(importance_scale$sig_label != "")

importance_scale$color <- ifelse(importance_scale$sig_label != "", "#b3cde3", "gray90")

p <- ggplot(importance_scale, aes(x=name, y=`%IncMSE`,fill=color)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_identity() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 2), limits = c(0, 85)) +
  labs(title = "CUE explained by environmental factors overall", y= "Increase in MSE (%)")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust=0, size=12),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none"
  )
print(p)

FigS2 <- p +
  geom_text(data = importance_scale, aes(x = name, y = `%IncMSE`, label = sig_label), nudge_y = 1) +
  annotate('text',x=10, y= 68, label = "Mean of squared residuals: 0.003",color = "black",
           family ="Arial",size=4,fontface="plain",hjust=0)+
  annotate('text',x=10, y=60, label = "Var explained: 87.4%",color = "black",
           family ="Arial",size=4,fontface="plain",hjust=0)
FigS2

ggsave("figures/figures/FigS2.tiff", width = 6.5, height = 5, dpi = 300)
#random forest topsoil____________________________________________________________________________________________________
Top <- data %>%
  filter(depth_class=="0-30")

set.seed(629)
rfP_top <- rfPermute(CUE ~ MAP + MAT + PS + TS + AI + Elevation +
                  Slope + BD + Moisture + pH + CEC + Sand + Clay + Silt + RootDepth +
                  Bedrock + AGB + BGB + LAI + Shannon_EVI + GPP,
                  data=Top,ntree=1000, mtry=10, nrep = 9, num.cores = 15,
                  importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE, rsq=TRUE)

summary(rfP_top)

mse_top <- data.frame(
  Trees = 1:rfP_top$rf$ntree,
  MSE = rfP_top$rf$mse
)
mse_top <- mse_top %>%
  mutate(Layer="Top")

importance_scale_top <- data.frame(importance(rfP_top,scale = T),check.names = F)
importance_scale_top 

importance_scale_top <- importance_scale_top[order(importance_scale_top$`%IncMSE`,decreasing = T),]
importance_scale_top

importance_scale_top$name <- rownames(importance_scale_top)
importance_scale_top$name <- factor(importance_scale_top$name, levels = importance_scale_top$name)
  
importance_scale_top$sig_label <- with(importance_scale_top, 
                                   ifelse(`%IncMSE.pval` < 0.001, '***',
                                   ifelse(`%IncMSE.pval` < 0.01, '**',
                                   ifelse(`%IncMSE.pval` < 0.05, '*', ''))))

num_sig_top <- sum(importance_scale_top$sig_label != "")

importance_scale_top$color <- ifelse(importance_scale_top$sig_label != "", "#b3cde3", "gray90")

p_top <- ggplot(importance_scale_top, aes(x=name, y=`%IncMSE`,fill=color)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_identity() +
    theme_classic() +
    scale_y_continuous(expand = c(0, 2), limits = c(0, 80)) +
    labs(y= "Increase in MSE (%)")+
    theme(panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 9, color = "black"),
          axis.text = element_blank(),
          legend.position = "none"
          )
print(p_top)

Fig2a <- p_top +
  geom_text(data = importance_scale_top, aes(x = name, y = `%IncMSE`, label = sig_label), nudge_y = 1)+
  annotate('text',x=10, y= 63, label = "Mean of squared residuals: 0.003",color = "black",
            family ="Arial",size=3,fontface="plain",hjust=0)+
  annotate('text',x=10, y=55, label = "Var explained: 85.3%",color = "black",
            family ="Arial",size=3,fontface="plain",hjust=0)
Fig2a
#random forest subsoil______________________________________________________________________________________________________________________
Sub <- data %>%
  filter(depth_class!="0-30")

set.seed(629)
rfP_sub <- rfPermute(CUE ~ MAP + MAT + PS + TS + AI + Elevation + 
                       Slope + BD + Moisture + pH + CEC + Sand + Clay + Silt + RootDepth +
                       Bedrock + AGB + BGB + LAI + Shannon_EVI + GPP,
                   data=Sub,ntree=1000, mtry=10, nrep = 9, num.cores = 15,
                   importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE, rsq=TRUE)

summary(rfP_sub)

mse_sub <- data.frame(
  Trees = 1:rfP_sub$rf$ntree,
  MSE = rfP_sub$rf$mse
)
mse_sub <- mse_sub %>%
  mutate(Layer="Sub")

importance_scale_sub <- data.frame(importance(rfP_sub,scale = T),check.names = F)
importance_scale_sub 

importance_scale_sub <- importance_scale_sub[order(importance_scale_sub$`%IncMSE`,decreasing = T),]
importance_scale_sub

importance_scale_sub$name <- rownames(importance_scale_sub)
importance_scale_sub$name <- factor(importance_scale_sub$name, levels = importance_scale_sub$name)

importance_scale_sub$sig_label <- with(importance_scale_sub, 
                                   ifelse(`%IncMSE.pval` < 0.001, '***',
                                          ifelse(`%IncMSE.pval` < 0.01, '**',
                                                 ifelse(`%IncMSE.pval` < 0.05, '*', ''))))

num_sig_sub <- sum(importance_scale_sub$sig_label != "")

importance_scale_sub$color <- ifelse(importance_scale_sub$sig_label != "", "#b3cde3", "gray90")

p_sub <- ggplot(importance_scale_sub, aes(x=name, y=`%IncMSE`,fill=color)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_identity() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 2), limits = c(0, 48)) +
  labs(y= "Increase in MSE (%)")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0,size=12),
        axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9, color = "black"),
        legend.position = "none")
print(p_sub)

Fig2b <- p_sub +
  geom_text(data = importance_scale_sub, aes(x = name, y = `%IncMSE`, label = sig_label), nudge_y = 1)+
  annotate('text',x=10, y= 35, label = "Mean of squared residuals: 0.003",color = "black",
           family ="Arial",size=3,fontface="plain", hjust=0)+
  annotate('text',x=10, y=30, label = "Var explained: 88.6%",color = "black",
           family ="Arial",size=3,fontface="plain",hjust=0)
Fig2b

above_row <- ggarrange(Fig2a, Fig2b,
                       ncol = 2, 
                       labels = c("a) Topsoil", "b) Subsoil"), 
                       font.label = list(size = 14, family = "ARL", face = "plain")
)
above_row

data_rf <- rbind(mse_top,mse_sub)
data_rf <- rbind(data_rf, mse_overall)

data_rf <- data_rf %>%
  mutate(Layer = factor(Layer, levels = c("Overall","Top","Sub"), labels = c("Overall","Topsoil","Subsoil")))

FigS11 <- ggplot(data_rf, aes(x=Trees, y=MSE))+
  facet_wrap(~Layer)+
  geom_line()+
  theme_bw()+
  coord_cartesian(ylim = c(0.0028,0.007), expand = T)+
  labs(x="Number of trees", y="Mean square error")+
  theme(axis.text = element_text(size=8, color="black", family="Arial"))
FigS11

ggsave("figures/figures/FigS11.tiff", width = 7, height = 3.5, dpi = 300)
##path analysis=========================================================================================
LVs <- c("topography","climate","season","texture","vegetation","soil","CUE")

path_matrix <- matrix(0,nrow = length(LVs), ncol = length(LVs), dimnames = list(LVs, LVs))

path_matrix["topography", c("climate","season","texture")] <- 1
path_matrix["climate", c("season","vegetation","CUE")] <- 1
path_matrix["season", c("texture","soil","vegetation","CUE")] <- 1
path_matrix["texture", c("soil","CUE")] <- 1
path_matrix["vegetation", c("soil","CUE")] <- 1
path_matrix["soil",c("CUE")] <- 1

path_matrix <- t(path_matrix)
blocks <- list(
  Topography = c("Elevation", "Slope"),
  climate = c("MAT","MAP"),
  season = c("PS","TS"),
  texture = c("Clay","Sand"),
  Vegetation = c("Shannon_EVI","BGB","GPP"),
  Soil = c("BD","pH","Moisture"),
  CUE = c("CUE")
)
modes <- c("A","A","A","A","B","B","A")
#path analysis topsoil____________________________________________________________________________________________
data_top <- data %>%
  filter(depth_class == "0-30") %>%
  select(-depth_class)

set.seed(629)
pls_top <- plspm(data_top, path_matrix, blocks, modes, scaled = T)
summary(pls_top)

df1 <- pls_top$inner_model$climate 
df2 <- pls_top$inner_model$season
df3 <- pls_top$inner_model$texture
df4 <- pls_top$inner_model$vegetation
df5 <- pls_top$inner_model$soil
df6 <- pls_top$inner_model$CUE
significance_top <- rbind(df1, df2)
significance_top <- rbind(significance_top,df3)
significance_top <- rbind(significance_top,df4)
significance_top <- rbind(significance_top,df5)
significance_top <- rbind(significance_top,df6)
significance_top <- significance_top %>%
  data.frame() %>%
  mutate(P = Pr...t..)

significance_top <- significance_top %>%
  mutate(Significance = case_when(
    P <= 0.001 ~ "***",
    P <= 0.01 ~ "**",
    P <= 0.05 ~ "*",
    TRUE ~ ""
  ))

innerplot(pls_top,txt.col = "black",colpos = "#6890c4BB", colneg = "#f9675dBB")

top_effect <- data.frame(
  Predictor = sub("-> CUE$","",pls_top$effects$relationships[grepl("-> CUE$",pls_top$effects$relationships)]),
  Effect = pls_top$effects$total[grepl("-> CUE$", pls_top$effects$relationships)],
  stringsAsFactors = F
)

top_effect <- top_effect %>%
  mutate(label = sprintf("%.2f",Effect),
         text_y=ifelse(Effect > 0,
                       Effect - 0.02,
                       Effect + 0.02
         ))

top_effect$Predictor <- trimws(top_effect$Predictor)

top_effect$Predictor <- factor(top_effect$Predictor, levels = c("topography","climate","season","texture","vegetation","soil"),
                               labels = c("Topography","Climatic factors","Climatic seasonality","Soil texture","Vegetation","Soil properties"))

Fig2d <- ggplot(top_effect, aes(x=Predictor,y=Effect))+
  geom_hline(yintercept = 0, linewidth=0.5) +
  geom_col(width = 0.5, fill=c("#FF6633","#FFCC66","#660099","#6699CC","#99CC66","#336699"),alpha=.3)+
  theme_classic()+
  labs(x="",y="Standardized total effects")+
  scale_y_continuous(limits = c(-0.5,0.5), expand = c(0,0))+
  scale_x_discrete(labels = label_wrap_gen(width = 5))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust=-2, size=8),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linetype=1, color="black", size=0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none"
  )

Fig2d
#path analysis subsoil_______________________________________________________________________________________________
data_sub <- data %>%
  filter(depth_class != "0-30") %>%
  select(-depth_class)

set.seed(629)
pls_sub <- plspm(data_sub, path_matrix, blocks, modes, scaled = T)
summary(pls_sub)

df1 <- pls_sub$inner_model$climate 
df2 <- pls_sub$inner_model$season
df3 <- pls_sub$inner_model$texture
df4 <- pls_sub$inner_model$vegetation
df5 <- pls_sub$inner_model$soil
df6 <- pls_sub$inner_model$CUE
significance_sub <- rbind(df1, df2)
significance_sub <- rbind(significance_sub,df3)
significance_sub <- rbind(significance_sub,df4)
significance_sub <- rbind(significance_sub,df5)
significance_sub <- rbind(significance_sub,df6)
significance_sub <- significance_sub %>%
  data.frame() %>%
  mutate(P = Pr...t..)

significance_sub <- significance_sub %>%
  mutate(Significance = case_when(
    P <= 0.001 ~ "***",
    P <= 0.01 ~ "**",
    P <= 0.05 ~ "*",
    TRUE ~ ""
  ))

innerplot(pls_sub,txt.col = "black",colpos = "#6890c4BB", colneg = "#f9675dBB")

sub_effect <- data.frame(
  Predictor = sub("-> CUE$","",pls_sub$effects$relationships[grepl("-> CUE$",pls_sub$effects$relationships)]),
  Effect = pls_sub$effects$total[grepl("-> CUE$", pls_sub$effects$relationships)],
  stringsAsFactors = F
)

sub_effect <- sub_effect %>%
  mutate(label = sprintf("%.2f",Effect),
         text_y=ifelse(Effect > 0,
                       Effect - 0.02,
                       Effect + 0.02
         ))

sub_effect$Predictor <- trimws(sub_effect$Predictor)

sub_effect$Predictor <- factor(sub_effect$Predictor, levels = c("topography","climate","season","texture","vegetation","soil"),
                               labels = c("Topography","Climatic factors","Climatic seasonality","Soil texture","Vegetation","Soil properties"))

Fig2f <- ggplot(sub_effect, aes(x=Predictor,y=Effect))+
  geom_hline(yintercept = 0, linewidth=0.5) +
  geom_col(width = 0.5, fill=c("#FF6633","#FFCC66","#660099","#6699CC","#99CC66","#336699"),alpha=.3)+
  theme_classic()+
  labs(x="",y="Standardized total effects")+
  scale_y_continuous(limits = c(-0.7,0.35), expand = c(0,0),breaks = c(-0.7, -0.35, 0, 0.35))+
  scale_x_discrete(labels = label_wrap_gen(width = 5))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust=-2, size=8),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linetype=1, color="black", size=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linetype=1, size=.5, lineend = 1),
        legend.position = "none"
  )

Fig2f

below_row <- ggarrange(Fig2d, Fig2f,
                       ncol = 2, 
                       labels = c("d) Topsoil", "f) Subsoil"), 
                       font.label = list(size = 14, family = "ARL", face = "plain")
)
below_row

Fig2 <- ggarrange(above_row, below_row,
                       ncol = 1, 
                       font.label = list(size = 10, family = "ARL", face = "plain")
)
Fig2

ggsave("figures/figures/Fig2.tiff", width = 10, height = 10, dpi = 300)
##here MAT with CUE=========================================================================================
r1 <- Top %>%
  dplyr::select(depth_class,MAT,CUE)

r2 <- Sub %>%
  dplyr::select(depth_class,MAT,CUE)

data_MAT <- rbind(r1,r2)

data_MAT <- data_MAT %>%
  mutate(Type=case_when(depth_class=="0-30" ~ "Topsoil",
                        TRUE ~"Subsoil")) %>%
  mutate(Type=factor(Type, levels = c("Topsoil","Subsoil")),
         CUE=as.numeric(CUE)) 

ann_text <- data.frame(
  Type = c("Topsoil","Subsoil"),
  MAT=c(-5, -5),
  CUE=c(0.6,0.6),
  label=c("R^2 == 0.01*\"*\"",
          "R^2 == 0.17*\"***\"")
)

ann_text$Type=factor(ann_text$Type, levels = c("Topsoil","Subsoil"))

FigS4 <- ggplot(data=data_MAT,aes(x=MAT,y=CUE))+
  facet_wrap(~Type)+
  theme_bw()+
  geom_point(aes(color=Type,fill=Type),alpha=0.3, size=3,shape=21,stroke=0.3,color="white")+
  stat_smooth(aes(color=Type,fill=Type),method = "lm")+
  scale_fill_manual(values = c('#2A83C7','#E88B33'))+
  scale_color_manual(values = c('#2A83C7','#E88B33'))+
  scale_x_continuous(breaks = c(-5,5,15,25))+
  coord_cartesian(xlim=c(-6,26), expand = T)+
  labs(x="MAT (˚C)")+
  geom_text(
    data = ann_text,
    aes(x=MAT, y=CUE, label = label),
    parse = T,
    hjust = 0,
    size = 4
  )+
  theme(legend.position = "none",
        strip.text = element_text(color="black",size=12,  family="Arial"),
        axis.text=element_text(color="black",size=12,  family="Arial"),   
        axis.title=element_text(color="black",size=12, family="Arial"))
FigS4
ggsave("figures/figures/FigS4.tiff", width = 7.5, height = 4, dpi = 300)


