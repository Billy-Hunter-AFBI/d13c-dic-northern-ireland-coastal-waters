############################################################
# WFD_D13C_Figures.R
#
# Create publication figures for:
# - Main manuscript
# - Supplementary material
#
# Billy Hunter
############################################################

############################################################
# Libraries
############################################################

library(tidyverse)
library(sf)
library(patchwork)
library(ggrepel)

############################################################
# Directories
############################################################

data_dir    <- "data/raw"
spatial_dir <- "data/spatial/Europe"
fig_dir     <- "figures/main"

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# Import data
############################################################

df <- read_csv(
  file.path(
    data_dir,
    "2020_WFD_del13C_trial.csv"
  ),
  show_col_types = FALSE
)

names(df) <- make.names(names(df))

############################################################
# Import Europe shapefile
############################################################

europe <- st_read(
  file.path(
    spatial_dir,
    "Europe.shp"
  ),
  quiet = TRUE
)

############################################################
# Remove records lacking coordinates
############################################################

df_map <- df %>%
  filter(
    !is.na(Long),
    !is.na(Lat)
  )

############################################################
# Convert stations to sf object
############################################################

stations_sf <- st_as_sf(
  df_map,
  coords = c("Long", "Lat"),
  crs = 4326
)

############################################################
# Output folder
############################################################

dir.create(
  "Figures",
  showWarnings = FALSE
)

############################################################
# Create Helper Fuction to add model stats to plots
############################################################

add_model_stats <- function(response, predictor, data){
  
  mod <- lm(
    as.formula(
      paste(response, "~", predictor)
    ),
    data = data
  )
  
  r2 <- round(
    summary(mod)$r.squared,
    2
  )
  
  pval <- signif(
    coef(summary(mod))[2,4],
    2
  )
  
  paste0(
    "R² = ", r2,
    "\np = ", pval
  )
  
}

############################################################
# FIGURE 1
# Spatial distribution and environmental controls of δ13C-DIC
#
# A = Station map
# B = Temperature
# C = Chlorophyll
# D = DIN
# E = Phosphate
# F = Silicate
# G = Ammonia
# H = DIC
# I = Salinity
############################################################

library(sf)
library(ggplot2)
library(ggrepel)
library(patchwork)

############################################################
# Helper function for plot annotations
############################################################

add_model_stats <- function(response, predictor, data){
  
  mod <- lm(
    as.formula(
      paste(response, "~", predictor)
    ),
    data = data
  )
  
  r2 <- round(summary(mod)$r.squared, 2)
  
  p <- signif(
    coef(summary(mod))[2,4],
    2
  )
  
  paste0(
    "R² = ",
    r2,
    "\nP = ",
    p
  )
  
}

############################################################
# Panel A
# Sampling stations
############################################################

map_panel <- ggplot() +
  
  geom_sf(
    data = europe,
    fill = "grey90",
    colour = "black",
    linewidth = 0.2
  ) +
  
  geom_sf(
    data = stations_sf,
    aes(colour = D13c),
    size = 3
  ) +
  
  geom_text_repel(
    data = df,
    aes(
      Long,
      Lat,
      label = STATION
    ),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  
  coord_sf(
    xlim = c(-7.3,-5.0),
    ylim = c(54.3,55.5),
    expand = FALSE
  ) +
  
  scale_colour_viridis_c(
    name = expression(delta^13*C)
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "right"
  )

############################################################
# Generic plotting function
############################################################

make_panel <- function(xvar,
                       xlabel){
  
  ggplot(
    df,
    aes_string(
      x = xvar,
      y = "D13c"
    )
  ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      hjust = -0.1,
      vjust = 1.1,
      size = 4,
      label = add_model_stats(
        "D13c",
        xvar,
        df
      )
    ) +
    
    theme_bw() +
    
    labs(
      x = xlabel,
      y = expression(delta^13*C)
    )
  
}

############################################################
# Environmental panels
############################################################

p_temp <- make_panel(
  "Temp",
  "Temperature (°C)"
)

p_chl <- make_panel(
  "CHL",
  expression(Chlorophyll~(mg~L^-1))
)

p_din <- make_panel(
  "INORG.N",
  expression(DIN~(mu*mol~L^-1))
)

p_phos <- make_panel(
  "PHOSPHATE",
  expression(Phosphate~(mu*mol~L^-1))
)

p_sil <- make_panel(
  "SILICA",
  expression(Silicate~(mu*mol~L^-1))
)

p_ammonia <- make_panel(
  "AMMONIA",
  expression(Ammonia~(mu*mol~L^-1))
)

p_dic <- make_panel(
  "DIC_molar",
  expression(DIC~(mu*mol~kg^-1))
)

p_sal <- make_panel(
  "Sal",
  "Salinity"
)

############################################################
# Panel tags
############################################################

map_panel  <- map_panel  + labs(tag = "A")

p_temp     <- p_temp     + labs(tag = "B")
p_chl      <- p_chl      + labs(tag = "C")

p_din      <- p_din      + labs(tag = "D")
p_phos     <- p_phos     + labs(tag = "E")
p_sil      <- p_sil      + labs(tag = "F")

p_ammonia  <- p_ammonia  + labs(tag = "G")
p_dic      <- p_dic      + labs(tag = "H")
p_sal      <- p_sal      + labs(tag = "I")

############################################################
# Figure layout
############################################################

figure1 <-
  
  (map_panel  | p_temp    | p_chl)  /
  (p_din      | p_phos    | p_sil)  /
  (p_ammonia  | p_dic     | p_sal)

############################################################
# Global styling
############################################################

figure1 <- figure1 &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

############################################################
# Export
############################################################

ggsave(
  "Figures/Figure1_Main.png",
  figure1,
  width = 15,
  height = 13,
  dpi = 600
)

############################################################
# Supplementary Figure S1
#
# Relationships between dissolved inorganic carbon (DIC)
# concentration and environmental variables across the
# twelve Northern Ireland coastal monitoring stations.
#
# Panels display fitted linear regression models relating
# DIC concentration (µmol kg^-1) to chlorophyll-a,
# dissolved inorganic nitrogen (DIN), phosphate,
# silicate, temperature and salinity. Blue lines indicate
# fitted regressions and shaded areas represent 95%
# confidence intervals. Statistical annotations show the
# coefficient of determination (R²) and associated
# p-values for each model.
#
# Figure exported as:
# Figures/FigureS1_DIC_Relationships.png
############################################################

vars_dic <- c(
  "CHL",
  "INORG.N",
  "PHOSPHATE",
  "SILICA",
  "Temp",
  "Sal"
)

plots_dic <- list()

for(v in vars_dic){
  
  plots_dic[[v]] <-
    
    ggplot(
      df,
      aes_string(
        x = v,
        y = "DIC_molar"
      )
    ) +
    
    geom_point(size = 3) +
    
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.2,
      size = 4,
      label = add_model_stats(
        "DIC_molar",
        v,
        df
      )
    ) +
    
    labs(
      x = case_when(
        v == "CHL" ~ "Chlorophyll (mg L⁻¹)",
        v == "INORG.N" ~ "DIN (µmol L⁻¹)",
        v == "PHOSPHATE" ~ "Phosphate (µmol L⁻¹)",
        v == "SILICA" ~ "Silica (µmol L⁻¹)",
        v == "Temp" ~ "Temperature (°C)",
        v == "Sal" ~ "Salinity"
      ),
      y = expression(DIC~(mu*mol~kg^-1))
    ) +
    
    theme_bw()
  
}

figS2 <-
  
  (plots_dic[["CHL"]] |
     plots_dic[["INORG.N"]] |
     plots_dic[["PHOSPHATE"]]) /
  
  (plots_dic[["SILICA"]] |
     plots_dic[["Temp"]] |
     plots_dic[["Sal"]]) +
  
  plot_annotation(
    tag_levels = "A"
  )

figS2 <- figS2 &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

ggsave(
  "Figures/Supplementary_FigureS2_DIC.png",
  figS2,
  width = 12,
  height = 8,
  dpi = 600
)

############################################################
# Supplementary Figure S2
#
# Sensitivity analysis excluding station WF24.
#
# Relationships between δ13C-DIC and the strongest
# environmental predictors identified in the full dataset
# were re-evaluated following removal of station WF24 to
# assess the influence of this observation on model
# performance and statistical significance.
#
# Panels display fitted linear regression models with
# corresponding 95% confidence intervals. Statistical
# annotations show the coefficient of determination (R²)
# and associated p-values for each model.
#
# Figure exported as:
# Figures/FigureS2_WF24_Sensitivity.png
############################################################

df_noWF24 <- df %>%
  filter(STATION != "WF24")

############################################################
# δ13C vs Temperature
############################################################

p_temp_wf24 <- ggplot(
  df_noWF24,
  aes(Temp, D13c)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.2,
    size = 4,
    label = add_model_stats(
      "D13c",
      "Temp",
      df_noWF24
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Temperature (°C)",
    y = expression(delta^13*C)
  )

############################################################
# δ13C vs DIN
############################################################

p_din_wf24 <- ggplot(
  df_noWF24,
  aes(INORG.N, D13c)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.2,
    size = 4,
    label = add_model_stats(
      "D13c",
      "INORG.N",
      df_noWF24
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = expression(DIN~(mu*mol~L^-1)),
    y = expression(delta^13*C)
  )

############################################################
# δ13C vs DIC
############################################################

p_dic_wf24 <- ggplot(
  df_noWF24,
  aes(DIC_molar, D13c)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.2,
    size = 4,
    label = add_model_stats(
      "D13c",
      "DIC_molar",
      df_noWF24
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = expression(DIC~(mu*mol~kg^-1)),
    y = expression(delta^13*C)
  )

############################################################
# δ13C vs Chlorophyll
############################################################

p_chl_wf24 <- ggplot(
  df_noWF24,
  aes(CHL, D13c)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.2,
    size = 4,
    label = add_model_stats(
      "D13c",
      "CHL",
      df_noWF24
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = expression(Chlorophyll~(mg~L^-1)),
    y = expression(delta^13*C)
  )

############################################################
# Assemble Figure S2
############################################################

figS3 <-
  
  (p_temp_wf24 | p_din_wf24) /
  (p_dic_wf24  | p_chl_wf24) +
  
  plot_annotation(
    tag_levels = "A"
  )

figS3 <- figS3 &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

ggsave(
  filename = "Figures/Supplementary_FigureS3_WF24_Removed.png",
  plot = figS3,
  width = 10,
  height = 8,
  dpi = 600
)
