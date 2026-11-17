############################################################
# DIC Stable Isotope Trial
# Billy Hunter
#
# Purpose:
# Explore relationships between DIC, δ13C,
# chlorophyll and nutrients using simple linear models.
#
# Sample size is small (~11-12 stations), therefore:
# - Simple linear regression only
# - No multiple regression
# - No interaction terms
############################################################
library(tidyverse)
library(broom)

#------------------------------------------------------------
# Project directories
#------------------------------------------------------------

data_dir      <- "data/raw"
fig_dir       <- "figures/main"
diag_dir      <- "figures/diagnostics"
results_dir   <- "results"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------
# Import data
#------------------------------------------------------------

df <- read_csv(
  file.path(data_dir, "2020_WFD_del13C_trial.csv"),
  show_col_types = FALSE
)

#------------------------------------------------------------
# Clean names
#------------------------------------------------------------

names(df) <- make.names(names(df))

#------------------------------------------------------------
# Convert variables to numeric
#------------------------------------------------------------

vars <- c(
  "D13c",
  "DIC_molar",
  "CHL",
  "INORG.N",
  "PHOSPHATE",
  "SILICA",
  "AMMONIA",
  "Sal",
  "Temp"
)

df[vars] <- lapply(df[vars], as.numeric)

#------------------------------------------------------------
# Full dataset
#------------------------------------------------------------

df_all <- df

#------------------------------------------------------------
# Sensitivity dataset (WF24 removed)
#------------------------------------------------------------

df_noWF24 <- df %>%
  filter(STATION != "WF24")

#------------------------------------------------------------
# Create folder for figures
#------------------------------------------------------------

dir.create(
  "Regression_Plots",
  showWarnings = FALSE
)

#------------------------------------------------------------
# Create Folder for Diagnostic Plots
#------------------------------------------------------------

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


#------------------------------------------------------------
# Correlation matrix
#------------------------------------------------------------
cor_matrix <- cor(
  df[, vars],
  use = "pairwise.complete.obs"
)

round(cor_matrix, 2)

write.csv(
  cor_matrix,
  "Correlation_Matrix.csv"
)

#------------------------------------------------------------
# Define Hypotheses
#------------------------------------------------------------

models <- tribble(
  ~Response, ~Predictor,
  
  "D13c",      "DIC_molar",
  "D13c",      "CHL",
  "D13c",      "INORG.N",
  "D13c",      "SILICA",
  "D13c",      "Temp",
  "D13c",      "Sal",
  
  "DIC_molar", "CHL",
  "DIC_molar", "INORG.N",
  "DIC_molar", "PHOSPHATE",
  "DIC_molar", "SILICA",
  "DIC_molar", "Temp",
  "DIC_molar", "Sal",
  
  "CHL",       "INORG.N",
  "CHL",       "PHOSPHATE",
  "CHL",       "SILICA",
  "CHL",       "Temp",
  "CHL",       "Sal"
)

############################################################
# Run all simple linear regressions
############################################################

run_lm_analysis <- function(data, output_dir){
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(output_dir, "Regression_Plots"),
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(output_dir, "Diagnostics"),
    showWarnings = FALSE
  )
  
  ##########################################################
  # Model summaries
  ##########################################################
  
  results <- list()
  
  for(i in 1:nrow(models)){
    
    response  <- models$Response[i]
    predictor <- models$Predictor[i]
    
    mod <- lm(
      as.formula(
        paste(response, "~", predictor)
      ),
      data = data
    )
    
    stats <- glance(mod)
    coefs <- tidy(mod)
    
    results[[i]] <- data.frame(
      Response  = response,
      Predictor = predictor,
      N         = nobs(mod),
      Slope     = coefs$estimate[2],
      R2        = stats$r.squared,
      Adj_R2    = stats$adj.r.squared,
      P_value   = coefs$p.value[2]
    )
  }
  
  results_df <- bind_rows(results)
  
  write.csv(
    results_df,
    file.path(
      output_dir,
      "Linear_Model_Summary.csv"
    ),
    row.names = FALSE
  )
  
  ##########################################################
  # Regression plots
  ##########################################################
  
  for(i in 1:nrow(models)){
    
    response  <- models$Response[i]
    predictor <- models$Predictor[i]
    
    mod <- lm(
      as.formula(
        paste(response, "~", predictor)
      ),
      data = data
    )
    
    r2 <- round(
      summary(mod)$r.squared,
      3
    )
    
    p <- signif(
      coef(summary(mod))[2,4],
      3
    )
    
    g <- ggplot(
      data,
      aes_string(
        x = predictor,
        y = response
      )
    ) +
      geom_point(size = 3) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        colour = "blue"
      ) +
      theme_bw() +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        hjust = 1.1,
        vjust = 1.2,
        label = paste0(
          "R² = ", r2,
          "\np = ", p
        )
      ) +
      labs(
        x = predictor,
        y = response
      )
    
    ggsave(
      filename = file.path(
        output_dir,
        "Regression_Plots",
        paste0(
          response,
          "_vs_",
          predictor,
          ".png"
        )
      ),
      plot = g,
      width = 5,
      height = 4,
      dpi = 300
    )
    
    ########################################################
    # Diagnostics
    ########################################################
    
    png(
      filename = file.path(
        output_dir,
        "Diagnostics",
        paste0(
          response,
          "_vs_",
          predictor,
          "_diagnostics.png"
        )
      ),
      width = 1800,
      height = 1600,
      res = 250
    )
    
    par(
      mfrow = c(2,2),
      mar = c(5,5,3,2)
    )
    
    plot(mod)
    
    mtext(
      paste0(
        response,
        " ~ ",
        predictor,
        "   R² = ",
        r2,
        "   p = ",
        p
      ),
      outer = TRUE,
      line = -2,
      cex = 1.2
    )
    
    dev.off()
    
  }
  
}



############################################################
# Plot all models
############################################################

datasets <- list(
  All_Stations = df_all,
  No_WF24 = df_noWF24
)

for(dataset_name in names(datasets)) {
  
  current_df <- datasets[[dataset_name]]
  
  output_folder <- file.path(
    "Regression_Plots",
    dataset_name
  )
  
  dir.create(
    output_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for(i in 1:nrow(models)) {
    
    response <- models$Response[i]
    predictor <- models$Predictor[i]
    
    formula_text <- paste(response, "~", predictor)
    
    mod <- lm(
      as.formula(formula_text),
      data = current_df
    )
    
    r2 <- round(summary(mod)$r.squared, 2)
    
    p <- signif(
      summary(mod)$coefficients[2,4],
      3
    )
    
    g <- ggplot(
      current_df,
      aes_string(
        x = predictor,
        y = response
      )
    ) +
      geom_point(size = 3) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        colour = "blue"
      ) +
      theme_bw(base_size = 12) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "n = ", nrow(current_df),
          "\nR² = ", r2,
          "\np = ", p
        ),
        hjust = 1.1,
        vjust = 1.5
      ) +
      labs(
        x = predictor,
        y = response
      )
    
    ggsave(
      filename = file.path(
        output_folder,
        paste0(
          response,
          "_vs_",
          predictor,
          ".png"
        )
      ),
      plot = g,
      width = 5,
      height = 4,
      dpi = 300
    )
    
  }
  
}


############################################################
# Generate diagnostic plots for all models
############################################################

datasets <- list(
  All_Stations = df_all,
  No_WF24 = df_noWF24
)

for(dataset_name in names(datasets)) {
  
  current_df <- datasets[[dataset_name]]
  
  output_dir <- file.path(
    "Diagnostics",
    dataset_name
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for(i in 1:nrow(models)) {
    
    response  <- models$Response[i]
    predictor <- models$Predictor[i]
    
    mod <- lm(
      as.formula(
        paste(response, "~", predictor)
      ),
      data = current_df
    )
    
    r2 <- round(
      summary(mod)$r.squared,
      3
    )
    
    pval <- signif(
      coef(summary(mod))[2,4],
      3
    )
    
    png(
      filename = file.path(
        output_dir,
        paste0(
          response,
          "_vs_",
          predictor,
          "_diagnostics.png"
        )
      ),
      width = 1800,
      height = 1600,
      res = 250
    )
    
    par(
      mfrow = c(2,2),
      mar = c(5,5,3,2)
    )
    
    plot(mod)
    
    mtext(
      paste0(
        dataset_name,
        "\n",
        response,
        " ~ ",
        predictor,
        "   R² = ",
        r2,
        "   p = ",
        pval
      ),
      outer = TRUE,
      line = -2,
      cex = 1.3
    )
    
    dev.off()
    
  }
  
}

