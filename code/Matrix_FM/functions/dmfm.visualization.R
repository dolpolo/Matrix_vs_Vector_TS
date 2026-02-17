################################################################################
########################## DESCRIPTIVE EA ANALYSIS #############################
################################################################################

# ==============================================================================
# GDP COMMUNALITIES IN EA 
# ==============================================================================
GDP_communality <- function(country_codes) {
  
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(zoo)
  
  # Create folder for saving graphs
  dir.create("Figures", showWarnings = FALSE)
  
  # Define recession periods
  recession_df <- data.frame(
    start = as.Date(c("2008-01-01", "2011-07-01", "2020-01-01")),
    end   = as.Date(c("2009-06-30", "2013-03-31", "2020-06-30"))
  )
  
  all_gdp <- list()
  
  for (cc in country_codes) {
    data_path <- paste0("Data/", cc, "/Original/", cc, "DataQ_LT.xlsx")
    gdp_label <- paste0("GDP_", cc)
    
    df <- read_excel(data_path)
    
    if (cc == "EA") {
      EA_corr_data <- dplyr::select(df, -Time)
      correlation_matrix <- cor(EA_corr_data, use = "pairwise.complete.obs")
      melted_corr <- melt(correlation_matrix)
      
      # Crop long labels for better readability
      label_crop <- function(x) {
        sapply(x, function(v) substr(v, 1, 15))
      }
      
      p_heat <- ggplot(melted_corr, aes(x = label_crop(Var1), y = label_crop(Var2), fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
        theme_bw(base_size = 12) +
        labs(title = "Correlation Heatmap of EA Variables", x = NULL, y = NULL) +
        theme(
          axis.text.x = element_text(angle = 90, size = 6, hjust = 1),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(hjust = 0.5)
        )
      
      ggsave("Figures/EA_Correlation_Heatmap.png", p_heat, width = 10, height = 8, bg = "white")
    }
    
    gdp_idx <- which(tolower(colnames(df)) == tolower(gdp_label))
    gdp_series <- df[[gdp_idx]] * 100
    date_series <- df[[1]]
    
    gdp_df <- data.frame(
      date = as.Date(date_series),
      gdp = gdp_series,
      country = cc
    )
    
    all_gdp[[cc]] <- gdp_df
  }
  
  full_gdp <- bind_rows(all_gdp)
  
  # Filter out EA for plots and set country order
  gdp_filtered <- full_gdp %>% filter(country != "EA")
  gdp_filtered$country <- factor(gdp_filtered$country, levels = c("DE", "FR", "IT", "ES"))
  
  # Line plot with recession shading
  p1 <- ggplot(gdp_filtered, aes(x = date, y = gdp, color = country)) +
    geom_rect(data = recession_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey", alpha = 0.2) +
    geom_line(linewidth = 1) +
    labs(title = "GDP Time Series By Country", x = "Date", y = "GDP") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  ggsave("Figures/GDP_TimeSeries_SelectedCountries.png", p1, width = 10, height = 6, bg = "white")
  
  # Faceted plot for DE, FR, IT, ES with fixed layout
  target_countries <- c("DE", "FR", "IT", "ES")
  gdp_faceted <- gdp_filtered %>% filter(country %in% target_countries)
  
  p2 <- ggplot(gdp_faceted, aes(x = date, y = gdp)) +
    geom_rect(data = recession_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey", alpha = 0.2) +
    geom_line(color = "steelblue", linewidth = 1) +
    facet_wrap(~ country, ncol = 2, scales = "free_y") +
    labs(title = "GDP Time Series By Country", x = "Date", y = "GDP (Index)") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("Figures/GDP_TimeSeries_Facet_Selected.png", p2, width = 10, height = 7, bg = "white")
  
  # Comovement heatmap
  wide_gdp <- gdp_filtered %>%
    filter(!is.na(gdp)) %>%
    pivot_wider(names_from = country, values_from = gdp)
  
  gdp_matrix <- as.matrix(wide_gdp[,-1])
  rownames(gdp_matrix) <- wide_gdp$date
  
  comov_matrix <- cor(gdp_matrix, use = "pairwise.complete.obs")
  comov_melt <- melt(comov_matrix)
  colnames(comov_melt) <- c("Country1", "Country2", "Correlation")
  
  # Set factor levels for ordering in heatmap
  comov_melt$Country1 <- factor(comov_melt$Country1, levels = c("DE", "FR", "IT", "ES"))
  comov_melt$Country2 <- factor(comov_melt$Country2, levels = c("DE", "FR", "IT", "ES"))
  
  p3 <- ggplot(comov_melt, aes(x = Country1, y = Country2, fill = Correlation)) +
    geom_tile(color = "black", linewidth = 0.3) +
    geom_text(aes(label = round(Correlation, 2)), size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1), name = "Correlation") +
    labs(title = "GDP Correlation Across Countries", x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5))
  
  ggsave("Figures/GDP_Correlation_Heatmap.png", p3, width = 8, height = 6, bg = "white")
  
  print(p_heat)
  print(p1)
  print(p2)
  print(p3)
  
  
  # Export for inspection
  assign("gdp_results", list(
    data = all_gdp,
    combined = full_gdp,
    comovement = comov_matrix
  ), envir = .GlobalEnv)
}


# ==============================================================================
# GDP COMMUNALITIES IN EA 
# ==============================================================================
plot_lagged_heatmap_grid <- function(lagged_corrs, focus_countries = c("DE", "FR", "IT", "ES")) {
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  
  # Estrae solo i lag da 1 a 12
  all_corrs <- bind_rows(lapply(names(lagged_corrs), function(lag_label) {
    lag_num <- as.integer(gsub("Lag_", "", lag_label))
    if (lag_num >= 1 && lag_num <= 12) {
      mat <- lagged_corrs[[lag_label]]
      df <- melt(mat)
      df$Lag <- lag_num
      return(df)
    } else {
      return(NULL)
    }
  }))
  
  colnames(all_corrs) <- c("Country1", "Country2", "Correlation", "Lag")
  all_corrs$Country1 <- as.character(all_corrs$Country1)
  all_corrs$Country2 <- factor(as.character(all_corrs$Country2), levels = c("DE", "FR", "IT", "ES"))
  
  # Heatmap per ogni paese
  make_plot <- function(ref_country) {
    df_plot <- all_corrs %>% filter(Country1 == ref_country)
    
    ggplot(df_plot, aes(x = Lag, y = Country2, fill = Correlation)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-1, 1), name = "Corr") +
      scale_x_continuous(breaks = c(3, 6, 9, 12)) +
      labs(title = ref_country, x = "Lag", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
  }
  
  plots <- lapply(focus_countries, make_plot)
  names(plots) <- focus_countries
  
  # Estrai legenda da un plot
  get_legend_plot <- function() {
    df_plot <- all_corrs %>% filter(Country1 == "DE")
    p <- ggplot(df_plot, aes(x = Lag, y = Country2, fill = Correlation)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-1, 1), name = "Correlation") +
      scale_x_continuous(breaks = c(3, 6, 9, 12)) +
      labs(x = "Lag", y = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
    cowplot::get_legend(p)
  }
  
  shared_legend <- get_legend_plot()
  
  # Griglia 2x2
  plot_grid_top <- cowplot::plot_grid(
    plots$DE, plots$FR, plots$IT, plots$ES,
    ncol = 2, align = "hv"
  )
  
  final_plot <- cowplot::plot_grid(
    plot_grid_top,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.08)
  )
  
  dir.create("Figures/LaggedByCountry", showWarnings = FALSE)
  ggsave("Figures/LaggedByCountry/Lagged_Correlation_Grid_SharedLegend.png",
         final_plot, width = 14, height = 10, bg = "white")
  
  print(final_plot)
}


GDP_lagged_correlation <- function(wide_gdp, max_lag = 12) {
  library(reshape2)
  library(ggplot2)
  library(gridExtra)
  
  # Controllo colonne numeriche e rimuove righe NA
  gdp_matrix <- as.matrix(wide_gdp[ , -1])
  rownames(gdp_matrix) <- as.character(wide_gdp$date)
  
  lagged_corr_list <- list()
  
  for (lag in 0:max_lag) {
    lagged_data <- gdp_matrix
    for (j in 1:ncol(gdp_matrix)) {
      lagged_data[, j] <- dplyr::lag(gdp_matrix[, j], lag)
    }
    
    # Solo righe complete
    valid_idx <- complete.cases(gdp_matrix, lagged_data)
    corr_mat <- cor(gdp_matrix[valid_idx, , drop = FALSE],
                    lagged_data[valid_idx, , drop = FALSE],
                    use = "pairwise.complete.obs")
    
    lagged_corr_list[[paste0("Lag_", lag)]] <- corr_mat
  }
  
  # Facoltativo: plot heatmap per ogni lag
  plots <- lapply(names(lagged_corr_list), function(name) {
    melted <- melt(lagged_corr_list[[name]])
    ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           limits = c(-1, 1), name = "Corr") +
      labs(title = name, x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Salva figura combinata (opzionale)
  ggsave("Figures/GDP_Lagged_Correlation_Grid.png",
         marrangeGrob(plots, nrow = 3, ncol = 3), width = 15, height = 10)
  
  return(lagged_corr_list)
}


# ==============================================================================
# GDP BY COUNTRY CONTRIBUTION
# ==============================================================================
plot_euro_area_gdp <- function(file_path = "Data/GDP_millionsEA.xlsx",
                               sheet_name = "GDP",
                               quarter = "2024-Q4",
                               special_countries = c("Italy", "Germany", "France", "Spain", "Ireland",
                                                     "Greece", "Portugal","Netherlands", "Belgium", "Austria"),
                               show_table = TRUE) {
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(knitr)
  library(kableExtra)
  
  # 1. Importa e pulisci dati
  gdp_data_raw <- read_excel(file_path, sheet = sheet_name)
  
  gdp_data <- gdp_data_raw %>%
    dplyr::select(Country = 1, value = !!sym(quarter)) %>%
    mutate(
      Country = case_when(
        Country == "Germany (until 1990 former territory of the FRG)" ~ "Germany",
        Country == "Euro area" ~ NA_character_,
        TRUE ~ Country
      ),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(!is.na(Country), !Country %in% c(
      "European Union - 27 countries (from 2020)",
      "Euro area – 20 countries (from 2023)")
    )
  
  # 2. Mappa base
  world <- ne_countries(scale = "large", returnclass = "sf")
  europe <- world %>% filter(region_un == "Europe")
  
  euro_area_countries <- c(
    "Austria", "Belgium", "Croatia", "Cyprus", "Estonia", "Finland",
    "France", "Germany", "Greece", "Ireland", "Italy", "Latvia",
    "Lithuania", "Luxembourg", "Malta", "Netherlands", "Portugal",
    "Slovakia", "Slovenia", "Spain"
  )
  
  # 3. Calcolo share
  gdp_total_ea <- gdp_data %>%
    filter(Country %in% euro_area_countries) %>%
    summarise(total = sum(value, na.rm = TRUE)) %>%
    pull(total)
  
  gdp_data <- gdp_data %>%
    mutate(GDP_Share_EA = (value / gdp_total_ea) * 100)
  
  # 4. Join mappa + dati
  europe_gdp <- europe %>%
    left_join(gdp_data, by = c("name" = "Country")) %>%
    mutate(
      EA_Member = name %in% euro_area_countries,
      SpecialCountry = name %in% special_countries,
      FillColor = case_when(
        EA_Member ~ GDP_Share_EA,
        TRUE ~ NA_real_
      )
    )
  
  # 5. Plot
  gdp_plot <- ggplot() +
    theme_void() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5,
                                family = "serif", margin = margin(b = 10)),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "aliceblue", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      legend.position = "none",
      plot.caption = element_text(size = 9, hjust = 1, family = "serif"),
      plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
    ) +
    
    # Mappa
    geom_sf(data = europe, fill = "grey95", color = "black", size = 0.1) +
    geom_sf(data = europe_gdp %>% filter(!EA_Member), fill = "grey85", color = "black", size = 0.1) +
    geom_sf(data = europe_gdp %>% filter(EA_Member),
            aes(fill = FillColor), color = "black", size = 0.2) +
    geom_sf(data = europe_gdp %>% filter(SpecialCountry),
            fill = NA, color = "red", size = 0.9) +
    geom_sf_text(data = europe_gdp %>% filter(SpecialCountry & EA_Member),
                 aes(label = paste0(round(GDP_Share_EA, 1), "%")),
                 size = 5.2, color = "black", fontface = "bold", family = "serif") +
    
    # Colori (senza legenda)
    scale_fill_gradientn(
      colours = c("white", "#d4e8f7", "#a9d0f5", "#77b8eb", "#449fe1", "#1177cc"),
      limits = c(0, 30),
      na.value = "grey80",
      guide = "none"
    ) +
    
    coord_sf(xlim = c(-25, 45), ylim = c(33, 72), expand = FALSE) +
    labs(title = "Euro Area GDP Composition")
  
  print(gdp_plot)
  
  ggsave("Figures/euro_area_gdp_map.pdf", plot = last_plot(), width = 9, height = 7)
  
  # 6. Tabella
  gdp_table <- gdp_data %>%
    filter(Country %in% euro_area_countries) %>%
    arrange(desc(GDP_Share_EA)) %>%
    dplyr::select(Country, `GDP (Million €)` = value, `Contribution (%)` = GDP_Share_EA)
  
  if (show_table) {
    return(kable(gdp_table, format = "latex", booktabs = TRUE, digits = 2,
                 caption = paste("GDP Contribution by Country in the Euro Area –", quarter)) %>%
             kable_styling(latex_options = c("hold_position", "striped", "scale_down")))
  } else {
    return(list(map = gdp_plot, table = gdp_table))
  }
}


################################################################################
########################## DESCRIPTIVE DMFM ANALYSIS ###########################
################################################################################

# ==============================================================================
# FACTORS VISUALIZATION
# ==============================================================================

plot_factors_tensor <- function(Factors, title = "Euro Area Business Cycle Factor ", free_y = TRUE, smooth = FALSE,
                                start_date = as.Date("2000-01-01"), freq = "month") {
  T <- dim(Factors)[1]
  k1 <- dim(Factors)[2]
  k2 <- dim(Factors)[3]
  
  # Costruzione sequenza temporale
  time_index <- seq.Date(from = start_date, by = freq, length.out = T)
  
  # Tabella delle recessioni (puoi anche passarla come argomento opzionale)
  recession_df <- data.frame(
    start = as.Date(c("2008-01-01", "2011-07-01", "2020-01-01")),
    end   = as.Date(c("2009-06-30", "2013-03-31", "2020-06-30"))
  )
  
  # Costruzione DataFrame completo
  df <- data.frame()
  for (i in 1:k1) {
    for (j in 1:k2) {
      df_temp <- data.frame(
        time = 1:T,
        Date = time_index,
        value = Factors[, i, j],
        row_factor = paste0("F_row_", i),
        col_factor = paste0("F_col_", j),
        component = paste0("F(", i, ",", j, ")")
      )
      df <- rbind(df, df_temp)
    }
  }
  
  # Colori personalizzati
  n_components <- length(unique(df$component))
  col_palette <- viridis::viridis(n_components, option = "D")
  
  # Plot con fasce grigie di recessione
  p <- ggplot(df, aes(x = Date, y = value, color = component)) +
    geom_rect(data = recession_df, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.5) +
    {if (smooth) geom_smooth(se = FALSE, linewidth = 0.8) else geom_line(linewidth = 0.8)} +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ component, scales = if (free_y) "free_y" else "fixed") +
    scale_color_manual(values = col_palette) +
    theme_light(base_size = 14) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = title,
      x = "Date",
      y = "Value"
    )
  
  return(p)
}

# ==============================================================================
# LIKELIHOOD CONVERGENCE
# ==============================================================================

plot_em_llk <- function(history_df) {
  # Plots log-likelihood evolution across EM iterations ( Use .
  # Input:
  #   - history_df: data frame with columns `iter` and `llk_new`
  # Output:
  #   - A ggplot object
  
  library(ggplot2)
  
  ggplot(history_df, aes(x = iter, y = llk_new)) +
    geom_line(color = "#0072B2") +
    geom_point(color = "#D55E00") +
    labs(
      title = "Log-likelihood Convergence in EM Alghoritm",
      x = "Iteration",
      y = "Log-Likelihood"
    ) +
    theme_minimal()
}


# ==============================================================================
# EVALUATION OF PERFORMANCES
# ==============================================================================

# nowcast_list <- roll_nowcast$M1
# gdp_idx = res$quarterly_start_idx

compute_rmsfe <- function(nowcast_list, Y_true, gdp_idx, DatesM) {
  # nowcast_list: lista con nomi-date (es. "2024-01-01"), ciascun elemento è un vettore di p1 paesi
  # Y_true: array [T, p1, p2] destandardizzato
  # gdp_idx: indice colonna GDP
  # DatesM: vettore delle date mensili corrispondenti a Y_true (es. res$DatesM)
  
  dates_vec <- as.Date(names(nowcast_list))  # vettore di date
  n <- length(dates_vec)
  p1 <- dim(Y_true)[2]
  rmsfe_mat <- matrix(NA, nrow = p1, ncol = n)
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    # Calcola il mese all’interno del trimestre
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1 # Fine trimestre ( wirdoooo)
    
    if (t_gdp_true > dim(Y_true)[1]) next  # fuori dai dati
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        rmsfe_mat[j, i] <- (forecast - actual)^2
      }
    }
  }
  
  # Media per paese
  rmsfe_vec <- sqrt(rowMeans(rmsfe_mat, na.rm = TRUE))
  names(rmsfe_vec) <- paste0("Country_", 1:p1)
  return(rmsfe_vec)
}

# =============================================================================
# Boxplot RMSFME
# =============================================================================

compute_rmsfe_long <- function(nowcast_list, Y_true, gdp_idx, DatesM, month_label = "M1") {
  dates_vec <- as.Date(names(nowcast_list))
  p1 <- dim(Y_true)[2]
  result_df <- data.frame()
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1
    if (t_gdp_true > dim(Y_true)[1]) next
    
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        err <- forecast - actual
        result_df <- rbind(result_df, data.frame(
          Country = paste0("Country_", j),
          Date = date_t,
          Month = month_label,
          Error = err,
          RMSFE = sqrt(err^2)
        ))
      }
    }
  }
  return(result_df)
}



# =============================================================================
# Plot aggiornamenti del nowcast 
# =============================================================================


plot_nowcast_vs_true_gdp <- function(roll_nowcast, Y_true, gdp_idx, DatesM) {
  p1 <- dim(Y_true)[2]  # numero paesi
  
  # Nomi coerenti: "Country_1", ..., "Country_p1"
  country_names <- paste0("Country_", 1:p1)
  
  # Funzione ausiliaria per costruire il data.frame da ciascuna lista di nowcast
  build_df <- function(nc_list, label) {
    df <- data.frame(
      Date = as.Date(names(nc_list)),
      do.call(rbind, nc_list)
    )
    colnames(df)[-1] <- country_names
    df$Month <- label
    return(df)
  }
  
  # Costruzione tabelle per ciascun mese
  df_M1 <- build_df(roll_nowcast$M1, "M1")
  df_M2 <- build_df(roll_nowcast$M2, "M2")
  df_M3 <- build_df(roll_nowcast$M3, "M3")
  
  # Unione di tutti i nowcast
  df_all <- bind_rows(df_M1, df_M2, df_M3) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "Nowcast") %>%
    mutate(Country = factor(Country, levels = country_names))
  
  # Serie reale del GDP (solo nei mesi in cui esce: ogni 3 mesi)
  df_true <- data.frame(
    Date = DatesM,
    Y_true_gdp = Y_true[ , , gdp_idx]
  )
  colnames(df_true)[-1] <- country_names
  df_true <- df_true %>%
    slice(seq(3, n(), by = 3)) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "True")
  
  # Merge tra nowcast e serie reale
  plot_df <- left_join(df_all, df_true, by = c("Date", "Country"))
  
  # Plot
  ggplot(plot_df, aes(x = Date)) +
    geom_line(aes(y = Nowcast, color = Month), linewidth = 1) +
    geom_point(aes(y = True), color = "black", size = 2, shape = 16, na.rm = TRUE) +
    facet_wrap(~ Country, scales = "free_y") +
    labs(
      title = "Nowcast del GDP vs Valori Reali per Paese",
      x = "Data", y = "GDP",
      color = "Mese del Trimestre"
    ) +
    scale_color_manual(values = c("M1" = "steelblue", "M2" = "orange", "M3" = "darkgreen")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}

# =============================================================================
# Plot Nowcast Pre and Post Covid
# =============================================================================

plot_nowcast_vs_true_gdp_filtered <- function(roll_nowcast, Y_true, gdp_idx, DatesM, 
                                              month_filter = c("M1", "M2", "M3"),
                                              start_date,
                                              end_date) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  p1 <- dim(Y_true)[2]
  country_names <- paste0("Country_", 1:p1)
  
  # Helper per costruire il dataframe da ogni lista di nowcast
  build_df <- function(nc_list, label) {
    df <- data.frame(
      Date = as.Date(names(nc_list)),
      do.call(rbind, nc_list)
    )
    colnames(df)[-1] <- country_names
    df$Month <- label
    return(df)
  }
  
  # Costruzione dei dataframe per ciascun mese richiesto
  df_all <- list()
  if ("M1" %in% month_filter) df_all <- append(df_all, list(build_df(roll_nowcast$M1, "M1")))
  if ("M2" %in% month_filter) df_all <- append(df_all, list(build_df(roll_nowcast$M2, "M2")))
  if ("M3" %in% month_filter) df_all <- append(df_all, list(build_df(roll_nowcast$M3, "M3")))
  
  df_nowcast <- bind_rows(df_all) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "Nowcast") %>%
    mutate(Country = factor(Country, levels = country_names)) %>%
    filter(Date >= start_date, Date <= end_date)
  
  # Serie reale del GDP ogni 3 mesi
  df_true <- data.frame(
    Date = DatesM,
    Y_true_gdp = Y_true[ , , gdp_idx]
  )
  colnames(df_true)[-1] <- country_names
  df_true <- df_true %>%
    slice(seq(3, n(), by = 3)) %>%
    pivot_longer(cols = all_of(country_names), names_to = "Country", values_to = "True")
  
  # Merge nowcast + true GDP
  plot_df <- left_join(df_nowcast, df_true, by = c("Date", "Country"))
  
  # Plot finale
  ggplot(plot_df, aes(x = Date)) +
    geom_line(aes(y = Nowcast, color = Month), linewidth = 1) +
    geom_line(aes(y = True, group = Country), color = "black", linewidth = 0.8, na.rm = TRUE) +
    geom_point(aes(y = True), color = "black", size = 2, shape = 16, na.rm = TRUE) +
    facet_wrap(~ Country, scales = "free_y") +
    labs(
      title = "Nowcast del GDP vs Valori Reali",
      x = "Data", y = "GDP",
      color = "Mese del Trimestre"
    ) +
    scale_color_manual(values = c("M1" = "steelblue", "M2" = "orange", "M3" = "darkgreen")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}


# ============================================================================== 
# RMSFE pre post covid 
# ==============================================================================

compute_rmsfe_period <- function(nowcast_list, Y_true, gdp_idx, DatesM, start_date, end_date) {
  dates_vec <- as.Date(names(nowcast_list))
  p1 <- dim(Y_true)[2]
  rmsfe_mat <- matrix(NA, nrow = p1, ncol = length(dates_vec))
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    if (date_t < start_date || date_t > end_date) next
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1
    if (t_gdp_true > dim(Y_true)[1]) next
    
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        rmsfe_mat[j, i] <- (forecast - actual)^2
      }
    }
  }
  sqrt(rowMeans(rmsfe_mat, na.rm = TRUE))
}

compute_rmsfe_long_period <- function(nowcast_list, Y_true, gdp_idx, DatesM, month_label = "M1",
                                      start_date, end_date) {
  dates_vec <- as.Date(names(nowcast_list))
  p1 <- dim(Y_true)[2]
  df <- data.frame()
  
  for (i in seq_along(dates_vec)) {
    date_t <- dates_vec[i]
    if (date_t < start_date || date_t > end_date) next
    t <- which(DatesM == date_t)
    if (length(t) == 0) next
    
    m_trimestre <- (as.integer(format(date_t, "%m")) - 1) %% 3 + 1
    t_gdp_true <- t + (3 - m_trimestre) + 1
    if (t_gdp_true > dim(Y_true)[1]) next
    
    for (j in 1:p1) {
      forecast <- nowcast_list[[i]][j]
      actual <- Y_true[t_gdp_true, j, gdp_idx]
      if (!is.na(actual) && !is.na(forecast)) {
        err <- forecast - actual
        df <- rbind(df, data.frame(
          Country = paste0("Country_", j),
          Date = date_t,
          Month = month_label,
          Error = err,
          RMSFE = sqrt(err^2)
        ))
      }
    }
  }
  return(df)
}
