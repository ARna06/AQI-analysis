#Main Workflow--modified for KPSS validation

library(vars)
library(urca)
library(tseries)
library(MASS)
library(lubridate)
library(dplyr)
library(tidyverse)
library(igraph)
library(ggraph)
library(grid)
#Workflow---

#Function to load data
setup <- function(Dataset){
  Dataset$time <- as.POSIXct(
    paste(Dataset$date, Dataset$hour),
    format = "%Y-%m-%d %H",
    tz = "Asia/Kolkata"
  )
  
  Dataset <-Dataset %>%
    mutate(
      date = ymd(Dataset$date),
      datetime = date+hours(Dataset$hour)
    )%>%
    arrange(datetime)
  
  df_pm25<<-Dataset %>%
    dplyr::select(ext_temp,pm25,pm10,ext_humi,so2,no2) %>%
    na.omit()
  
  return(df_pm25)
}

#0. Check for %missingness in each column
#0.1 Function to check missingness
check_missingness <- function(file_path,threshold=20) {
  
  # Read CSV
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  count<-0
  # Total number of rows
  n <- nrow(data)
  
  # Compute missing counts and percentages
  missing_count <- sapply(data, function(x) sum(is.na(x)))
  
  missing_percent <- sapply(data, function(x) {
    round(mean(is.na(x)) * 100, 2)
  })
  
  # Create summary table
  result <- data.frame(
    Column = names(data),
    Missing_Count = missing_count,
    Missing_Percentage = missing_percent
  )
  flagged_cols <- names(missing_percent[missing_percent > threshold])
  if(length(flagged_cols)>0){
    count=count+1
    }
  return(count)
  #return(result)
}
##Function to identify missing cols
get_high_missing_cols <- function(
    data,
    threshold = 0.20
){
  
  missing_prop <- colMeans(is.na(data))
  
  names(
    missing_prop[
      missing_prop > threshold
    ]
  )
}
#1. Test for stationarity using ADF
check_stationarity <- function(dataset, param, lag) {
  x <- dataset[[param]]
  adf_result <<- ur.df(x, type = "trend", lags = lag)
  s <- summary(adf_result)
  
  test_statistic <- s@teststat[1]        # ADF tau statistic for drift
  critical_value <- s@cval["tau3", "5pct"]    # 5% critical value
  #cat(test_statistic, critical_value, "\n")
  is_stationary <- test_statistic < critical_value
  
  return(is_stationary)
}

#2. Test whether residuals from ADF regression are white noise
##We use Ljung-Box test
find_lag_threshold <- function(dataset, param){
  pthreshold <- 0.05
  x <- dataset[[param]]
  for (lag in 1:400){
    ### Fitting trend
    adf_temp <- ur.df(x, type = "trend", lags = lag)
    res <<- residuals(adf_temp@testreg)
    
    # This is an unrelated Lag. 
    test_result <<- Box.test(res, lag = 50, type = "Ljung-Box")
    if (test_result$p.value > pthreshold) {
      return(lag)
    }
  }
}

#3. Validate stationarity using KPSS
check_kpss_stationarity <- function(col, alpha = 0.05) {
      # Remove NA values
      series <- na.omit(col)
      # Perform KPSS test
      test <- kpss.test(series)
      # Decision rule
      return(test$p.value>alpha)
}
#4. Difference if necessary
#Function
AttainStationarity <- function(dataset, ParamNames) {
  difference_counts <- rep(0, length(ParamNames))
  names(difference_counts) <- ParamNames
  for (param in ParamNames) {
    #lag <- find_lag_threshold(dataset, param)
    # cat("lag threshold of", param, lag, "\n")
    difference_counter <- 0
    
    while(!check_kpss_stationarity(dataset[[param]])) {
      # cat(param, "is not stationary. Differencing again...\n")
      
      dataset <- dataset %>%
        mutate(!!sym(param) := !!sym(param) - lag(!!sym(param))) %>%
        filter(!is.na(!!sym(param)))
      
      difference_counter <- difference_counter + 1
      #lag <- find_lag_threshold(dataset, param)
      # cat("After differencing", difference_counter, "times, new lag threshold of", param, lag, "\n")
    }
    
    # cat(param, "is stationary after differencing", difference_counter, "times. Final lag threshold:", lag, "\n\n")
    difference_counts[param] <- difference_counter
  }
  
  # Return the joined vector: dataset + difference_counts
  return(list(dataset = dataset, difference_counts = difference_counts))
}

#5. Fit pairwise VAR models and
#6. Get Granger Causality results pairwise for each pair of variables

granger_pair <- function(var1, var2, data, lag.max = 48, ic = "AIC(n)") {
  
  df <- data %>%
    dplyr::select(all_of(c(var1, var2))) %>%
    na.omit()
  
  # Select optimal lag
  p <- VARselect(df, lag.max = lag.max, type = "const")$selection[[ic]]
  
  # Fit VAR
  fit <- VAR(df, p = p, type = "const")
  
  # Granger causality both directions
  gc_1_causes_2 <- causality(fit, cause = var1)$Granger  # var1 → var2
  gc_2_causes_1 <- causality(fit, cause = var2)$Granger  # var2 → var1
  
  tibble(
    cause    = c(var1,  var2),
    effect   = c(var2,  var1),
    opt_lag  = p,
    F_stat   = c(gc_1_causes_2$statistic, gc_2_causes_1$statistic),
    df1      = c(gc_1_causes_2$parameter[1], gc_2_causes_1$parameter[1]),
    df2      = c(gc_1_causes_2$parameter[2], gc_2_causes_1$parameter[2]),
    p_value  = c(gc_1_causes_2$p.value,    gc_2_causes_1$p.value),
    sig      = symnum(
      c(gc_1_causes_2$p.value, gc_2_causes_1$p.value),
      cutpoints = c(0, .001, .01, .05, .1, 1),
      symbols   = c("***", "**", "*", ".", " ")
    ) %>% as.character()
  )
}

##Draw Granger-causality diagrams
plot_granger_graph <- function(
    gc_tbl,
    title = "Granger Causality Network"
) {
  
  
  # -----------------------------------
  # Pretty labels
  # -----------------------------------
  
  label_map <- c(
    "ext_temp" = "ext_temp",
    "ext_humi" = "ext_humi",
    "pm10"     = "PM[10]",
    "pm25"     = "PM[2.5]",
    "so2" = "SO[2]",
    "no2" = "NO[2]"
  )
  
  # -----------------------------------
  # Filter significant edges
  # -----------------------------------
  
  edges <- gc_tbl %>%
    filter(p_value <= 0.10) %>%
    mutate(
      line_type = ifelse(
        p_value < 0.05,
        "solid",
        "dashed"
      ),
      
      edge_width = case_when(
        p_value < 0.001 ~ 2.5,
        p_value < 0.01  ~ 2.0,
        p_value < 0.05  ~ 1.5,
        TRUE            ~ 0.8
      )
    )
  
  # -----------------------------------
  # Create graph
  # -----------------------------------
  
  g <- graph_from_data_frame(
    edges,
    directed = TRUE
  )
  
  # -----------------------------------
  # Apply labels
  # -----------------------------------
  
  V(g)$label <- label_map[V(g)$name]
  
  # -----------------------------------
  # Plot graph
  # -----------------------------------
  
  ggraph(g, layout = "circle") +
    
    geom_edge_arc(
      aes(
        width    = edge_width,
        linetype = line_type
      ),
      
      colour = "black",
      
      arrow = arrow(
        type   = "closed",
        length = unit(2.5, "mm")
      ),
      
      end_cap   = circle(13, "mm"),
      start_cap = circle(13, "mm"),
      
      strength = 0.12,
      
      show.legend = TRUE
    ) +
    
    scale_edge_width_identity() +
    
    scale_edge_linetype_manual(
      name = "Significance",
      
      values = c(
        "solid"  = "solid",
        "dashed" = "dashed"
      ),
      
      labels = c(
        "solid"  = "p < 0.05",
        "dashed" = "0.05 < p < 0.10"
      )
    ) +
    
    # Nodes
    
    geom_node_point(
      size   = 22,
      shape  = 21,
      fill   = "white",
      colour = "black",
      stroke = 1.2
    ) +
    
    # Labels
    
    geom_node_text(
      aes(label = label),
      
      parse = TRUE,
      
      size = 4.2,
      
      fontface = "bold",
      
      colour = "black"
    ) +
    
    labs(title = title) +
    
    theme_void(base_size = 13) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face  = "bold",
        size  = 16,
        margin = margin(b = 10)
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        face = "bold"
      ),
      
      plot.background = element_rect(
        fill   = "white",
        colour = NA
      ),
      
      plot.margin = margin(
        30, 30, 30, 30
      )
    )
}
#7. Fit a multivariate VAR model
fit_multivar_var <- function(data,
                             vars,
                             lag.max = 48,
                             ic = "AIC(n)") {
  
  df <- data %>%
    dplyr::select(all_of(vars)) %>%
    na.omit()
  
  p <- VARselect(
    df,
    lag.max = lag.max,
    type = "const"
  )$selection[[ic]]
  
  fit <- VAR(df,
             p = p,
             type = "const")
  
  return(list(
    fit = fit,
    lag = p
  ))
}

#8. Find FEVD and IRF 
##IRF
run_irf <- function(fit,
                    impulse,
                    response,
                    filename,
                    output_dir,
                    n.ahead = 48,
                    ortho = TRUE,
                    boot = TRUE,
                    ci = 0.95){
  
  irf_res <- irf(
    fit,
    impulse = impulse,
    response = response,
    n.ahead = n.ahead,
    ortho = ortho,
    boot = boot,
    ci = ci,
    runs=1000
  )
  
  # Create output directory if absent
  dir.create(output_dir,
             recursive = TRUE,
             showWarnings = FALSE)
  
  # Save plot
  png(
    filename = paste0(
      output_dir, "/",
      filename, "_",
      impulse, "_to_", response,
      ".png"
    ),
    width = 1200,
    height = 800
  )
  
  plot(irf_res)
  
  dev.off()
  
  return(irf_res)
}

##FEVD
run_fevd <- function(fit,
                     variable,
                     filename,
                     output_dir,
                     n.ahead = 48){
  
  fevd_res <- fevd(
    fit,
    n.ahead = n.ahead
  )
  
  variable_fevd <- as.data.frame(fevd_res[[variable]])
  
  # Create output directory
  dir.create(output_dir,
             recursive = TRUE,
             showWarnings = FALSE)
  
  # Save CSV
  write.csv(
    variable_fevd,
    paste0(
      output_dir, "/",
      filename, "_",
      variable,
      "_fevd.csv"
    ),
    row.names = FALSE
  )
  
  # Save FEVD plot
  png(
    filename = paste0(
      output_dir, "/",
      filename, "_",
      variable,
      "_fevd.png"
    ),
    width = 1200,
    height = 800
  )
  
  plot(fevd_res)
  
  dev.off()
  
  return(variable_fevd)
}

###Bringing all functions together
#Parent directory
dir<-"/Users/rahulkonar/Documents/Stat Comprehensive/AQMS_Data"
# Get all CSV files beginning with "file"
csv_files <- list.files(
  path = dir,
  pattern = "file_221.csv",
  full.names = TRUE
)

#Parameter names
param_names<-c("ext_temp","ext_humi","pm10","pm25","so2","no2")

#Running entire workflow
for(fname in csv_files){
  full_path <- fname #dir+file_path
  #missingness condition
  #if(check_missingness(full_path)==0){
    data<-read.csv(full_path) #read the data if missing values<20%
    data<-setup(data) #data preprocessing for appropriate format
    #lags<-find_lag_threshold(data,param) #finding appropriate lag for ADF
    #cheat code--1 
    ##use kpss only for checking stationarity (in which case I don't need to find lag)
    df_pm25<-AttainStationarity(data,param_names)
    # --- After AttainStationarity ---
    stationary_data <- df_pm25$dataset
    
    # 5 & 6. Pairwise Granger Causality
    pairs <- combn(param_names, 2, simplify = FALSE)
    
    gc_results <- map_dfr(pairs, function(p) {
      tryCatch(
        granger_pair(p[1], p[2], stationary_data, lag.max = 48),
        error = function(e) {
          message("Granger failed for ", p[1], " ~ ", p[2], ": ", e$message)
          NULL
        }
      )
    })
    
    # Save Granger results table
    station_name <- tools::file_path_sans_ext(basename(fname))
    out_dir <- file.path(dir, "outputs", station_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(gc_results,
              file.path(out_dir, paste0(station_name, "_granger.csv")),
              row.names = FALSE)
    
    # Plot and save Granger causality network
    granger_plot <- plot_granger_graph(gc_results)
    ggsave(
      filename = file.path(out_dir, paste0(station_name, "_granger_network.png")),
      plot     = granger_plot,
      width    = 10,
      height   = 8,
      dpi      = 150
    )
  

#     # 7. Fit multivariate VAR
#     var_result <- fit_multivar_var(stationary_data, param_names, lag.max = 48)
#     fit <- var_result$fit
#     message("Optimal lag for multivariate VAR: ", var_result$lag)
#     
#     # 8a. IRF — all impulse → response combinations
#     irf_pairs <- expand.grid(impulse  = param_names,
#                              response = param_names,
#                              stringsAsFactors = FALSE) %>%
#       filter(impulse != response)   # skip self-effects (optional)
#     
#     irf_dir <- file.path(out_dir, "irf")
#     
#     pwalk(irf_pairs, function(impulse, response) {
#       tryCatch(
#         run_irf(
#           fit        = fit,
#           impulse    = impulse,
#           response   = response,
#           filename   = station_name,
#           output_dir = irf_dir,
#           n.ahead    = 48
#         ),
#         error = function(e)
#           message("IRF failed ", impulse, "->", response, ": ", e$message)
#       )
#     })
#     
#     # 8b. FEVD — one plot per response variable
#     fevd_dir <- file.path(out_dir, "fevd")
#     
#     walk(param_names, function(v) {
#       tryCatch(
#         run_fevd(
#           fit        = fit,
#           variable   = v,
#           filename   = station_name,
#           output_dir = fevd_dir,
#           n.ahead    = 48
#         ),
#         error = function(e)
#           message("FEVD failed for ", v, ": ", e$message)
#       )
#     })
#     
#     message("Done: ", station_name)
#     
#   } 
#else {
    #message("Skipping ", basename(fname), " — high missingness in: ",
            #paste(check_missingness(fname), collapse = ", "))
  #}
 }

 

problematic_files <- c(
  "file_121.csv",
  "file_125.csv",
  "file_127.csv",
  "file_137.csv",
  "file_149.csv",
  "file_16.csv",
  "file_161.csv",
  "file_162.csv",
  "file_163.csv",
  "file_165.csv",
  "file_17.csv",
  "file_177.csv",
  "file_181.csv",
  "file_189.csv",
  "file_199.csv",
  "file_202.csv",
  "file_207.csv",
  "file_21.csv",
  "file_211.csv",
  "file_212.csv",
  "file_213.csv",
  "file_214.csv",
  "file_224.csv",
  "file_226.csv",
  "file_228.csv",
  "file_229.csv",
  "file_237.csv",
  "file_239.csv",
  "file_240.csv",
  "file_241.csv",
  "file_248.csv",
  "file_26.csv",
  "file_265.csv",
  "file_287.csv",
  "file_299.csv",
  "file_30.csv",
  "file_306.csv",
  "file_307.csv",
  "file_52.csv",
  "file_53.csv",
  "file_70.csv",
  "file_72.csv",
  "file_92.csv"
)

# -----------------------------------
# Create full paths
# -----------------------------------

problematic_paths <- file.path(
  dir,
  problematic_files
)

for(fname in problematic_paths){
  
  full_path <- fname
  data <- read.csv(fname)
  data <- setup(data)
  # -----------------------------------
  # Identify columns with high missingness
  # -----------------------------------
  high_missing_cols <- get_high_missing_cols(data)
  # -----------------------------------
  # Identify problematic columns
  # -----------------------------------
  if(length(high_missing_cols) > 0){
    
    message(
      "Dropping columns in ",
      basename(fname),
      ": ",
      paste(high_missing_cols,
            collapse = ", ")
    )
    
    data <- data %>%
      select(-all_of(high_missing_cols))
  }
  
  
 
  
  
  # -----------------------------------
  # Determine remaining parameters
  # -----------------------------------
  
  current_params <- intersect(
    param_names,
    colnames(data)
  )
  
  # Need at least 2 variables
  if(length(current_params) < 2){
    
    message(
      "Skipping ",
      basename(fname),
      " — insufficient variables after dropping columns."
    )
    
    next
  }
  
  # -----------------------------------
  # Remove remaining NA rows
  # -----------------------------------
  
  data <- na.omit(data)
  
  # -----------------------------------
  # Stationarity transformation
  # -----------------------------------
  
  df_pm25 <- AttainStationarity(
    data,
    current_params
  )
  
  stationary_data <- df_pm25$dataset
  
  # -----------------------------------
  # Pairwise Granger causality
  # -----------------------------------
  
  pairs <- combn(
    current_params,
    2,
    simplify = FALSE
  )
  
  gc_results <- map_dfr(
    pairs,
    function(p){
      
      tryCatch(
        granger_pair(
          p[1],
          p[2],
          stationary_data,
          lag.max = 48
        ),
        
        error = function(e){
          
          message(
            "Granger failed for ",
            p[1],
            " ~ ",
            p[2],
            ": ",
            e$message
          )
          
          NULL
        }
      )
    }
  )
  
  # -----------------------------------
  # Output directory
  # -----------------------------------
  
  station_name <- tools::file_path_sans_ext(
    basename(fname)
  )
  
  out_dir <- file.path(
    dir,
    "outputs",
    station_name
  )
  
  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # -----------------------------------
  # Save Granger table
  # -----------------------------------
  
  write.csv(
    gc_results,
    
    file.path(
      out_dir,
      paste0(
        station_name,
        "_granger.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  # -----------------------------------
  # Plot Granger graph
  # -----------------------------------
  
  granger_plot <- plot_granger_graph(
    gc_results
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        station_name,
        "_granger_network.png"
      )
    ),
    
    plot = granger_plot,
    
    width = 6,
    height = 5,
    
    dpi = 300
  )
  
  message(
    "Done: ",
    station_name
  )
}
