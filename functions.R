tidy_ts <- function(data, bands, timestamps, coverage, scale_factor = 1) {
  
  fun <- function(x) {
    names(bands) <- bands
    ts_bands <- lapply(bands, function(b){
      x %>% 
        dplyr::select(starts_with(b)) %>% 
        as.numeric() * scale_factor
    }) %>% 
      tibble::as_tibble()
    
    tibble::tibble(Index = timestamps) %>% 
      dplyr::bind_cols(ts_bands) 
  }
  
  data %>% 
    purrrlyr::by_row(..f = fun, .collate = "list") %>% 
    dplyr::transmute(longitude = data$x,
                     latitude = data$y,
                     start_date = timestamps[1],
                     end_date = tail(timestamps, 1),
                     label_id = as.integer(data$label_id),
                     label = data$label,
                     coverage = factor(coverage),
                     time_series = .out)  
}

tibble_to_twdtwTimeSeries <- function (data.tb){
  
  zoo.ls <- data.tb$time_series %>%
    purrr::map (function (ts) {
      df <- data.frame (ts)
      return (zoo::zoo (df[,2:ncol(df),drop=FALSE], df[,1]))
    })
  
  labels.fc <-  as.factor (data.frame (dplyr::select(data.tb, label))[,1])
  
  ts.tw <-  methods::new("twdtwTimeSeries", timeseries = zoo.ls, labels = labels.fc)
  return (ts.tw)
}

sits_fromTWDTW_time_series <- function (patterns, coverage){
  # get the time series from the patterns
  tb.lst <- purrr::map2 (patterns@timeseries, patterns@labels, function (ts, lab) {
    # tranform the time series into a row of a sits table
    ts.tb <- zoo::fortify.zoo(ts)
    # store the sits table in a list
    mylist        <- list()
    mylist [[1]]  <- tibble::as_tibble (ts.tb)
    # add the row to the sits table
    row   <- tibble::tibble(longitude    = 0.00,
                            latitude     = 0.00,
                            start_date   = ts.tb[1,"Index"],
                            end_date     = ts.tb[nrow(ts.tb),"Index"],
                            label        = as.character (lab),
                            coverage     = coverage,
                            time_series  = mylist)
    return (row)
  })
  # create a sits table to store the result
  patterns.tb <- sits_table()
  patterns.tb <- tb.lst %>%
    purrr::map_df (function (row) {
      dplyr::bind_rows (patterns.tb, row)
    })
  return (patterns.tb)
}


twdtw_parallel <- function(x, y, resample = TRUE, length = NULL, weight.fun = NULL,
                           dist.method = "Euclidean", step.matrix = symmetric1, n = NULL,
                           span = NULL, min.length = 0, theta = 0.5, from = NULL, to = NULL, 
                           by = NULL, breaks = NULL, overlap = 0.5, cl = detectCores(), 
                           chunk_size, rep_chunk, ...){
  
  out <- list(classification = list(), proc_time = tibble::tribble(~n, ~time, ~cores))
  
  
  fun_twdtw <- function(data, y, resample, length, weight.fun,
                        dist.method, step.matrix, n,
                        span, min.length, theta, from, to, by, overlap) {

    zoo.ls <- data$time_series %>%
      purrr::map (function (ts) {
        df <- data.frame (ts)
        return (zoo::zoo (df[,2:ncol(df),drop=FALSE], df[,1]))
      })
    
    ts <- methods::new("twdtwTimeSeries", timeseries = zoo.ls, labels = as.factor(data$label))
    
    res <- dtwSat::twdtwApply(x = ts, y = y, resample = resample, length = length,
                              weight.fun = weight.fun, dist.method = dist.method, 
                              step.matrix = step.matrix, n = n, span = span, 
                              min.length = min.length, theta = theta, keep = TRUE)
    # plotAlignments(res)
    res <- dtwSat::twdtwClassify(res, from = from, to = to, by = by, overlap = 0)
    
    bind_rows(lapply(as.list(res), function(d) {
      d = d[[1]]
      if(nrow(d) != 1){
        return(tibble::tibble(twdtw_distance = Inf, twdtw_label = labels(y)[1])) 
      }
      tibble::tibble(twdtw_distance = d$distance, twdtw_label = as.character(d$label))
    })) %>% bind_cols(data, .)
    
  }
  
  # Step 0: Create Clusters
  cluster <- multidplyr::create_cluster(cores = cl)
  
  for(i in seq_along(rep_chunk)){
    
    k  <- rep_chunk[i]
    
    chunk <- chunk_size$start[k]:chunk_size$end[k]
    
    cat("\n Processing repetition ",i,"/",length(rep_chunk),". Chunk size ", length(chunk))
    
    # Step 1: Setup Clusters dependencies
    multidplyr::cluster_library(cluster, "tidyverse")
    multidplyr::cluster_library(cluster, "purrr") 
    multidplyr::cluster_library(cluster, "dtwSat") 
    multidplyr::cluster_assign_value(cluster, "fun_twdtw", fun_twdtw) 
    multidplyr::cluster_assign_value(cluster, "y", y) 
    multidplyr::cluster_assign_value(cluster, "resample", resample) 
    multidplyr::cluster_assign_value(cluster, "length", length) 
    multidplyr::cluster_assign_value(cluster, "weight.fun", weight.fun) 
    multidplyr::cluster_assign_value(cluster, "dist.method", dist.method)  
    multidplyr::cluster_assign_value(cluster, "step.matrix", step.matrix) 
    multidplyr::cluster_assign_value(cluster, "n", n) 
    multidplyr::cluster_assign_value(cluster, "span", span) 
    multidplyr::cluster_assign_value(cluster, "min.length", min.length) 
    multidplyr::cluster_assign_value(cluster, "theta", theta) 
    multidplyr::cluster_assign_value(cluster, "from", from)  
    multidplyr::cluster_assign_value(cluster, "to", to) 
    multidplyr::cluster_assign_value(cluster, "by", by) 
    multidplyr::cluster_assign_value(cluster, "overlap", overlap)  
    
    # Step 2: Run Parallelized Code 
    start <- proc.time()
    res <- x %>% 
      dplyr::slice(chunk) %>%
      multidplyr::partition(cluster = cluster) %>%
      dplyr::do(fun_twdtw(., y, resample, length, weight.fun, dist.method, 
                          step.matrix, n, span, min.length, theta, from, to, 
                          by, overlap)) %>% 
      dplyr::collect() %>%
      dplyr::ungroup() 
    time_elapsed_parallel <- proc.time() - start 
    
    out$classification[[i]] <- res
    out$proc_time <- bind_rows(out$proc_time, 
                               tibble::tibble(n = length(chunk), 
                                              time = time_elapsed_parallel[3], 
                                              cores = detectCores()))
    
    gc()
  }
  
  # Step 3: Stop and clean clluster
  parallel::stopCluster(cl = cluster)
  parallel::setDefaultCluster(cl = cluster)
  
  out
  
}

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(coef(x)[2], digits = 2),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));                 
}

get_py_F1 <- function(true, pred){
  rPython::python.assign("y_true", true)
  rPython::python.assign("y_pred", pred)
  rPython::python.exec("from sklearn.metrics import f1_score")
  rPython::python.exec("f1 = f1_score(y_true, y_pred, average='weighted')")
  rPython::python.get("f1")  
}

get_clusters <- function(x, labels = unique(x$label), nsamples = 100, path = "./",
                         bands = names(x$time_series[[1]]), type = "hierarchical", 
                         k = 2, distance = "DTW", seed = NULL, ...){
  
  proc_time = tibble::tribble(~n, ~time, ~cores)
  
  if(!is.null(seed))
    set.seed(seed)
  
  cl_fun <- function(data, bands, type, k, distance, ...){
    
    ts <- lapply(data$time_series, function(x) as.matrix(x[, bands[-1], drop = FALSE]) )
    
    clusters  <- dtwclust::tsclust(series = ts, type = type, k = k, distance = distance, ...)
    
    clusters
  }
  
  # Create parallel workers
  cl <- parallel::makeCluster(parallel::detectCores())
  invisible(parallel::clusterEvalQ(cl, library(dtwclust)))
  doParallel::registerDoParallel(cl)
  
  pb <- progress::progress_bar$new(total = length(labels))
  
  for(i in seq_along(labels)){
    start <- proc.time()
    clust_dtw <- training_ts %>% 
      filter(label == labels[i]) %>% 
      dplyr::sample_n(size = nsamples) %>% 
      cl_fun(., bands, type, k, distance, ...)
    time_elapsed_parallel <- proc.time() - start 
    
    saveRDS(clust_dtw, file = paste0(path,"/cluster_",nsamples,"_",distance,"_",labels[i],".rds"))
    
    proc_time <- dplyr::bind_rows(proc_time, 
                                  tibble::tibble(n = length(clust_dtw@cluster), 
                                                 time = time_elapsed_parallel[3], 
                                                 cores = detectCores()))
    
    rm("clust_dtw")
    gc()
    
    pb$tick()
    
  }
  
  # Stop parallel workers
  parallel::stopCluster(cl)
  
  # Return to sequential computations. This MUST be done if stopCluster() was called
  foreach::registerDoSEQ()
  
  proc_time
  
}


filter_ts <- function(tbl, width = 3, max_diff = 0.4){
  
  if(width < 2)
    stop("width smaller than 2")
  
  # get NDVI diff   
  df_ts <- tbl %>% 
    .$NDVI %>% 
    diff()
  
  # Check if there are consecutive values that are equal within a window 
  tag1 <- any(rollapply(df_ts, width = width - 1, FUN = function(x) all(x == 0)))
  
  # Check if there are consecutive values with large differences 
  tag2 <- any(df_ts > max_diff)
  
  return(tag = ifelse(tag1 | tag2, FALSE, TRUE))
  
}

read_clusters_to_twdtw <- function(path, pattern, k = 1){
  filenames <- dir(path = path, pattern = pattern, full.names = TRUE)
  labels_names <- gsub(".rds", "", gsub(pattern, "", basename(filenames)))
  names(filenames) <- labels_names
  clust_dtw <- lapply(filenames, read_rds)
  
  ts_list <- do.call("twdtwTimeSeries", do.call("c", lapply(labels_names, function(i){
    # Select three largest clusters 
    cl = order(clust_dtw[[i]]@clusinfo$size, decreasing = TRUE)[1:k]
    lapply(cl, function(j){
      I <- which(clust_dtw[[i]]@cluster %in% j)
      ts_list <- lapply(clust_dtw[[i]]@datalist[I], function(x) zoo::zoo(x, timestamps)) 
      labels <- rep(paste0(i, "_", j), length(ts_list))
      new("twdtwTimeSeries", timeseries = ts_list, labels = labels)
    })
  })))
  
  ts_list
  
}


get_labels_cluster_classification <- function(x, samples_label){
  x$classification[[1]] %>%
           dplyr::mutate(twdtw_label = gsub("_.", "", twdtw_label)) %>% 
           dplyr::left_join(samples_label, by = c("twdtw_label" = "label")) %>% 
           dplyr::left_join(samples_label, by = c("label" = "label"), suffix = c("_pred", "_true")) %>% 
           dplyr::select(label, twdtw_label, label_id_pred, label_id_true)
}

