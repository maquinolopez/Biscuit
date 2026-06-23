###########################################
###    Setting chronologies using 
###   tie points using iterations        
###    By Marco A. Aquino-Lopez         
###  BISCUIT: 
### Bayesian Integration of Sequences 
### Using Chronology Intersection Tie-points
###########################################


# Function which takes the depth vector and matrix of MCMC results
density_plot <- function(depths, tar_mt, xlabel, ylabel = "proxy units",
                         add = FALSE, axis = TRUE, flip = FALSE,
                         breaks = 150, range_probs = c(0.001, 0.999),
                         normalize = c("column", "global"),
                         max_alpha = 0.75, density_power = 0.7,
                         summary_probs = c(0.05, 0.5, 0.95),
                         show_summary = TRUE,
                         line_col = "black",
                         interval_col = rgb(0.75, 0, 0, 0.85)) {
  normalize <- match.arg(normalize)
  tar_mt <- as.matrix(tar_mt)
  
  if (length(depths) != ncol(tar_mt)) {
    stop("'depths' length must match the number of columns in 'tar_mt'.")
  }
  
  finite_values <- tar_mt[is.finite(tar_mt)]
  if (length(finite_values) == 0) {
    stop("'tar_mt' has no finite values to plot.")
  }
  
  value_limits <- as.numeric(quantile(finite_values,
                                      probs = range_probs,
                                      na.rm = TRUE,
                                      names = FALSE))
  if (!all(is.finite(value_limits)) || diff(value_limits) <= 0) {
    value_limits <- range(finite_values, na.rm = TRUE)
  }
  
  if (!all(is.finite(value_limits)) || diff(value_limits) <= 0) {
    padding <- max(abs(value_limits[1]), 1) * 0.01
    value_limits <- c(value_limits[1] - padding, value_limits[1] + padding)
  }
  
  if (length(breaks) == 1) {
    value_breaks <- seq(value_limits[1], value_limits[2], length.out = breaks + 1)
  } else {
    value_breaks <- sort(unique(breaks))
  }
  
  if (length(value_breaks) < 2) {
    stop("'breaks' must define at least two bin edges.")
  }
  
  if (length(depths) == 1) {
    depth_width <- 1
    depth_edges <- c(depths[1] - depth_width / 2, depths[1] + depth_width / 2)
  } else {
    mids <- depths[-length(depths)] + diff(depths) / 2
    depth_edges <- c(depths[1] - diff(depths)[1] / 2,
                     mids,
                     tail(depths, 1) + tail(diff(depths), 1) / 2)
  }
  
  hist_counts <- matrix(0, nrow = length(value_breaks) - 1, ncol = ncol(tar_mt))
  for (i in seq_len(ncol(tar_mt))) {
    values <- tar_mt[, i]
    values <- values[is.finite(values)]
    values <- values[values >= value_breaks[1] & values <= tail(value_breaks, 1)]
    if (length(values) > 0) {
      hist_counts[, i] <- hist(values, breaks = value_breaks, plot = FALSE)$counts
    }
  }
  
  if (normalize == "column") {
    hist_scale <- apply(hist_counts, 2, max)
    hist_scale[hist_scale == 0] <- NA_real_
    density_level <- sweep(hist_counts, 2, hist_scale, "/")
  } else {
    hist_scale <- max(hist_counts)
    density_level <- hist_counts / hist_scale
  }
  density_level[!is.finite(density_level)] <- 0
  density_level <- density_level ^ density_power
  
  main_line <- colMeans(tar_mt, na.rm = TRUE)
  summary_lines <- apply(tar_mt, 2, quantile,
                         probs = summary_probs,
                         na.rm = TRUE,
                         names = FALSE)
  
  if (!add) {
    if (flip) {
      plot(NA,
           xlim = value_limits,
           ylim = range(depth_edges),
           xlab = ylabel,
           ylab = xlabel,
           xaxt = ifelse(axis, "s", "n"),
           yaxt = ifelse(axis, "s", "n"))
    } else {
      plot(NA,
           xlim = range(depth_edges),
           ylim = value_limits,
           xlab = xlabel,
           ylab = ylabel,
           xaxt = ifelse(axis, "s", "n"),
           yaxt = ifelse(axis, "s", "n"))
    }
  }
  
  for (i in seq_len(ncol(tar_mt))) {
    alpha <- max_alpha * density_level[, i]
    keep <- alpha > 0
    if (!any(keep)) {
      next
    }
    
    cols <- gray(1 - density_level[keep, i], alpha = alpha[keep])
    if (flip) {
      rect(xleft = value_breaks[-length(value_breaks)][keep],
           xright = value_breaks[-1][keep],
           ybottom = depth_edges[i],
           ytop = depth_edges[i + 1],
           col = cols,
           border = NA)
    } else {
      rect(xleft = depth_edges[i],
           xright = depth_edges[i + 1],
           ybottom = value_breaks[-length(value_breaks)][keep],
           ytop = value_breaks[-1][keep],
           col = cols,
           border = NA)
    }
  }
  
  if (show_summary) {
    median_index <- which.min(abs(summary_probs - 0.5))
    lower_index <- 1
    upper_index <- length(summary_probs)
    
    if (flip) {
      lines(summary_lines[lower_index, ], depths, col = interval_col, lwd = 1.2)
      lines(summary_lines[upper_index, ], depths, col = interval_col, lwd = 1.2)
      lines(summary_lines[median_index, ], depths, col = line_col, lwd = 1.6)
    } else {
      lines(depths, summary_lines[lower_index, ], col = interval_col, lwd = 1.2)
      lines(depths, summary_lines[upper_index, ], col = interval_col, lwd = 1.2)
      lines(depths, summary_lines[median_index, ], col = line_col, lwd = 1.6)
    }
  } else {
    if (flip) {
      lines(main_line, depths, col = line_col, lwd = 1.4)
    } else {
      lines(depths, main_line, col = line_col, lwd = 1.4)
    }
  }
}

# target function as preparation
target_density <-function(tar_ages,tar){
  tar = t(tar)
  # Initialize an empty list to store the density objects for each depth
  kde_list <- c()
  
  # Loop through each column (depth) to calculate the kernel density
  for (col in 1:length(tar_ages)) {
    values_at_depth <- tar[,col]
    # Compute kernel density and add to the list
    kde <- density(values_at_depth, kernel = 'gaussian' )
    kde_list <- c(kde_list , as.numeric(kde['bw'] ) )
    # print(kde$bw)
  }

  return(kde_list)
}



Biscuit <- function(Input, Target, folder ='~/Documents/Biscuit/',
                   n_sections = TRUE,n_tie_points, sampling_date,
                   thin=25,burn=1e+3,iters=2.5e+3,
                   shape_acc = 1.5,  mean_acc = 10,
                   strength_mem = 10, mean_mem = .5,
                   d_by=1, 
                   target_age_model = "Bacon", # "Bacon" or "Plum"
                   run_target = FALSE,
                   target_density_method = c("KDE", "GMM"),
                   gmm_function_url = "https://raw.githubusercontent.com/maquinolopez/MISO_Mixture_Inference/main/GMM_funR.R",
                   gmm_max_components = 20,
                   gmm_model_names = "V",
                   gmm_verbose = FALSE,
                   bacon_ask = FALSE,
                   bacon_suggest = FALSE
                   ){
  target_density_method <- match.arg(toupper(target_density_method), c("KDE", "GMM"))
  source(paste0(folder,'/twalk.R'))
  
  # Prepare data 
  colnames <- c("depth", "value")
  # Read CSV into input as dataframe
  input_p  <- read.csv(paste(folder,Input,'-',Target,'/',Input,"_proxy.csv",sep='') , header = TRUE, col.names = colnames)
  target_p <- read.csv(paste(folder,Input,'-',Target,'/',Target,"_proxy.csv",sep='') , header = TRUE, col.names = colnames)
  
  # check number of tie-points 
  # if n_sections is true, the number of sections will be four
  # else, the number of sections will be equal to n_tie_points
  if (n_sections){
    n_sections = as.integer(n_tie_points * 4 )
  }

  # Load run and return the ages at the depths of the tie points
  # it return a matrix. 
  get_ages <- function(tiepoints,target) {
    if (target_age_model == "Bacon"){
    library(rbacon)
    Bacon(target,run=run_target,coredir=paste0(folder,Input,'-',Target,'/'),ask = bacon_ask, suggest = bacon_suggest)
    agedepth() 
    
    depths <- tiepoints$tar_tie
    df <- matrix(NA, nrow = length(depths), ncol = length(Bacon.Age.d(depths[1])))
    cat('\n Getting the age posterior distributions')
    pb <- txtProgressBar(min = 0, max = length(depths), style = 3, char = "=") 
    for (i in 1:length(depths)) {
      df[i, ] <- as.numeric(Bacon.Age.d(depths[i]))
      setTxtProgressBar(pb, i)
    }
    close(pb)  # Close the progress bar
    return(df)
    }
    if (target_age_model == "Plum"){
    library(rplum)
    Plum(target,run=run_target,coredir=paste0(folder,Input,'-',Target,'/'),ask = bacon_ask, suggest = bacon_suggest)
    agedepth() 
    
    depths <- tiepoints$tar_tie
    df <- matrix(NA, nrow = length(depths), ncol = length(Bacon.Age.d(depths[1])))
    cat('\n Getting the age posterior distributions')
    pb <- txtProgressBar(min = 0, max = length(depths), style = 3, char = "=") 
    for (i in 1:length(depths)) {
      df[i, ] <- as.numeric(Bacon.Age.d(depths[i]))
      setTxtProgressBar(pb, i)
    }
    close(pb)  # Close the progress bar
    return(df)
    }
  }


  # Function to select tie points by user
  select_tiepoints <- function(input, target, n) {
    n = 2*n
    par(mfrow = c(1, 1))  # Set up a single plot
    
    # Calculate the y-offset for the target time series
    y_offset <- max(input[,2]) +  min(target[2]) + 0.2  # Adjust the offset for better visualization
    
    # Calculate the y-axis limits that encompass both series
    y_min <- min(min(input[,2]), min(target[,2] + y_offset))
    y_max <- max(max(input[,2]), max(target[,2] + y_offset))
    target[,2] = y_offset + target[,2]
    # Plot both time series with adjusted y-axis limits
    plot(input, type = "l", main = "Time Series Comparison", xlab = "Time", ylab = "Value", col = "blue", ylim = c(y_min-1, y_max+1), yaxt = "n")
    lines(target, col = "red")
    
    legend("topleft", legend = c("Input", "Target"),  col = c(rgb(0,0,1,.9), rgb(1,0,0,.9)),  lty = 1, lwd = 2, bty = "n")
    
    tie_input <- c()
    tie_target <- c()
    for (i in 1:n) {
      if (i %% 2 == 1) {
        cat(paste("Select tiepoint", (i+1)/2, " in the input (blue): Click on the plot to choose a point.\n"))
        point <- locator(n = 1, type = "p", col = "green", pch = 16)
        ik = which( input[,1]>=point$x )[1]
        tie_input <- c(tie_input , input[ik,1] )    
      }else{
        cat(paste("Select tiepoint", (i)/2, " in the target (red): Click on the plot to choose a point.\n"))
        point <- locator(n = 1, type = "p", col = "green", pch = 16)
        ik = which( target[,1]>=point$x )[1]
        tie_target <- c(tie_target , target[ik,1] )    
      }
      
    }
    tiepoints = data.frame(inp_tie = tie_input, tar_tie = tie_target)
    return(tiepoints)
  }
  
  tie_points <-  select_tiepoints(input_p, target_p, n_tie_points)
  
  # Directory where you want tox search for the file
  directory_path <- paste0(folder,Input,'-',Target,'/',Target)
  matching_files <- list.files(path = directory_path, pattern = paste0(Target,'.csv$'), full.names = TRUE)
  tar_input_file <- read.table(matching_files[1], header = TRUE, sep = ",")  # load input file

  if (length(colnames(tar_input_file))<=5){
    if ( max(tar_input_file$depth) < max(input_p$depth)){
      cat('\nThe target chronology is smaller than the proxy record.')
      cat('\nWe can only align in the windows of the target chronology.')

      }
  }else{
    if ( max(tar_input_file$depth.cm.) < max(input_p$depth) ){
      cat('\nThe target chronology is smaller than the proxy record.')
      cat('\nWe can only align in the windows of the target chronology.')
      break
      }
  }
  

  tar <- get_ages(tie_points ,Target)

  # Find files matching the pattern "_ages.txt" in the directory
  matching_files <- list.files(path = directory_path, pattern = "_ages.txt$", full.names = TRUE)

  if (length(matching_files) > 0) {
    file_path <- matching_files[1]
    target_gdm <- read.table(file_path, header = TRUE, sep = "\t")  
  } else {
    cat("Age depth model for the target not found.\n")
    break
  }
  
  target_p$ages_tar <- approx(x = target_gdm$depth,y = target_gdm$mean,xout = target_p$depth)$y

  if (target_density_method == "KDE") {
    kde_list <- target_density(tie_points$tar_tie,tar)
  } else {
    source(gmm_function_url)
    gmm_fit <- build_best_gmm(t(tar),
                              maxComponents = gmm_max_components,
                              modelNames = gmm_model_names,
                              verbose = gmm_verbose)
    gmm_weights <- lapply(gmm_fit$models, function(model) {
      weights <- as.numeric(model$parameters$pro)
      if (length(weights) == 0) {
        weights <- rep(1 / length(model$parameters$mean), length(model$parameters$mean))
      }
      weights
    })
  }
  # create initial sections
  x_sections <- seq(input_p$depth[1],tail(input_p$depth,1),length.out=n_sections)
  # if surface was taken
  x_sections <- seq(0,tail(input_p$depth,1),length.out=n_sections)
  n_slopes = length(x_sections) - 1
  b_length = x_sections[2]-x_sections[1]
  sampling_date_sd = .1
  
  # calculate shapes 
  # parameter for prior of alphas 
  scale_acc = mean_acc / shape_acc
  
  # parameters for memory 
  m_alpha = strength_mem * mean_mem 
  m_beta =  strength_mem * (1 - mean_mem) 
  

  #####################
  ########## Functions ##########
  #####################
  # support function
  # Calcualtes the lower and upper bounds of the parameter state
  sampling_date = 1950 - sampling_date
  low <- c(sampling_date - 3*sampling_date_sd,
           0, # mem
           rep(0, n_sections-1)) # m_s
  up <- c( sampling_date + 3*sampling_date_sd,  #tao0
           1, # mem,
           rep(Inf, n_sections-1)) # alphas
  mx_age <- 100000
  
  # Function which checks if the parameter state is within the range
  # of the parameter state.
  supp <- function(params){
    # tau0 <- param[1]
    # w <- param[2]
    alp <- alphas(params,n_slopes)
    ifelse(all( c(params > low, params < up) ) & 
             all( alp > 0 ), 
           return(TRUE), 
           return(FALSE)
    )
  }
  
  # Function which given the parameter state and a set of depths returns the estimated ages
  # 
  ttime <- function(x, param,b_length){
    t0 <- param[1]
    # w <- param[2]
    slopes <- param[-c(1,2)]
    t_is <- c(t0, t0 + cumsum(slopes * b_length))
    return( approx(x_sections,t_is,x)$y )
  }
  
  
  # Function which generates from the slopes and w (memory),
  # This is necesary because the alphas are distributed as gamma distribution
  # and MCMC estimates the slopes and w instead of alphas.
  alphas <- function(param, n_slopes){
    w <- param[2]
    slopes <- param[-c(1,2)]
    alf <- ( slopes[-n_slopes] - w * slopes[-1] ) / (1 - w)
    return(c(alf, tail(slopes,1)) )
  }
  
  
  # Prior densities
  logprior <-  function(param){
    tao0 <- param[1]
    w <- param[2]
    # slopes <- param[-c(1,2)]
    T0 <- dnorm(tao0, mean = sampling_date, sd = sampling_date_sd,log = TRUE)
    # memory value
    dmem <- dbeta(w, m_alpha, m_beta, log = TRUE)
    # denisity for acc rates 
    dalphas <- dgamma(alphas(param,n_slopes), shape = shape_acc, scale = scale_acc, log = TRUE)
    #
    d <- c(T0, dmem, dalphas) 
    return(sum(d))
  }
  
  log_density_floor <- log(.Machine$double.xmin)
  
  log_sum_exp <- function(x, floor_value = log_density_floor) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(floor_value)
    }
    max_x <- max(x)
    value <- max_x + log(sum(exp(x - max_x)))
    ifelse(is.finite(value), value, floor_value)
  }
  
  if (target_density_method == "KDE") {
    n_kde <- length(tar[1,])
    sd_convertor <- 1/(1.06*n_kde^(-1/5))

    l_target_kernelR <- function(d, new_x, kde) {
      sd_kernel <- sd_convertor * kde
      target_values <- tar[d,]
      keep <- is.finite(target_values) & is.finite(new_x) &
        is.finite(sd_kernel) & sd_kernel > 0
      if (!any(keep)) {
        return(log_density_floor)
      }
      log_sum_exp(dnorm(x = new_x,
                        mean = target_values[keep],
                        sd = sd_kernel,
                        log = TRUE))
    }
    
    loglikelihood_uqR <- function(params){
      t_times <- ttime(tie_points$inp_tie,params,b_length)
      sum(vapply(seq_along(t_times), function(k) {
        l_target_kernelR(d = k, new_x = t_times[k], kde = kde_list[[k]])
      }, numeric(1)))
    }
  } else {
    l_target_gmmR <- function(new_x, means, variances, weights) {
      sds <- sqrt(variances)
      keep <- is.finite(means) & is.finite(sds) & sds > 0 &
        is.finite(weights) & weights > 0
      if (!any(keep)) {
        return(log_density_floor)
      }
      log_sum_exp(log(weights[keep]) + dnorm(new_x,
                                             mean = means[keep],
                                             sd = sds[keep],
                                             log = TRUE))
    }
    
    loglikelihood_uqR <- function(params){
      t_times <- ttime(tie_points$inp_tie,params,b_length)
      sum(vapply(seq_along(t_times), function(k) {
        l_target_gmmR(new_x = t_times[k],
                      means = gmm_fit$means[[k]],
                      variances = gmm_fit$variances[[k]],
                      weights = gmm_weights[[k]])
      }, numeric(1)))
    }
  }
  
  
  
  # objective function
  obj <- function(param){ 
    - (logprior(param) + loglikelihood_uqR(param)  ) 
    }

  sampler <-  function(rerun=FALSE){
    # generate to0
    t0 <- runif(1, low[1],up[1]) 
    # memory value
    w <-  runif(n = 1, min = .3, max = 1)
    years_sam <- tar[,sample(ncol(tar),1)]
    if (length(tie_points$inp_tie)>1){
      tmp_mean <- lm(years_sam ~  tie_points$inp_tie)$coefficients[2]  
    }else{
      tmp_mean <- (years_sam - sampling_date) /  tie_points$inp_tie
    }
    
    
    ms <-  rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )  
    # rgamma(1, scale = m_slopes/100 , shape = 100 )
    for (i in 1:(n_slopes-1)){
      alpha =  rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )  
      #rgamma(1,scale = m_slopes/100, shape = 100 )
      ms <- c(ms, w * ms[1] + (1-w) * alpha )
    }
    return(c(t0, w, ms))
  }
  
  
  initial_search <- function(){
    x1 <- sampler()
    while(!(supp(x1))){
      x1 <- sampler(rerun = TRUE)  
    }
    return(x1)
  }
  
  
  #### run twalk ####
  start.time <- Sys.time()
  message("\nSearching for initial values...")
  x1 <- initial_search()
  x2 <- initial_search()
  message("\nInitiating the t-walk process...\n")
  # burn = burn * length(x1)
  
  
  output <- Runtwalk(iters,Obj = obj,dim = length(x1),x0 = x1,xp0 = x2,Supp =supp,
                      thinning = length(x1)*thin,burnin= length(x1) * burn )
                     #thinning = thin,burnin= burn )
  
  end_time <- Sys.time()
  total_time = end_time - start.time 
  print(total_time)
  
  # twalk state
  lastx = tail(output$output, 1)
  lastxp = tail(output$outputp, 1)
  # Write the matrix to a CSV file
  # write.csv(matrix_to_save <- matrix(c(lastx, lastxp), nrow = 2,byrow = TRUE), paste0(folder,"twalkstate.csv"), row.names=FALSE)#, col.names=FALSE)

  c <- output$output[-1,]
  energy <- output$Us[-1]
  energy2 <- output$Ups[-1]
  # write.csv(c,paste0(folder,Input,'-',Target,'_',n_sections,".out"), row.names=FALSE)#, col.names=FALSE)
  
  iat <- as.integer(IAT(output,to=output$Tr))
  
  # 

  if (!(iat < 2) ){
    cat("\n================== IAT DIAGNOSTIC ==================\n")
    cat("IAT Value:", iat, "\n")
    cat("Interpretation: The chain exhibits high correlation among successive samples.\n")
    cat("Recommendation: Consider increasing the thinning value and rerunning the chain.\n")
    cat("\n====================================================\n")

  }
  
  ######## Plot results######
  ##
  layout(matrix(c(1,1,1,2,2,2,3,3,3,
                  1,1,1,2,2,2,3,3,3,
                  1,1,1,2,2,2,3,3,3,
                  4,4,4,4,4,4,4,4,4,
                  4,4,4,4,4,4,4,4,4,
                  4,4,4,4,4,4,4,4,4,
                  4,4,4,4,4,4,4,4,4,
                  4,4,4,4,4,4,4,4,4
                  ), 8, 9, byrow = TRUE))
  # plot energy
  plot(energy,type = 'l',col="grey40")
  # lines(energy2,col="grey60")
  
  
  ## Plot posterior and prior of accumulations
  { 
    d_m <- density(output$output[,-c(1,2)])
    # Get y-values of the gamma curve across the specified x-range
    y_gamma <- dgamma(seq(from=0, to=mean_acc +  3 * (shape_acc*scale_acc), length.out=1000),
                      shape = shape_acc, scale = scale_acc)
    
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_gamma),max(d_m$y))) )
    
    # Plot gamma curve with adjusted y-axis limits
    curve(dgamma(x, shape = shape_acc, scale = scale_acc), 
          from=0, to=mean_acc +  3*(shape_acc*scale_acc),  yaxt="n",
          xlab="Accumulation Rates (yr/cm)", ylab="Density", 
          main=" ", ylim=y_range)
    
    # Add density plot
    lines(d_m, col='blue')
  }
  
  ## plot memory prior vs posteriors
  {
    d_w <- density(output$output[,2])
    
    y_beta <- dbeta(seq(from=0, to=1, length.out=1000),
                    m_alpha, m_beta)
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_beta),max(d_w$y))) )
    
    curve(dbeta(x, m_alpha, m_beta), from=0, to=1, ylim=y_range, 
          xlab="memory", ylab="Density", main="",yaxt="n")
    lines(d_w,col='blue')
  }
  
  # Initialize matrix to store assembles
  xs <- seq(input_p[1,1],tail(input_p[,1],1),length.out = 1000)
  t_mat <- matrix(NA, nrow = nrow(c), ncol = length(xs))

  # Loop through each row of c
  for(i in 1:nrow(c)) {
    # Apply tau to breaks with this row of c as params
    t_row <- ttime(xs, c[i,],b_length )
    # Store result in matrix
    t_mat[i,] <- t_row
  }
  
  t_mat2 <- matrix(NA, nrow = nrow(c), ncol = length(x_sections))
  
  # Loop through each row of c
  for(i in 1:nrow(c)) {
    # Apply tau to breaks with this row of c as params
    t_row <- ttime(x_sections, c[i,],b_length )
    # Store result in matrix
    t_mat2[i,] <- t_row
  }
  
  output$age_models = t_mat2
  output$elbows = x_sections
  output$age_means = colMeans(t_mat2)
  output$tie_points
  
  # Quantile functions  
  q5 <- function(x) quantile(x, probs = 0.05)
  q50 <- function(x) quantile(x, probs = 0.50) 
  q95 <- function(x) quantile(x, probs = 0.95)
  
  # Apply to each column
  output$age_min <- apply(t_mat2, 2, q5)
  output$age_medium <- apply(t_mat2, 2, q50)
  output$age_max <- apply(t_mat2, 2, q95)
  
  
  
  
  # plot the age-depth models 
  {
      density_plot(xs,tar_mt = t_mat ,
                   xlabel = 'Depth',ylabel ='Age',flip = F) 
      # Get column quarantines
      quants <- apply(t_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
      lines(  xs, quants[1,], col = "red")
      lines(  xs, quants[2,], col = "red")
      lines(  xs, quants[3,], col = "red")
      length_x = .3 *( tail(input_p[,1],1) - input_p[1,1] )

      # add posterior samples from the target
      for (i in 1:length(tie_points$tar_tie)){
        d_kde <- density(tar[i,], bw = "nrd0")
        lines(tie_points$inp_tie[i] + length_x * d_kde$y,d_kde$x,col=rgb(0,0,1,.8))
        lines(tie_points$inp_tie[i] - length_x * d_kde$y,d_kde$x,col=rgb(0,0,1,.8))
      }
      # add the proxy plot
      target_p$norm_val <- (target_p$value - min(target_p$value))/sd(target_p$value)
      input_p$ages_inp <- approx(xs,quants[2,],input_p$depth)$y
      input_p$norm_val <- (input_p$value - min(input_p$value))/sd(input_p$value)
      d_length <- .01 * (tail(input_p$depth,1) - input_p$depth[1])
      # plot the proxies
      lines(d_length * input_p$norm_val,input_p$ages_inp ,col=rgb(0,0,1,.4))
      lines(d_length * target_p$norm_val,target_p$ages_tar,col=rgb(0,1,.01,.4))
      legend("topleft", legend = c("Input", "Target"), 
             col = c(rgb(0,0,1,.4), rgb(0,1,.01,.4)), 
             lty = 1, lwd = 2, bty = "n")
      
  }
  
  # Define the sequence of x values for interpolation
  x_values <- seq(min(xs), max(xs), by = d_by)
  
  # Interpolate the minimum, medium, and maximum ages
  age_min_interp <- approx(x = output$elbows, y = output$age_min, xout = x_values)
  age_medium_interp <- approx(x = output$elbows, y = output$age_medium, xout = x_values)
  age_max_interp <- approx(x = output$elbows, y = output$age_max, xout = x_values)
  age_mean_interp <- approx(x = output$elbows, y = output$age_means, xout = x_values)
  
  # Create a dataframe with the interpolated values and x values
  interpolated_df <- data.frame(
    Depth = x_values,
    Min = age_min_interp$y,
    Medium = age_medium_interp$y,
    Max = age_max_interp$y
  )
  
  output$ages <- interpolated_df
  
  
  
  
  
  message("\n We're all out of BISCUITs - your integrated chronology is now complete! ")
  return(output)
}



# 
n_tie_points=2

Target = 'GB0220'
Input = 'GB0619'


# Check run time.
start.time1 = Sys.time()

out = Biscuit(Input=Input,Target=Target,'~/Documents/Biscuit/',
                   n_sections = TRUE, n_tie_points, sampling_date = 2020.1,
                   thin=10,burn=3e+3,iters=3e+3,
                   shape_acc = 1.5,  mean_acc = 20,
                   strength_mem = 20, mean_mem = .5,
                   run_target = FALSE,
                   target_age_model = "Bacon",
                   target_density_method = "KDE"
                  )

end.time1 = Sys.time()
total_time1 = end.time1 - start.time1
print(total_time1)
