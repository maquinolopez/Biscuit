###########################################
###    Setting chronologies using 
###   tie points using iterations        
###    By Marco A. Aquino-Lopez         
###  BISCUIT: 
### Bayesian Integration of Sequences 
### Using Chronology Intersection Tie-points
###########################################

density_plot <- function(depths, tar_mt,xlabel,ylabel = "proxy units",add = FALSE,axis=TRUE,flip=FALSE){
  depths_bw <- depths[2] - depths[1] 
  t_max <- quantile(tar_mt, probs = c(0.999))
  t_min <- min(tar_mt)
  if (!add){
    if(axis){
      if (flip){
        plot(colMeans(tar_mt),depths , type='l', xlab=ylabel, ylab=xlabel, xlim = c(t_min,t_max) , ylim=c( depths[1],tail(depths,1)),col=rgb(0,0,0,1))  
      }else{
        plot(depths,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, ylim = c(t_min,t_max) ,  xlim=c( depths[1],tail(depths,1)),col=rgb(0,0,0,1))   
      }
      
    }else{
      if (flip){
        plot(colMeans(tar_mt),depths , type='l', ylab=xlabel, xlab=ylabel, xlim = c(t_min,t_max), ylim=c( depths[1],tail(depths,1)),col=rgb(0,0,0,1),xaxt = 'n')   
      }else{
        plot(depths,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, ylim = c(t_min,t_max), xlim=c( depths[1],tail(depths,1)),col=rgb(0,0,0,1),xaxt = 'n')    
      }
    }
  }else{
    if (flip){
      lines(colMeans(tar_mt),depths ,col=rgb(0,0,0,1)) 
    }else{
      lines(depths,colMeans(tar_mt) ,col=rgb(0,0,0,1))  
    }
  }
  
  for(i in 1:ncol(tar_mt) ){
    h <- hist(tar_mt[,i], plot = FALSE,breaks=150)
    cols <- gray(1-h$counts/max(h$counts),alpha = .4)
    # Plot non-zero rects 
    if (flip){
      rect(ybottom = depths[i], ytop = depths[i]+depths_bw,
           xleft = h$breaks[-1], #head(h$breaks, -1),
           xright = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
    }else{
      rect(xleft = depths[i], xright = depths[i]+depths_bw,
           ybottom = h$breaks[-1], #head(h$breaks, -1),
           ytop = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
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
                   # suggest_tie = TRUE, 
                   # tolerance_preproc = 10,
                   run_target = FALSE
                   ){
  source(paste0(folder,'/twalk.R'))
  
  # Prepare data
  colnames <- c("depth", "value")
  # Read CSV into input as dataframe
  input_p  <- read.csv(paste(folder,Input,'-',Target,'/',Input,"_proxy.csv",sep='') , header = TRUE, col.names = colnames)
  target_p <- read.csv(paste(folder,Input,'-',Target,'/',Target,"_proxy.csv",sep='') , header = TRUE, col.names = colnames)
  
  if (n_sections){
    n_sections = as.integer(n_tie_points * 4 )
  }

  # Load run and return the ages at the depths of the tie points
  # it return a matrix. 
  get_ages <- function(tiepoints,target) {
    library(rplum)
    Plum(target,run=run_target,coredir=paste0(folder,Input,'-',Target,'/'))
    agedepth() 
    
    depths <- tiepoints$tar_tie
    df <- matrix(NA, nrow = length(depths), ncol = length(Bacon.Age.d(depths[1])))
    cat('\n Getting the age posterior distributions')
    pb <- txtProgressBar(min = 0, max = length(depths), style = 3, char = "=") 
    for (i in 1:length(depths)) {
      df[i, ] <- as.numeric(Bacon.Age.d(depths[i]))
      setTxtProgressBar(pb, i)
    }
    return(df)
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
    
    # rainbow_colors <- rep(rgb(.8,.5,1,.4) , length(intersections_input)) #rainbow(length(intersections_input),alpha=.4)
    # 
    # for (i in 1:length(intersections_input)) {
    #   points(input[input[,1]==intersections_input[i],1],
    #          input[input[,1]==intersections_input[i],2], col=rainbow_colors[i], pch=16)
    # }
    # for (i in 1:length(intersections_target)) {
    #   points(target[target[,1]==intersections_target[i],1],
    #          target[target[,1]==intersections_target[i],2], col=rainbow_colors[i], pch=16)
    # }
    # break
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
  
  

  # Pre-analysis to suggest tie points
  # input_diff <- diff(min_max_normalize(input_p[,2]))
  # input_depth_diff  <- input_p[which(input_diff < -1),1]
  # 
  # target_diff <-  diff(min_max_normalize(target_p[,2]))
  # target_depth_diff  <- input_p[which(target_diff < -1),1]
  # 
  # distance_matrix <- matrix(nrow = length(input_depth_diff), ncol = length(target_depth_diff))
  # for (i in 1:length(input_depth_diff)) {
  #   for (j in 1:length(target_depth_diff)) {
  #     distance_matrix[i, j] <- abs(input_depth_diff[i] - target_depth_diff[j])
  #   }
  # }
  # 
  # Find the minimum distance in the matrix
  # min_distance <- min(distance_matrix)
  
  # # Initialize a list to store the intersections
  # intersections_input <- c()
  # intersections_target <- c()
  # intersections_distan <- c()
  # 
  # 
  # # Loop over the matrix to find values within the minimum distance plus tolerance
  # for (i in 1:nrow(distance_matrix)) {
  #   for (j in 1:ncol(distance_matrix)) {
  #     if (distance_matrix[i, j] <= min_distance + tolerance_preproc) {
  #       # Append the intersection information as a list to the intersections list
  #       intersections_input <- c(intersections_input , input_depth_diff[i])
  #       intersections_target <- c(intersections_target , target_depth_diff[j])
  #       intersections_distan <- c(intersections_distan , distance_matrix[i, j])
  #     }
  #   }
  # }
  
  # repeat to get the >1
  # input_depth_diff  <- input_p[which(input_diff > 1),1]
  # target_depth_diff  <- input_p[which(target_diff > 1),1]
  # distance_matrix <- matrix(nrow = length(input_depth_diff), ncol = length(target_depth_diff))
  # for (i in 1:length(input_depth_diff)) {
  #   for (j in 1:length(target_depth_diff)) {
  #     distance_matrix[i, j] <- abs(input_depth_diff[i] - target_depth_diff[j])
  #   }
  # }
  # 
  # Find the minimum distance in the matrix
  # min_distance <- min(distance_matrix)
  
  # Loop over the matrix to find values within the minimum distance plus tolerance
  # for (i in 1:nrow(distance_matrix)) {
  #   for (j in 1:ncol(distance_matrix)) {
  #     if (distance_matrix[i, j] <= min_distance + tolerance_preproc) {
  #       # Append the intersection information as a list to the intersections list
  #       intersections_input <- c(intersections_input , input_depth_diff[i])
  #       intersections_target <- c(intersections_target , target_depth_diff[j])
  #       intersections_distan <- c(intersections_distan , distance_matrix[i, j])
  #     }
  #   }
  # }
  # 
  # print(intersections_target)
  # print(intersections_input)
  # Get tie points
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

  kde_list <- target_density(tie_points$tar_tie,tar)
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
  
  ########## Functions ##########
  #####################

  # get the necessary functions
  # support function
  sampling_date = 1950 - sampling_date
  low <- c(sampling_date - 3*sampling_date_sd,
           0, # mem
           rep(0, n_sections-1)) # m_s
  up <- c( sampling_date + 3*sampling_date_sd,  #tao0
           1, # mem,
           rep(Inf, n_sections-1)) # alphas
  mx_age <- 100000
  
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
  
  # Age model functions
  ttime <- function(x, param,b_length){
    t0 <- param[1]
    # w <- param[2]
    slopes <- param[-c(1,2)]
    t_is <- c(t0, t0 + cumsum(slopes * b_length))
    return( approx(x_sections,t_is,x)$y )
  }
  
  
  # alphas
  alphas <- function(param, n_slopes){
    w <- param[2]
    slopes <- param[-c(1,2)]
    alf <- ( slopes[-n_slopes] - w * slopes[-1] ) / (1 - w)
    return(c(alf, tail(slopes,1)) )
  }
  
  
  # Prior density
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
  
  
  # likelihood  
  loglikelihood_uq <- function(params){
    t_times <- ttime(tie_points$inp_tie,params,b_length)
    ll <- 0 
    for (k in 1:length(t_times)){
      kde <- as.double(kde_list[[k]]['bw'])
      # lk <- l_target_kernel(d = k, new_x = t_times[k],kde = kde)
      lk <- l_target_kernelC(d = k, new_x = t_times[k],kde = kde)
      ll <- ll + lk # d, y,kde_list
    }
    return(ll)
  }

  n_kde <- length(tar[1,])
  sd_convertor <- 1/(1.06*n_kde^(-1/5))
  # Define a function to get the log density for a given depth and value y
  l_target_kernel <- function(d, new_x,kde) {
    # return(sum(dnorm(new_x , tar[d,],sd_convertor * kde$bw,log=T)) / n_kde )
    return(log( sum(dnorm(x = new_x , mean = tar[d,], sd = sd_convertor * kde )))  )
  }
  # This is the C implementation of the likelihood
  Rcpp::sourceCpp("~/GitHub/Biscut/targetDensity.cpp")
  # l_target_kernelC <- function(d, new_x,kde) {
  #   result <- lTargetKernelCpp(d, new_x, kde, tar[d,], sd_convertor)
  #   return( result )
  # }
  
  loglikelihood_uqC <- function(params){
    t_times <- ttime(tie_points$inp_tie,params,b_length)
    ll <-loglikelihood_uqCpp(params, t_times, kde_list,tar, sd_convertor )
    return(ll)
  }
  
  
  
  # objective function
  obj <- function(param){ 
    - (logprior(param) + loglikelihood_uqC(param)  ) 
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



out = Biscuit(Input=Input,Target=Target,'~/Documents/Biscuit/',
                   n_sections = TRUE, n_tie_points, sampling_date = 2020.1,
                   thin=50,burn=3e+3,iters=4e+3,
                   shape_acc = 1.5,  mean_acc = 20,
                   strength_mem = 20, mean_mem = .5,
                   run_target = FALSE)



