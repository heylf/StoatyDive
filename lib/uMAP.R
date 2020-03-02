
#require("umap")
require("uwot")
require("data.table") # for the shift function
require("zoo") # for rollmeans

# new libraries 
require("corrplot")
require("kernlab")

# for testing
# require("Rtsne")
# require("kohonen")

umap_main <- function(data_path, filename, lam, maximal_cluster_number, optimal_k, on_off_smoothing){
  
  print(data_path)
  print(filename)
  print(paste0("Lam: ", lam))
  print(paste0("Maximal k cluster: ", maximal_cluster_number))
  if ( optimal_k == -1  ){
    print("Find k optimal cluster")
  } else {
    print(paste0("Set k cluster to: ", optimal_k))
  }
  if ( on_off_smoothing ){
    print("Turn on smoothing") 
  }

  data_path <- "~/Documents/StoatyDiveResults/SLBP"
  filename <- "SLBP_rep2_sorted"
  lam <- 0.3
  maximal_cluster_number <- 15
  on_off_smoothing <- TRUE
  optimal_k <- 6
  
  ############
  ##  Input ##
  ############
  
  if (dir.exists(paste0(data_path, "/clustering_", filename)) == FALSE ){
    dir.create(paste0(data_path, "/clustering_", filename))
  }
  output_path <- paste0(data_path, "/clustering_", filename, "/")
  
  # get data
  data <- read.table(paste0(data_path, "/data_classification_", filename, ".tsv"), header=FALSE)
  data2 <- matrix(as.numeric(unlist(data)), nrow=nrow(data), ncol=ncol(data))
  rownames(data2) <- seq(1, nrow(data2))
  
  outputmatrix <- read.table(paste0(data_path, "/CV_tab_", filename , ".bed"), header=FALSE)
  
  ################
  ##  Normalize ##
  ################
  set.seed(123)
  
  min_max_noramlize <- function(x){
    x_new <- (x-min(x)) / (max(x)-min(x)) 
    return(x_new)
  }
  
  data_normalized <- t(apply(data2, 1, min_max_noramlize))
  
  # Set NaNs to zero
  data_normalized[which(is.na(data_normalized))] <- 0
  
  # Remove duplicates for uMAP
  duplicates <- which(duplicated(data_normalized))
  data_removed_duplicates <- data_normalized[which(duplicated(data_normalized) == FALSE),]
  
  #################
  ##  Smooth Out ##
  #################
  
  smoothing <- function(y, lambda, dim){
    data_points <- length(y)
    
    # Test if we have just a vector of constant values.
    # Lambda is part of the splines ridge regression.
    which(y != max(y))
    if ( length(which(y != max(y))) != 0 ){
      x <- c(1:data_points)
      y.smooth <- smooth.spline(x, y, spar=lambda)
      x_new <- seq(1, data_points, data_points/dim)
      y_new <- predict(y.smooth, x_new)$y
      y_new <- y_new/max(y_new)
      return(y_new)
    } 
    
    x_new <- seq(1, data_points, data_points/dim)
    y_new <- rep( max(y), length(x_new) )
    return(y_new)
  }
  
  smoothing_dim <- 150
  
  if ( ncol(data_removed_duplicates) > smoothing_dim ) {
    smoothing_dim <- ncol(data_removed_duplicates)
  }
    
  data_smoothed <- data_removed_duplicates
  if( on_off_smoothing ) {
    print("[NOTE] Smooth data")
    data_smoothed <- t(apply(data_removed_duplicates, 1, smoothing, lambda=lam, dim=smoothing_dim))
    
    # get rid of negative values
    data_smoothed[which(data_smoothed < 0)] = 0.0
    
    # get rid of diract impulses bigger than max (> 1.0)
    data_smoothed[which(data_smoothed > 1.0)] = 1.0
    
  }
  
  data_ready <- data_smoothed
  
  #######################
  ##  kmeans Functions ##
  #######################
  
  kmeansAIC <- function(fit){
    m <- ncol(fit$centers)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(D + 2*m*k)
  }
  
  # function to find optimal number of clusters
  kopt_means <- function(data, k.max){
    
    num_centroids <- c(1:k.max) 
    
    kmeans <- sapply(num_centroids, function(k){kmeans(data, centers = k, nstart = 100, iter.max = 10000)})
    
    # get total within-cluster sum of squares
    twss <- numeric(ncol(kmeans))
    for ( i in 1:ncol(kmeans) ){
      twss[i] <- kmeans[,i]$tot.withinss
    }
    
    # get percentage of variance explained
    perc <- numeric(ncol(kmeans))
    for ( i in 1:ncol(kmeans) ){
      perc[i] <- (kmeans[,i]$betweenss / kmeans[,i]$totss) * 100
    }
    
    aic <- apply(kmeans, 2, kmeansAIC)
    
    opt <- which.min(aic) - 1
    
    differences <- numeric()
    for( i in 1:(length(twss)-1) ) {
      differences <- c(differences, abs(twss[i] - twss[i+1]))
    }
    percantage_diff <- differences/sum(differences)
    
    if (  is.na(which( percantage_diff < 0.1 )[1] ) == FALSE ){
      if ( which( percantage_diff < 0.1 )[1] < opt ) {
        opt <- which( percantage_diff < 0.1 )[1] + 1
      }
    }
    
    pdf(paste0(output_path, "/kmeans_Optimization.pdf"), 15, 15)
    par(cex = 2.0, mfrow = c(2,2), cex = 1.5, family = 'serif')
    
    plot(num_centroids, twss,
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab="Total within Cluster Sum of Squares",
         xlim=c(1,k.max))
    abline(v = opt, lty =2)
    
    plot(num_centroids, perc,
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab=" ",
         ylim = c(0,100),
         xlim=c(1,k.max))
    text( (opt - 0.5) , perc[opt], round(perc[opt], digits = 1))
    abline(v = opt, lty =2)
    
    plot(num_centroids, log(aic),
         type="b", pch = 19, frame = FALSE,
         xlab="Number of Clusters",
         ylab="AIC (log)",
         xlim=c(1,k.max))
    abline(v = opt, lty =2)
    dev.off()
    
    return(c(opt, perc[opt]))
  }
  
  ###################################
  ### Adding Additional Features  ###
  ###################################
  
  AUC <- function(y){
    x <- c(1:length(y))
    return( sum( diff(x)*abs(rollmean(y,2)) ) ) 
  }
  
  Inflection <- function(y){
    infl <- c(FALSE, diff(diff(y)>0)!=0)
    return(length(which(infl == TRUE)))
  }
  
  MSE <- function(sim, obs){
    return( ( sim - obs)^2 )
  }
  
  ARC <- function(y){
    x <- c(1:length(y))
    return(sum(sqrt(diff(x)^2 + diff(y)^2)))
  }
  
  # Check for Pleateau
  number_of_constant_value <- apply(data_removed_duplicates, 1, function(x){ return(length(which(x == 1)))} )
  data_ready <- cbind(data_ready, number_of_constant_value)
  
  # Are under the curve
  auc <- apply(data_smoothed, 1, AUC)
  data_ready <- cbind(data_ready, auc)

  # Arc Length Feature
  arc_length <- apply(data_removed_duplicates, 1, ARC)
  data_ready <- cbind(data_ready, arc_length)
  
  # OLD Features (many useful for the future)
  # Number of inflexion points
  # data_smoothed_for_inflection <- t(apply(data_removed_duplicates, 1, smoothing, lambda=.3, dim=smoothing_dim))
  # inflections_points <- apply(data_smoothed, 1, Inflection)
  # data_ready <- cbind(data_ready, inflections_points)
  
  # # MSE Features
  # data_smoothed_for_mse <- t(apply(data_removed_duplicates, 1, smoothing, lambda=.8, dim=ncol(data_removed_duplicates)))
  # errors <- sapply(1:nrow(data_ready), function(r){ sum((data_removed_duplicates[r,]-data_smoothed_for_mse[r,])^2) })
  # mse <- errors/smoothing_dim
  # data_ready <- cbind(data_ready, mse)
  
  # Max Absolute Difference Feature
  # max_abs_diff <- apply( data_removed_duplicates, 1, function(x){max(abs(diff(x)))} )
  # data_ready <- cbind(data_ready, max_abs_diff)
  
  # Variance Features
  variance <- apply(data_removed_duplicates, 1, var)
  data_ready <- cbind(data_ready, variance)
  
  # #########
  # ## PCA ##
  # #########
  # 
  # # Just PCA for comparison
  # pdf(paste0(output_path,"/PCA.pdf"))
  # par(family = 'serif', cex = 1.5, mfrow=c(2, 2), mgp = c(2, 1, 0))
  # pca <- prcomp(data_ready)$x
  # plot(pca[,1], pca[,2], xlab=paste0("PC ",1), ylab=paste0("PC ",2))
  # plot(pca[,1], pca[,3], xlab=paste0("PC ",1), ylab=paste0("PC ",3))
  # plot(pca[,2], pca[,3], xlab=paste0("PC ",2), ylab=paste0("PC ",3))
  # dev.off()
  # 
  # ##########
  # ## tSNE ##
  # ##########
  # 
  # tsne <- Rtsne(data_ready, dims = 2, perplexity=100, verbose=TRUE, max_iter = 5000, eta=100, momentum=0.1, pca = TRUE)
  # pdf(paste0(output_path,"/tSNE.pdf"))
  # par(family = 'serif', cex = 1.5, mgp = c(2, 1, 0))
  # plot(tsne$Y[,1], tsne$Y[,2], xlab=paste0("Dim ",1), ylab=paste0("Dim ",2))
  # dev.off()
  # 
  # #########
  # ## SOM ##
  # #########
  # 
  # # Create the SOM Grid - you generally have to specify the size of the
  # # training grid prior to training the SOM. Hexagonal and Circular
  # # topologies are possible
  # som_grid <- somgrid(xdim = 10, ydim= 10, topo="hexagonal")
  # 
  # # Finally, train the SOM, options for the number of iterations,
  # # the learning rates, and the neighbourhood are available
  # set.seed(123)
  # 
  # som_model <- som(data_ready,
  #                  grid=som_grid,
  #                  rlen=5000,
  #                  alpha=c(0.05, 0.01),
  #                  radius=c(5, 0),
  #                  dist.fcts="euclidean")
  # 
  # pdf(paste0(output_path, "/SOM.pdf"))
  # par(family = 'serif')
  # plot(som_model, type="count")
  # dev.off()
  
  #############
  ### uMAP  ###
  #############
  
  d <- 2
  # custom.config = umap.defaults
  # custom.config$n_epochs = 5000
  # custom.config$n_components = d
  # custom.config$min_dist = 0.01
  # custom.config$n_neighbors = 5
  
  #data_umap = umap(data_ready, config = custom.config)
  #new_data <- data_umap$layout
  data_umap = uwot::umap(data_ready, n_epochs=5000, n_components=d, min_dist=0.01, n_neighbors=5)
  new_data <- data_umap
  
  # Remove constant value profiles
  constant_profiles <- apply( data_removed_duplicates, 1, function(x){return(length(which(x!=0)))} )
  non_constant_profiles <- which(constant_profiles!=0)
  constant_profiles <- which(constant_profiles==0)
  new_data_for_kmeans <- new_data[non_constant_profiles,]
  
  ###############
  ### kmeans  ###
  ###############

  returned_optimisation_values <- kopt_means(new_data_for_kmeans, maximal_cluster_number)
  optimal_num_centroids <- returned_optimisation_values[1]
  optimal_perc <- returned_optimisation_values[2]

  if ( optimal_k > 0 ){
    optimal_num_centroids <- optimal_k
  }

  ## use kmeans clustering with optimal number of centroids
  kmeans_obj <- kmeans(new_data_for_kmeans, centers=optimal_num_centroids, nstart = 10000, iter.max = 10000)
  clusters <- kmeans_obj$cluster

  num_clusters <- length(unique(clusters))

  kmeans_centers <- kmeans_obj$centers
  memberships <- matrix(nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids)
  for(i in 1:nrow(kmeans_centers)){
    memberships[,i] <- apply(new_data_for_kmeans, 1, function(x){dist(rbind(kmeans_centers[i,],x))})
  }
  prob_memberships <- memberships
  prob_memberships <- apply(memberships, 1, function(x){1-(x/max(x))})
  prob_memberships <- matrix(prob_memberships, nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids, byrow=TRUE)
  prob_memberships <- apply(prob_memberships, 1, function(x){x/sum(x)})
  prob_memberships <- matrix(prob_memberships, nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids, byrow=TRUE)
  
  fc_memberships <- prob_memberships
  
  # Add cluster of constant values
  if ( length(constant_profiles) != 0 ){
    num_clusters <- num_clusters + 1
    fc_add_constant_cluster <- rep(-1, nrow(new_data))
    fc_add_constant_cluster[non_constant_profiles] <- clusters
    fc_add_constant_cluster[constant_profiles] <- num_clusters
    fc_clusters <- fc_add_constant_cluster

    fc_memberships <- matrix(ncol=num_clusters, nrow=nrow(new_data))
    fc_memberships[non_constant_profiles,] <- cbind(prob_memberships, rep(0.0, nrow(prob_memberships)))
    fc_memberships[constant_profiles,] <- c(rep(0.0, optimal_num_centroids), 1.0)
    colnames(fc_memberships) <- c(1:num_clusters)
  }
  clusters <- fc_clusters

  pdf(paste0(output_path, "/fuzzy_clustering_overview.pdf"), width = 2, height = 30)
  par(family="serif")
  if( nrow(prob_memberships) > 200 ){
    corrplot(prob_memberships[1:200,], is.corr=FALSE, cl.lim=c(0.0, 1.0), cl.ratio=0.5)
  } else {
    corrplot(prob_memberships, is.corr=FALSE, cl.lim=c(0.0, 1.0))
  }
  dev.off()
 
  ############################
  ### Spectral Clustering  ###
  ############################
  # 
  # returned_optimisation_values <- kopt_means(new_data_for_kmeans, maximal_cluster_number)
  # optimal_num_centroids <- returned_optimisation_values[1]
  # optimal_perc <- returned_optimisation_values[2]
  # 
  # if ( optimal_k > 0 ){
  #   optimal_num_centroids <- optimal_k
  # }
  # 
  # optimal_num_centroids <- 5
  # 
  # ## use kmeans clustering with optimal number of centroids
  # spectral_obj <- NULL
  # attempt <- 1
  # while( is.null(spectral_obj) && attempt <= 100 ) {
  #   attempt <- attempt + 1
  #   try(
  #     spectral_obj <- kkmeans(new_data_for_kmeans, centers=optimal_num_centroids, kernel="besseldot", iterations=5000)
  #   )
  # } 
  # 
  # clusters <- spectral_obj[1:length(spectral_obj)]
  # 
  # num_clusters <- length(unique(clusters))
  # 
  # spectral_centers <- centers(spectral_obj)
  # memberships <- matrix(nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids)
  # for(i in 1:nrow(spectral_centers)){
  #   memberships[,i] <- apply(new_data_for_kmeans, 1, function(x){dist(rbind(spectral_centers[i,],x))})
  # }
  # prob_memberships <- memberships
  # prob_memberships <- apply(memberships, 1, function(x){1-(x/max(x))})
  # prob_memberships <- matrix(prob_memberships, nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids, byrow=TRUE)
  # prob_memberships <- apply(prob_memberships, 1, function(x){x/sum(x)})
  # prob_memberships <- matrix(prob_memberships, nrow=nrow(new_data_for_kmeans), ncol=optimal_num_centroids, byrow=TRUE)
  # 
  # fc_memberships <- prob_memberships
  # 
  # # Add cluster of constant values
  # if ( length(constant_profiles) != 0 ){
  #   num_clusters <- num_clusters + 1
  #   fc_add_constant_cluster <- rep(-1, nrow(new_data))
  #   fc_add_constant_cluster[non_constant_profiles] <- clusters
  #   fc_add_constant_cluster[constant_profiles] <- num_clusters
  #   fc_clusters <- fc_add_constant_cluster
  #   
  #   fc_memberships <- matrix(ncol=num_clusters, nrow=nrow(new_data))
  #   fc_memberships[non_constant_profiles,] <- cbind(prob_memberships, rep(0.0, nrow(prob_memberships)))
  #   fc_memberships[constant_profiles,] <- c(rep(0.0, optimal_num_centroids), 1.0)
  #   colnames(fc_memberships) <- c(1:num_clusters)
  # }
  # clusters <- fc_clusters
  # 
  # pdf(paste0(output_path, "/fuzzy_clustering_overview.pdf"), width = 5, height = 15)
  # if( nrow(prob_memberships) > 100 ){
  #   corrplot(prob_memberships[1:100,], is.corr=FALSE, cl.lim=c(0.0, 1.0))
  # } else {
  #   corrplot(prob_memberships, is.corr=FALSE, cl.lim=c(0.0, 1.0))
  # }
  # dev.off()
  # 
  #########################
  ##  Visualize Clusters ##
  #########################
  
  colors <- rainbow(num_clusters)
  cluster_colors <- unlist(lapply(clusters, function(x) { return(colors[x]) }))
  seq_clusters <- c(1:num_clusters)
  
  plot_uMAP <- function(c1, c2, data){
    par(mar=c(4.1, 4.1, 4.1, 4.1), xpd=TRUE, mgp = c(2, 1, 0))
    plot(data[,c1], data[,c2], xlab=paste0("Dim ",c1), ylab=paste0("Dim ",c2), col="white")
    legend("topright", inset=c(-0.2,0), col=colors[seq_clusters], legend=seq_clusters, pch=seq_clusters, title="Group")
    
    for(k in 1:num_clusters){
      datapoints <- which(clusters == k)
      points(data[datapoints,c1], data[datapoints,c2], col=colors[k], pch=k)
    }
  }
  
  pdf(paste0(output_path,"/uMAP.pdf"))
  par(family = 'serif', cex = 1.5)
  
  plot_uMAP(1, 2, new_data)
  
  if (d >= 3){
    plot_uMAP(1, 3, new_data)
    plot_uMAP(2, 3, new_data)
  }
  
  dev.off()
  
  ##########################
  ## Create Output Matrix ##
  ##########################
  
  # match duplicates with non duplicates 
  uMAP_peaks_corresponding_to_duplicate <- match(data.frame(t(data_normalized[duplicates,])), data.frame(t(data_removed_duplicates)))
  clusters_of_duplicates <- clusters[uMAP_peaks_corresponding_to_duplicate]
  
  # Generate cluster column
  cluster_col <- matrix(ncol = 1, nrow=nrow(data))
  cluster_col[which(duplicated(data_normalized) == FALSE), 1] <- clusters
  cluster_col[duplicates, 1] <- clusters_of_duplicates
  
  # Fill output matrix 
  summit_positions <- apply(data2, 1, which.max)
  outputmatrix <- cbind(outputmatrix, summit_positions)
  outputmatrix <- cbind(outputmatrix, rep(-1, nrow(outputmatrix)))
  outputmatrix[match(duplicates, outputmatrix[,13]), ncol(outputmatrix)] <- clusters_of_duplicates
  outputmatrix[match(rownames(data_removed_duplicates), outputmatrix[,13]), ncol(outputmatrix)] <- clusters
  
  write.table(outputmatrix, file=paste0(data_path, "/final_tab_", filename, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
  
  #########################
  ## Create Fuzzy Matrix ##
  #########################
  
  fuzzy_data_table <- matrix(ncol=ncol(fc_memberships), nrow=nrow(data))
  
  fuzzy_data_table[which(duplicated(data_normalized) == FALSE),] <- fc_memberships
  fuzzy_data_table[duplicates,] <- fc_memberships[uMAP_peaks_corresponding_to_duplicate,]
    
  colnames(fuzzy_data_table) <- paste0("cluster_", c(1:ncol(fc_memberships)))
  rownames(fuzzy_data_table) <- paste0("peak_", c(1:nrow(data)))
  
  write.table(round(fuzzy_data_table,4), file=paste0(data_path, "/fuzzy_clustering_table.tsv"), row.names = TRUE, col.names = NA, quote = FALSE, sep="\t")
  
  ########################################
  ## Generate Peak Profiles of Clusters ##
  ########################################
  
  print("[NOTE] Create cluster profiles")
  
  # Smoothed profiles
  for_smoothed_profiles <- t(apply(data_normalized, 1, smoothing, lambda=lam, dim=smoothing_dim))
  for_smoothed_profiles[which(for_smoothed_profiles < 0)] = 0.0
  for_smoothed_profiles[which(for_smoothed_profiles > 1.0)] = 1.0

  peak_length <- length(data_normalized[1,])
  unique_clusters <- unique(clusters)
  num_clusters <- length(unique_clusters)
  for( i in 1:num_clusters ){
    peaks <- which(cluster_col[,1] == unique_clusters[i])
    num_peaks <- length(peaks)
    
    if ( num_peaks >= 20 ){
      num_peaks <- 20
    }
    
    if( on_off_smoothing ) {
      pdf(paste0(output_path, "/cluster_smoothed", unique_clusters[i],".pdf"))
      par(family = 'serif', cex = 1.5, mgp = c(2, 1, 0))
      for( j in 1:num_peaks ){
        barplot(for_smoothed_profiles[peaks[j],], ylab = "Normalized Read Count", xlab = "Relative Nucleotide Position", 
               col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA)
        abline(v=round(smoothing_dim/2), lty="dashed")
        mtext("0", side=1, line=0.1, at=round(smoothing_dim/2)) 
      }
      dev.off()
    }
    
    pdf(paste0(output_path, "/cluster_", unique_clusters[i],".pdf"))
    par(family = 'serif', cex = 1.5, mgp = c(2, 1, 0))
    for( j in 1:num_peaks ){
      barplot(data_normalized[peaks[j],], ylab = "Normalized Read Count", xlab = "Relative Nucleotide Position", 
           col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA)
      abline(v=round(peak_length/2), lty="dashed")
      mtext("0", side=1, line=0.1, at=round(peak_length/2)) 
    }
    dev.off()
    
  }
  
  # Sort clusters and colors
  unique_clusters <- sort(unique_clusters)
  
  # Overviews
  #peak_selection <- c(1,1,1,1,1,1,1)
  peak_selection <- rep(1, num_clusters)
  peaks <- which(cluster_col[,1] == unique_clusters[1])
  
  xlab_position <- ceiling(num_clusters/2) 
  
  pdf(paste0(output_path, "/overview_clusters.pdf"), width = 15, height = 3)
  par(family = 'serif', cex = 1.5, mfrow=c(1, num_clusters), mar=c(4,3,2,0), mgp = c(2, 1, 0))
  barplot(data_normalized[peaks[peak_selection[1]],], ylab = "Normalized Read Count", xlab = "", 
          col = colors[1], pch=20, ylim=c(0,1), space=0, border=NA)
  abline(v=round(peak_length/2), lty="dashed")
  mtext("0", side=1, line=0.5, at=round(peak_length/2), family = 'serif', cex = 0.7)
  title(main = paste0("1: ", length(peaks)), line=0.9)
    
  for( i in 2:num_clusters ){
    peaks <- which(cluster_col[,1] == unique_clusters[i])
    par(mar=c(4,1,2,1), mgp = c(2, 1, 0))
    if ( i == xlab_position ){
      barplot(data_normalized[peaks[peak_selection[i]],], ylab = "", xlab = "Relative Nucleotide Position", 
              col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA, yaxt="n")
    } else {
      barplot(data_normalized[peaks[peak_selection[i]],], ylab = "", xlab = "", 
              col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA, yaxt="n")
    }
    abline(v=round(peak_length/2), lty="dashed")
    mtext("0", side=1, line=0.5, at=round(peak_length/2), family = 'serif', cex = 0.7)
    title(main = paste0(i, ": ", length(peaks)), line=0.9)
  }
  
  dev.off()
  
  # Average Profiles
  for_avaerage_profile <- data_normalized
  
  peaks <- which(cluster_col[,1] == unique_clusters[1])
  
  pdf(paste0(output_path, "/cluster_average_profiles.pdf"), width = 15, height = 3)
  par(family = 'serif', cex = 1.5, mfrow=c(1, num_clusters), mar=c(4,3,2,0), mgp = c(2, 1, 0))
  
  average_profile <- apply(for_avaerage_profile[peaks,], 2, mean)
  barplot(average_profile, ylab = "Normalized Read Count", xlab = "",
          col = colors[unique_clusters[1]], pch=20, ylim=c(0,1), space=0, border=NA)
  abline(v=round(peak_length/2), lty="dashed")
  mtext("0", side=1, line=0.1, at=round(peak_length/2),  family = 'serif', cex = 0.7)
  title(main = paste0("1: ", length(peaks)), line=0.9)
  
  for( i in 2:num_clusters ){
    peaks <- which(cluster_col[,1] == unique_clusters[i])
    average_profile <- apply(for_avaerage_profile[peaks,], 2, mean)
    
    par(mar=c(4,1,2,1), mgp = c(2, 1, 0))
    if ( i == xlab_position ){
      barplot(average_profile, ylab = "", xlab = "Relative Nucleotide Position", 
              col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA, yaxt="n")
    } else {
      barplot(average_profile, ylab = "", xlab = "", 
              col = colors[unique_clusters[i]], pch=20, ylim=c(0,1), space=0, border=NA, yaxt="n")
    }
    abline(v=round(peak_length/2), lty="dashed")
    mtext("0", side=1, line=0.5, at=round(peak_length/2), family = 'serif', cex = 0.7)
    title(main = paste0(i, ": ", length(peaks)), line=0.9)
  }
  
  dev.off()

}
  
#################
##  Parameters ##
#################

args = commandArgs(trailingOnly=TRUE)

if ( length(args) < 6 ) {
  print("[TEST] Not enough arguments.")
} else {
  umap_main(args[1], args[2], as.numeric(args[3]), as.numeric(args[4]), as.numeric(args[5]), args[6])
}
