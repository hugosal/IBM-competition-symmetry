#parametros a probar

# 1  estructura espacial
# 2 variacion inicial en tamanios
# 3 densisdad de plantas
# 4 asimetria
# 5 cosas de especies?

# en este quiero probar las condiciones que dice Weiner1990 sobre el efecto
# de las interactiones y su simetria en la variacion de 

source("circles_area_of_overlap.R")
#source("get_ppis.R")

seed_value <-26
set.seed(seed_value)

#esta funcion va a hacer que un conjunto de coordenadas xy se ajusten 
# de tarl manera que center_x y center_y esten en el centro del mundo
# y asi el wrapping no afecte el calculo de areas, y regresa todas con la nueva
# escala, con el centro como primer elemento
center_world_arround <- function(center_xy, ws, coords){
  row_number <- which(rownames(coords)==as.character(center_xy))
  x_offset <- (ws/2) - coords$x[row_number]
  y_offset <- (ws/2) - coords$y[row_number]
  coords$x <- (coords$x + x_offset) %% ws
  coords$y <- (coords$y + y_offset) %% ws
  coords
}

plot_plantcomm <- function(com, numbers = FALSE, main="", ws = NULL){
  colores <- col2rgb(1:4)
  colores <- apply(colores, 2, FUN = function(x)rgb(x[1]/255,
                                                    x[2]/255,
                                                    x[3]/255,
                                                    alpha = 0.5)  )
  plot(com$x, com$y, col=com$ft, type="n", main=main, xlab="x", ylab="y",
       xlim = c(min(com$x) - max(com$sz), max(ws, max(com$x) + max(com$sz))),
       ylim = c(min(com$y) - max(com$sz), max(ws, max(com$y) + max(com$sz))))
  
  if (! is.null(ws)){
    abline(v = 0, h = 0)
    abline(v = ws, h = ws)
  }
  for (j in 1:nrow(com)){
  plotrix::draw.circle(x = com$x[j], y = com$y[j], radius = com$sz[j], 
                       col = colores[com$ft[j]])
    }
  if (numbers){
    text(x = com$x, y = com$y, labels = 1:nrow(com))
  }
}

generate_initial_points <- function(N, ws){
  if (sqrt(N) == round(sqrt(N)) ){
    n_per_row <- sqrt(N)
    grid_poinst <- seq(from = 0, to = ws, 
                       by = ws/(n_per_row))[1:n_per_row]
    x <- rep(grid_poinst, n_per_row)
    y <- as.vector(sapply(grid_poinst, function(x)rep(x, n_per_row)))
    return(cbind(x, y))
  }else if(sqrt(N/2) == round(sqrt(N/2))){
    i <- sqrt(N/2)
    separation <- ws/i
    init_grid  <- seq(from = 0,
                      to = ws - (separation/2), 
                      by = separation/2)
    n <- length(init_grid)
    is_even <- (1:n) %% 2 == 0
    coords <- matrix(ncol = 2, nrow = 0)
    for (x in 1:n){
      this_x <- init_grid[x]
      if(is_even[x]){
        yes <- init_grid[!is_even]
      }else{
        yes <- init_grid[is_even]
      }
      coords <- rbind(coords, matrix(c(rep(this_x, length(yes)),
                                       yes), nrow = length(yes)))}
    return(coords)
  }else{
    pts <- read.csv("configurations/point_configurations.csv", 
                    header = FALSE, row.names = 1)
    if (N %in% rownames(pts)){
      coords <- pts[N == rownames(pts), ]
      coords <- coords[coords != ""]
      coords <- do.call(rbind, lapply(coords, function(x) as.numeric(strsplit(x, ":")[[1]])))
      return(coords)
    }else{
      stop("Configuration not found yet")}
  }
}

ws <- 20 # world size 

timesteps <- 20   # length of each run

n_reps <- 15 # number of LHS samples

uniform_LHS_sample_from_range <- function(lower, upper, n_samples){
  limits <- seq(from = lower, to = upper, length.out = n_samples + 1)
  sapply(1:n_samples, function(x){ runif(n = 1, 
                                         min = limits[x], 
                                         max = limits[x + 1])})
  }

# cada fila va a ser una fila de parametros a usar
LHS_param <- matrix(c(sample(uniform_LHS_sample_from_range(lower = 0, 
                                                           upper = 1, 
                                                           n_samples =n_reps)), #wrl reach, 0 - 1
                      sample(uniform_LHS_sample_from_range(lower = 0.5, 
                                                           upper = 1.5, 
                                                           n_samples =n_reps)),#max initial size,
                      sample(uniform_LHS_sample_from_range(lower = 0.01, 
                                                           upper = 0.05, 
                                                           n_samples = n_reps)), #pop densiti 4 - 20/ws**2
                      sample(uniform_LHS_sample_from_range(lower = 0, 
                                                           upper = 2, 
                                                           n_samples =n_reps)),# competition asymmetry parameter
                      sample(uniform_LHS_sample_from_range(lower = 1, 
                                                           upper = 5, 
                                                           n_samples =n_reps)), # max growing size
                      sample(uniform_LHS_sample_from_range(lower = 0, 
                                                           upper = 2, 
                                                           n_samples =n_reps))), # max growth rate
                    nrow = n_reps)

colnames(LHS_param) <- c("world_reachability", 
                         "max_initial_sz", 
                         "pop_density",
                         "comp_symmetry", 
                         "max_growing_size",
                         "max_growth_rate")

output_meassure <- numeric(n_reps)

for (rep in 1:nrow(LHS_param)){
  print(paste( "rep ", rep))
  
  # extract value of model parameters
  world_reachablity <- LHS_param[rep, 1]
  max_initial_size <- LHS_param[rep, 2]
  pop_density <- LHS_param[rep, 3]
  theta <- LHS_param[rep, 4]
  max_S <- LHS_param[rep, 5]
  max_grwth_rt <- LHS_param[rep, 6]
  
  # set population size
  initial.n <- round(ws**2 * pop_density)
  config_found <- c(4,5,8,9,10,13,15,16,17,18,20,24,25,26,29,32,34,35,36,37,40)
  
  if (! initial.n %in% config_found){
    print("adj")
    initial.n <- config_found[which.min(abs(initial.n - config_found))]
  }
  print(initial.n)
  
  # set plants initial coordinates
  coordinates <- generate_initial_points(N = initial.n, ws = ws)
  
  plantcomm = data.frame(
    x = coordinates[, 1],
    y = coordinates[, 2],
    ft = 1,
    sz = runif(initial.n, min = 0.5, max = max_initial_size))
  
  # modify plant coordinates accoding to spatial disarrangement parameter
  
  random_angles <- runif(nrow(plantcomm), 
                         min = 0, max = 2*pi)
  
  plantcomm$x <- plantcomm$x + cos(random_angles)*world_reachablity
  plantcomm$y <- plantcomm$y + sin(random_angles)*world_reachablity
  

  # plot_plantcomm(plantcomm, numbers = TRUE, ws = ws,
  #                main=paste("density", round(pop_density,2),
  #                           "reachabil", round(world_reachablity,2)))

  a <- (-4*max_grwth_rt)/(max_S**2)
  b <- (4*max_grwth_rt)/max_S
  # a < 0 and b > 0
  
  
  for (t in 1:timesteps){
    #plot_plantcomm(plantcomm, numbers = TRUE, ws = ws, main=t)
    print(paste( "t ", t))
  
  	n = nrow(plantcomm)
  	
  	dx = as.matrix(dist(plantcomm$x,diag=TRUE,upper=TRUE))
  	dx[dx > ws/2] = ws - dx[dx > ws/2]   		### this wraps the interaction effects
  	dy = as.matrix(dist(plantcomm$y,diag=TRUE,upper=TRUE))
  	dy[dy > ws/2] = ws - dy[dy > ws/2]
  	dists = sqrt(dx^2 + dy^2)
  
  	radius_sum <- outer(X = plantcomm$sz, Y = plantcomm$sz, FUN = "+")
    
  	# which plant area of effect overlap
  	disteffect <- dists < radius_sum
  	n_in_groups <- apply(disteffect, MARGIN = 1, FUN = sum)
  	
  	# this optimization trick is to coimpute first intersections
  	# with most circles, then subsets of those groups will be extracted
  	# from memory instead of computed again
  	most_to_least_crowded <- order(n_in_groups, decreasing = T)

  	intersect_to_compute <- (1:n)[most_to_least_crowded]
  	intersect_to_compute <- unname(intersect_to_compute[n_in_groups[most_to_least_crowded] > 1])
  	# compute unbound growth rate, if this plant can grow further, compute
  	# intersections areas. Otherwise, it will stay the same size, computing areas
  	# will be irrelevant and thus it is omitted
  	
  	cant_grow_more <- which(plantcomm$sz >= max_S)
  		
  	# initiate vector of resources obtained with the area of effect of each plant
  	resources_obtained_p_ind <- plantcomm$sz**2 * pi
  
  	intersections_list_memory <-list()
  
  	for(j in intersect_to_compute){
  	  
  	  if (j %in% cant_grow_more){
  	    if((a * plantcomm$sz[j]**2) + (b * plantcomm$sz[j]) > 1e-6){
  	      stop("Algo mal con optimizacion de crecimiento")
  	    }
  	    next
  	  }
  	  
  	  names_original_plantcomm <- which(disteffect[j, ])
  
  	  previously_seen <- unlist(lapply(seq_along(intersections_list_memory),
  	                            function(x, names_this_set){
  	                              set_memory <- intersections_list_memory[[x]]$"names_plant"
  	                              if (length(names_this_set) == length(intersect(names_this_set, set_memory))){
  	                                x
  	                              }
  	                            },
  	                            names_this_set = names_original_plantcomm))
  
  	  if (!is.null(previously_seen)){ # if this set of intersections is recorded 
  	    intersections <- intersections_list_memory[[previously_seen[1] ]]$"intersections"
  	    names_original_plantcomm_memory <- intersections_list_memory[[previously_seen[1] ]]$"names_plant"
  
  	    resources_obtained <- sapply(seq_along(intersections),
  	                                 function(x){
  	                                   contestants <- as.numeric(strsplit(names(intersections[x]), split = ":")[[1]])
  	                                   if (! j %in% names_original_plantcomm_memory[contestants]){
  	                                     return(0)
  	                                   }
  	                                   controled_percent <-plantcomm$sz[j]**theta/sum(plantcomm$sz[names_original_plantcomm_memory[contestants]]**theta)
  	                                   intersections[x] * controled_percent
  	                                 })
  
  	    resources_obtained_p_ind[j] <-sum(resources_obtained)
  
  
  	  }else{
  
    	  plantcomm_sub <- plantcomm[names_original_plantcomm, ]
  
    	  plantcomm_sub[, c(1,2)] <- center_world_arround(coords = plantcomm_sub[, c(1,2)],
    	                                                  center_xy = j, ws = ws)
    	  # graficame_esta(centers_x = plantcomm_sub$x,
    	  #                centers_y = plantcomm_sub$y,
    	  #                radii = plantcomm_sub$sz)
    	  
    	  if (length(names_original_plantcomm) > 8){
    	    print(paste(length(names_original_plantcomm),", grt progress ", plantcomm$sz[j]/max_S))
    	    }
  
    	  intersections <- unlist(Librino_N(centers_x = plantcomm_sub$x,
    	                                    centers_y = plantcomm_sub$y,
    	                                    radii = plantcomm_sub$sz), use.names = T)
    	  
    	  if (is.numeric(validate_Librino(intersections, radii = plantcomm_sub$sz))){
    	    stop("Something wrong with intersections")
    	  }
    	  
    	  intersections <- intersections[intersections > 1e-5]
    
    	  resources_obtained <- sapply(seq_along(intersections),
      	                               function(x){
      	                                 contestants <- as.numeric(strsplit(names(intersections[x]), split = ":")[[1]])
      	                                 if (! j %in% names_original_plantcomm[contestants]){
      	                                   return(0)
      	                                 }
      	                                 controled_percent <- plantcomm$sz[j]**theta/sum(plantcomm$sz[names_original_plantcomm[contestants]]**theta)
      	                                 intersections[x] * controled_percent
      	                               })
  
      	  resources_obtained_p_ind[j] <- sum(resources_obtained)
  
      	  if (length(names_original_plantcomm) > 4){
  
      	    intersections_list_memory[[length(intersections_list_memory) + 1]]<- list(
      	      "names_plant" = names_original_plantcomm,
      	      "intersections" = intersections)
      	    }
      	  }
    	  }
  	# plants increase in size asymptotically
  
  
  	if(!all(resources_obtained_p_ind <= plantcomm$sz**2*pi)){
  	  stop("Something wrong with resource partitions")
  	}
  
  	# growth rate according to  10.1093/oxfordjournals.aob.a086287
  	
  	ind_unbound_growth_rate <-  (a * plantcomm$sz**2) + (b * plantcomm$sz)
  	
  	growthratemodifier <- resources_obtained_p_ind/(plantcomm$sz**2*pi)
  
  	indgr <- growthratemodifier * ind_unbound_growth_rate 
  
  	indgr[indgr < 0] <- 0
  
  	plantcomm$sz = plantcomm$sz + indgr
  
  
  }
  output_meassure[rep] <- sd(plantcomm$sz)

}


 LHS_final <- as.data.frame(LHS_param)
 LHS_final$output_meassure <- output_meassure

 # 

# write.csv(LHS_final, file = paste("LHS_sampling_results_",
#                                  n_reps,
#                                  "reps_", seed_value,"seed", ".csv",
#                                  sep = ""))

cor(LHS_final$output_meassure, LHS_final[,1:5], )

cor.test(LHS_final$output_meassure, LHS_final$max_initial_sz)

#
#
# png("correlation_plots.png",   width = 18, height = 18, units = "cm",
#     res = 300)
# 
# par(mfrow=c(2,3))
# with(LHS_final, {
#   plot(pop_density, output_meassure,
#      ylab=bquote(sigma[S]), xlab = "D" , pch=19)
#   grid()
#   points(pop_density, output_meassure, pch=19)
#   text(x = 0.02, 0.3, "A", cex=1.5, col="red")
# 
#   plot(world_reachability, output_meassure,
#        ylab=bquote(sigma[S]), xlab = bquote(kappa) , pch=19)
#   grid()
#   points(world_reachability, output_meassure, pch=19)
#   text(x = 0.1, 0.3, "B", cex=1.5, col="red")
# 
#   plot(max_initial_sz, output_meassure,
#        ylab=bquote(sigma[S]), xlab = bquote(s[0]) , pch=19)
#   grid()
#   points(max_initial_sz, output_meassure, pch=19)
#   text(x = 0.65, 0.3, "C", cex=1.5, col="red")
# 
#   plot(comp_symmetry, output_meassure,
#        ylab=bquote(sigma[S]), xlab = bquote(theta) , pch=19)
#   grid()
#   points(comp_symmetry, output_meassure, pch=19)
#   text(x = 0.15, 0.3, "D", cex=1.5, col="red")
# 
#   plot(max_growing_size, output_meassure,
#        ylab=bquote(sigma[S]), xlab = bquote(theta) , pch=19)
#   grid()
#   points(max_growing_size, output_meassure, pch=19)
#   text(x = 0.15, 0.3, "max_growing_size", cex=1.5, col="red")
# 
#   plot(max_growth_rate, output_meassure,
#        ylab=bquote(sigma[S]), xlab = bquote(theta) , pch=19)
#   grid()
#   points(max_growth_rate, output_meassure, pch=19)
#   text(x = 0.15, 0.3, "max_growth_rate", cex=1.5, col="red")
# })
#  dev.off()
