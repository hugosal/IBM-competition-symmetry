#parametros a probar

# 1  estructura espacial
# 2 variacion inicial en tamanios
# 3 densisdad de plantas
# 4 asimetria
# 5 tasa de crecimiento

# en este quiero probar las condiciones que dice Weiner1990 sobre el efecto
# de las interactiones y su simetria en la variacion de 

library(CirclesIntersections)
source("auxiliary_functions.R")

#source("get_ppis.R")

seed_value <-26
set.seed(seed_value)


uniform_LHS_sample_from_range <- function(lower, upper, n_samples){
  limits <- seq(from = lower, to = upper, length.out = n_samples + 1)
  sapply(1:n_samples, function(x){ runif(n = 1, 
                                         min = limits[x], 
                                         max = limits[x + 1])})
}

ws <- 20 # world size 

timesteps <- 50   # length of each run

n_reps <- 8000 # number of LHS samples
  
# The  number of plants of a population that would have just enough resources is
intermediate_pop <- 16

# The maximum size of plants such that they have just enough resources is
max_S <- ws/(sqrt(intermediate_pop) * 2)  # ¡¡¡¡¡¡¡¡¡asumiendo que es un numero cuadrado!!!!!!!!!!!!!!!!!!!!!!!


###### lo que este dentro de esta seccion es solo para hacer el arbol complicandose
# vector_dec_from_bin <- numeric(2**4)
# binarynames <- sapply(X = 1:length(vector_dec_from_bin),
#                                         FUN = function(x){
#                              bin <- paste(rev(as.integer(intToBits(x))), collapse="") 
#                              substr(bin, start = nchar(bin) - 4 + 1 , stop = nchar(bin))})
# var_to_include <- lapply(binarynames, function(x) strsplit(x, "")[[1]]==1)

# for(vars_t_inc in seq_along(var_to_include)) {
 # print(paste('##################', vars_t_inc))

 # number_var_rthis_run <- sum(var_to_include[[vars_t_inc]])
  consider <- c(world_reachability = TRUE, max_initial_sz = TRUE, 
                n_overcrowding_plants = TRUE, comp_symmetry = TRUE)

# consider <- c(world_reachability = var_to_include[[vars_t_inc]][1], 
#              max_initial_sz = var_to_include[[vars_t_inc]][2], 
#              n_overcrowding_plants = var_to_include[[vars_t_inc]][3],
#              comp_symmetry = var_to_include[[vars_t_inc]][4], 
#              max_growth_rate = FALSE, indiviudal_var_growth_rate = FALSE)


null_values <- c(world_reachability = 0, max_initial_sz = 0.5, n_overcrowding_plants = 0, 
  comp_symmetry = 1)

LHS_param <- matrix(c(sample(uniform_LHS_sample_from_range(lower = 0, 
                                                           upper = 10, 
                                                           n_samples = n_reps)), # world reach kappa parameter
                      sample(uniform_LHS_sample_from_range(lower = 0.5, 
                                                           upper = max_S/2, 
                                                           n_samples = n_reps)), # max initial size,
                      sample(floor(uniform_LHS_sample_from_range(lower = -10, 
                                                                 upper = 15,  # lower and upper chosen for rounding error
                                                                 n_samples = n_reps))), # number of overcrowding plants
                      sample(c(uniform_LHS_sample_from_range(lower = 0,
                                                             upper = 1,  
                                                             n_samples = n_reps %/% 2),
                               uniform_LHS_sample_from_range(lower = 1,  
                                                             upper = 100, 
                                                             n_samples = n_reps %/% 2))) # competition asymmetry parameter, divided en two ranges
                      ),
                    nrow = n_reps)

colnames(LHS_param) <- c("world_reachability", 
                         "max_initial_sz", 
                         "n_overcrowding_plants",
                         "comp_symmetry")


for (con in seq_along(consider)){
  if ( !consider[con] ) {
    LHS_param[, con ] <- rep(null_values[con], n_reps)
  } 
}


results_over_time <- data.frame(re = numeric(n_reps * (timesteps + 1)),
                      t = numeric(n_reps * (timesteps + 1)),
                      size_mean = numeric(n_reps * (timesteps + 1)),
                      size_sd = numeric(n_reps * (timesteps + 1)), 
                      total_biomass = numeric(n_reps * (timesteps + 1)),
                      mean_competiton = numeric(n_reps * (timesteps + 1)),
                      mean_no_competitors = numeric(n_reps * (timesteps + 1)),
                      sd_no_competitors = numeric(n_reps * (timesteps + 1)))
counter <- 1
for (re in 1:nrow(LHS_param)){
  print(paste( "re ", re))
  # extract value of model parameters
  world_reachablity <- LHS_param[re, 1]
  max_initial_size <- LHS_param[re, 2]
  n_overcrowding_plants <- LHS_param[re, 3]
  theta <- LHS_param[re, 4]
  
  max_grwth_rt <- 0.1
  
  # set population size
  initial.n <- intermediate_pop + n_overcrowding_plants
  config_found <- c(4,5,9,10,13,16,17,20,25,29,36,49,64)
  
  if (! initial.n %in% config_found){
    print("adj")
    initial.n <- config_found[which.min(abs(initial.n - config_found))]
  }
  print(unname(initial.n))

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
  
  plantcomm$x <- (plantcomm$x + cos(random_angles) * world_reachablity) %% ws
  plantcomm$y <- (plantcomm$y + sin(random_angles) * world_reachablity) %% ws
  
  # plants have variation in their groth speed, which can be higher tan the preset
  # according to indiviudal_var_growth_rate
    
    a <- (4 * max_grwth_rt)/ max_S
    b <- (4 * max_grwth_rt)/ (max_S**2)

    # compute distance between points
    
    dx = as.matrix(dist(plantcomm$x,diag=TRUE,upper=TRUE))
    dx[dx > ws/2] = ws - dx[dx > ws/2]   		### this wraps the interaction effects
    dy = as.matrix(dist(plantcomm$y,diag=TRUE,upper=TRUE))
    dy[dy > ws/2] = ws - dy[dy > ws/2]
    dists = sqrt(dx^2 + dy^2)
    
    n <- nrow(plantcomm)
    
    critical_distance <- outer(plantcomm$sz, plantcomm$sz, "+")   #max_S * 2
    neighbours <- apply(dists < critical_distance, MARGIN = 1, sum) - 1
    mean_neighbours <- mean(neighbours)
    sd_neighbours <- sd(neighbours)
    
    # record data at time 0: before plants start interacting
    results_over_time[counter, ] <- c(re, 
                                      0, 
                                      mean(plantcomm$sz), 
                                      sd(plantcomm$sz), 
                                      sum(pi * plantcomm$sz**2), 
                                      NA,
                                      mean_neighbours,
                                      sd_neighbours)
    counter <- counter + 1
  
  for (t in 1:timesteps){
    #plot_plantcomm(plantcomm, numbers = TRUE, ws = ws, main=bquote(t==.(t)))
    #print(paste( "t ", t))
    	radius_sum <- outer(X = plantcomm$sz, 
  	                    Y = plantcomm$sz, FUN = "+")
    
  	# which plant area of effect overlap
  	disteffect <- dists < radius_sum
  	n_in_groups <- apply(disteffect, MARGIN = 1, FUN = sum)
  	
  	# this optimization trick is to compute first intersections
  	# with most circles, then subsets of those groups will be extracted
  	# from memory instead of computed again
  	most_to_least_crowded <- order(n_in_groups, decreasing = T)

  	intersect_to_compute <- (1:n)[most_to_least_crowded]
  	intersect_to_compute <- unname(intersect_to_compute[n_in_groups[most_to_least_crowded] > 1])
  	# compute unbound growth rate, if this plant can grow further, compute
  	# intersections areas. Otherwise, it will stay the same size, computing areas
  	# will be irrelevant and thus it is omitted
  	
  	cant_grow_more <- which(plantcomm$sz >= max_S )
  		
  	# initiate vector of resources obtained with the area of effect of each plant
  	resources_obtained_p_ind <- plantcomm$sz**2 * pi
  
  	intersections_list_memory <-list()
  
  	for(j in intersect_to_compute){
  	  
  	  if (j %in% cant_grow_more){
  	    if((a * plantcomm$sz[j]) - (b * plantcomm$sz[j]**2) > 1e-6){
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
      	    #print("usando hack de memoria")
  
      	    intersections_list_memory[[length(intersections_list_memory) + 1]] <- list(
      	      "names_plant" = names_original_plantcomm,
      	      "intersections" = intersections)
      	    }
      	  }
    	  }
  	# plants increase in size asymptotically
  
  
  	if(any( (resources_obtained_p_ind -  (pi * plantcomm$sz**2)) > 1e-10   )  ){
  	  stop("Something wrong with resource partitions")
  	}
  	
  	# 1 in competition_effect indicates a plants without competition
  	competition_effect <- resources_obtained_p_ind / (pi * plantcomm$sz**2)
  	
  	# logistic growth rate according to 10.1016/S0304-3800(98)00182-3
  	
  	indgr <-  (competition_effect * a * plantcomm$sz) - (b * plantcomm$sz**2)

  	indgr[indgr < 0] <- 0
  	
  	plantcomm$sz = plantcomm$sz + indgr
  	
  	critical_distance <- outer(plantcomm$sz, plantcomm$sz, "+")   #max_S * 2
  	neighbours <- apply(dists < critical_distance, MARGIN = 1, sum) - 1
  	mean_neighbours <- mean(neighbours)
  	sd_neighbours <- sd(neighbours)
  	
  	results_over_time[counter, ] <- c(re, 
  	                                  t, 
  	                                  mean(plantcomm$sz), 
  	                                  sd(plantcomm$sz), 
                                      sum(pi * plantcomm$sz**2), 
                                      mean(competition_effect),
  	                                  mean_neighbours,
  	                                  sd_neighbours)
  	counter <- counter + 1

  }
}

LHS_final <- as.data.frame(LHS_param)

results_over_time <- cbind(results_over_time,
                           data.frame(matrix(ncol = length(colnames(LHS_final[1, 1:4])),
                                            nrow = nrow(results_over_time),
                                            dimnames = list(NULL, colnames(LHS_final[1, 1:4])))))

results_over_time$coef_var <- results_over_time$size_sd/results_over_time$size_mean

for (case in unique(results_over_time$re)){

  results_over_time[results_over_time$re==case , 9:12] <- LHS_final[case, 1:4]
}


# write.csv(results_over_time, file = paste("tree", number_var_rthis_run,"_variables_",paste(names(consider)[consider], collapse= "_"),
#  "_",n_reps, "_reps.csv", sep = ""))
# }
write.csv(results_over_time, 
          file = paste("IBM_res_", paste(names(consider)[consider], collapse= "_"),
                       "_", n_reps, "_reps_", ws, "ws_", timesteps, "tsteps_", 
                       seed_value, "_seed.csv", sep = ""))
print("Success")

