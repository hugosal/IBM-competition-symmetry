#parametros a probar

# 1  estructura espacial
# 2 variacion inicial en tamanios
# 3 densisdad de plantas
# 4 asimetria
# 5 tasa de crecimiento

# en este quiero probar las condiciones que dice Weiner1990 sobre el efecto
# de las interactiones y su simetria en la variacion de 

source("circles_area_of_overlap.R")
#source("get_ppis.R")

seed_value <-23
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

plot_plantcomm <- function(com, numbers = FALSE, main="", ws = NULL, circle = TRUE){
  
  text_box <- function(x, y, labels, cex){
    sw   <- strwidth(labels)
    sh   <- strheight(labels)
    frsz <- 0
    
    if (strwidth(labels)>strwidth("1")){
      text(x, y, labels, font=2, col="white", cex=cex)
      text(x, y, labels, cex=cex*0.98, font=1)
    }else{
      text(x, y, labels, font=2, col="white", cex=cex)
      text(x, y, labels, cex=cex*0.98, font=1) # con  el problema del que no se veia bien
      # aca cambia el 0.98
    }
  }
  
  colores <- col2rgb(3:6)
  colores <- apply(colores, 2, FUN = function(x)rgb(x[1]/255,
                                                    x[2]/255,
                                                    x[3]/255,
                                                    alpha = 0.5)  )
  plot(com$x, com$y, col=com$ft, type="n", main=main, xlab="", ylab="",
       yaxt="n", xaxt="n", yaxs="i", xaxs="i", asp=1,
       xlim = c(0,ws),
       ylim = c(0,ws))
  
  
  
  for (j in 1:nrow(com)){
    
    if (circle){
      at_border <- c(com[j, 1] + com[j, 4] > ws | com[j, 1] - com[j, 4] < 0,
                     com[j, 2] + com[j, 4] > ws | com[j, 2] - com[j, 4] < 0)
    }else{
      at_border <- c(abs(com[j, 1]- ws) < 1e-4 | com[j, 1] < 1e-4,
                     abs(com[j, 2]- ws) < 1e-4 | com[j, 2] < 1e-4)
      
    }
    
    
    
    if (any(at_border)){
      for(k in list(c(0,0),
                    c(ws,0),
                    c(-ws,0),
                    c(0,ws),
                    c(0,-ws),
                    c(ws,ws),
                    c(-ws,ws)
      )){
        if (circle){
          plotrix::draw.circle(x = com$x[j]+k[1], y = com$y[j]+k[2], radius = com$sz[j], 
                               col = colores[com$ft[j]])
        }else{
          points(x = com$x[j]+k[1], y = com$y[j]+k[2], col=com$ft, pch=19)
        }
        
      }}
    else{
      if(circle){
        plotrix::draw.circle(x = com$x[j], y = com$y[j], radius = com$sz[j], 
                             col = colores[com$ft[j]])
      }else{
        points(com$x, com$y, col=com$ft, pch=19)
      }
    }
    offst <- ws*0.97
    letter_sz <- 0.9
    
    
    if (numbers){
      
      if (sum(at_border) == 0){
        if (circle){
          text_box(x = com$x[j], y = com$y[j], labels = j, cex=letter_sz)
        }else{
          text_box(x = com$x[j]+(ws-offst), y = com$y[j]+(ws-offst), 
                   labels = j, cex=letter_sz)
        }
        
      }else if (sum(at_border)==2){
        if (circle){text_box(x = ws-offst, y = ws-offst, labels = j, cex=letter_sz)
          text_box(x = ws-offst, y = offst, labels = j, cex=letter_sz)
          text_box(x = offst, y = offst, labels = j, cex=letter_sz)
          text_box(x = offst, y = ws-offst, labels = j, cex=letter_sz)
        }else{
          text_box(x = ws-offst, y = ws-offst, labels = j, cex=letter_sz)
          text_box(x = ws-offst, y = offst, labels = j, cex=letter_sz)
          text_box(x = offst, y = offst, labels = j, cex=letter_sz)
          text_box(x = offst, y = ws-offst, labels = j, cex=letter_sz)
          
        }
        
      }else if (at_border[2]){
        if (circle){
          text_box(x = com$x[j], ws-offst, labels = j, cex=letter_sz)
          text_box(x = com$x[j], offst, labels = j, cex=letter_sz)
        }else{
          text_box(x = com$x[j]+ws-offst, ws-offst, labels = j, cex=letter_sz)
          text_box(x = com$x[j]+ws-offst, offst, labels = j, cex=letter_sz)
        }
        
        
      }else if (at_border[1]){
        
        if (circle){
          text_box(x = ws-offst , y = com$y[j], labels = j, cex=letter_sz)
          text_box(x = offst , y = com$y[j], labels = j, cex=letter_sz)
          
        }else{
          text_box(x = ws-offst , y = com$y[j]+ws-offst, labels = j, cex=letter_sz)
          text_box(x = offst , y = com$y[j]+ws-offst, labels = j, cex=letter_sz)
          
        }
        
      }
    }}
  
  if (! is.null(ws)){
    polygon(c(-10, ws+50, ws+50,-1), c(0,0,-100,-100), col="white", border=NA)
    
    polygon(c(0, 0, -100,-100), c(0,ws+5,-20,-20), col="white" , border=NA)
    
    polygon(c(ws, ws, ws+5,ws+5), c(-60,ws+50,ws+50,-60), col="white", border=NA)
    
    polygon(c(-60, ws+50, ws+50,-10), c(ws,ws,ws+50,ws+50), col="white", border=NA)
    
    lines(x = c(0, ws), y = c(0, 0)  )
    lines(x = c(0, 0), y = c(ws, 0)  )
    lines(x = c(ws, 0), y = c(ws, ws)  )
    lines(x = c(ws, ws), y = c(ws, 0)  )
  }
}

generate_initial_points <- function(N, ws){
  if (sqrt(N) == round(sqrt(N)) ){
    n_per_row <- sqrt(N)
    grid_poinst <- seq(from = 0, to = ws, 
                       by = ws/(n_per_row))[1:n_per_row]
    x <- rep(grid_poinst, n_per_row)
    y <- as.vector(sapply(grid_poinst, function(x) rep(x, n_per_row)))
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
      return(coords / 20 * ws) # Scaled to fit ws
    }else{
      stop("Configuration not found yet")}
  }
}

uniform_LHS_sample_from_range <- function(lower, upper, n_samples){
  limits <- seq(from = lower, to = upper, length.out = n_samples + 1)
  sapply(1:n_samples, function(x){ runif(n = 1, 
                                         min = limits[x], 
                                         max = limits[x + 1])})
}

ws <- 20 # world size 

timesteps <- 50   # length of each run

# slf_thinning_limit <- 0.1 # if competition effect is stronger or equal to this
#                           plants will die

n_reps <- 900 # number of LHS samples

  
# The  number of plants of a population that would have just enough resources is
intermediate_pop <- 16

# The maximum size of plants such that they have just enough resources is
max_S <- ws/(sqrt(intermediate_pop) * 2)  # ¡¡¡¡¡¡¡¡¡asumiendo que es un numero cuadrado!!!!!!!!!!!!!!!!!!!!!!!

# cada fila va a ser una fila de parametros a usar


###### lo que este dentro de esta seccion es solo para hacer el arbol complicandose
vector_dec_from_bin <- numeric(2**4)
binarynames <- sapply(X = 1:length(vector_dec_from_bin),
                                        FUN = function(x){
                             bin <- paste(rev(as.integer(intToBits(x))), collapse="") 
                             substr(bin, start = nchar(bin) - 4 + 1 , stop = nchar(bin))})
var_to_include <- lapply(binarynames, function(x) strsplit(x, "")[[1]]==1)

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
                                                           n_samples = n_reps)), # wrl reach,
                      sample(uniform_LHS_sample_from_range(lower = 0.5, 
                                                           upper = max_S/2, 
                                                           n_samples = n_reps)), # max initial size,
                      sample(floor(uniform_LHS_sample_from_range(lower = -10, 
                                                                 upper = 15,  # lower and upper chosen for rounding error
                                                                 n_samples = n_reps))), # number of overcrowding plants
                      sample(c(uniform_LHS_sample_from_range(lower = 0,
                                                             upper = 1,  
                                                             n_samples = n_reps/2),
                               uniform_LHS_sample_from_range(lower = 1,  
                                                             upper = 100, 
                                                             n_samples = n_reps/2))) # competition asymmetry parameter, divided en two ranges
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
    
    critical_distance <- max_S * 2
    
    mean_potential_competitors <- mean(apply(dists < critical_distance, MARGIN = 1, sum) - 1)
    sd_potential_competitors <- sd(apply(dists < critical_distance, MARGIN = 1, sum) - 1)
    
    # record data at time 0: before plants start interacting
    results_over_time[counter, ] <- c(re, 
                                      0, 
                                      mean(plantcomm$sz), 
                                      sd(plantcomm$sz), 
                                      sum(pi * plantcomm$sz**2), 
                                      NA,
                                      mean_potential_competitors,
                                      sd_potential_competitors)
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
  	
  	################### Mortality by self thinning. Plants will die if the
  	# effect of competition is stronger than a thresholf
  	
  	#dead_by_self_thinning <- competition_effect <= slf_thinning_limit
  	
  	#if (any(dead_by_self_thinning)) {
  	#  print("overcrowding dead")
   # 	  plantcomm <- plantcomm[!dead_by_self_thinning, ]
   # 	  rownames(plantcomm) = 1:nrow(plantcomm)
   # 	  a <- a[!dead_by_self_thinning]
   # 	  b <- b[!dead_by_self_thinning]}

  
  	
  	####################

  	results_over_time[counter, ] <- c(re, 
  	                                  t, 
  	                                  mean(plantcomm$sz), 
  	                                  sd(plantcomm$sz), 
                                      sum(pi * plantcomm$sz**2), 
                                      mean(competition_effect),
  	                                  mean_potential_competitors,
  	                                  sd_potential_competitors)
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

