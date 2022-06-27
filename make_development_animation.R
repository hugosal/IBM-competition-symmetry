# script to make a gif animation of the development of a plant pop given a 
# set of population attribute parameters
gc()

source("circles_area_of_overlap.R")



library(animation)



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

ws <- 20 # world size 

timesteps <- 20 # length of each run

# The  number of plants of a population that would have just enough resources is
intermediate_pop <- 16

# The maximum size of plants such that they have just enough resources is
max_S <- ws/(sqrt(intermediate_pop) * 2)  # ?????????asumiendo que es un numero cuadrado!!!!!!!!!!!!!!!!!!!!!!!


world_reachablity <- 5
max_initial_size <- 2
n_overcrowding_plants <- 5
theta <-  1


seed_value <-31
set.seed(seed_value)

max_grwth_rt <- 0.1

y_lim_coefvar <- 0.5


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
  ft = rep(1, initial.n), 
  sz = runif(initial.n, min = 0.5, max = max_initial_size))

# modify plant coordinates accoding to spatial disarrangement parameter

random_angles <- runif(nrow(plantcomm), 
                       min = 0, max = 2*pi)

plantcomm$x <- (plantcomm$x + cos(random_angles) * world_reachablity) %% ws
plantcomm$y <- (plantcomm$y + sin(random_angles) * world_reachablity) %% ws

# plot_plantcomm(plantcomm, numbers = TRUE, ws = ws)
a <- (4 * max_grwth_rt)/ max_S
b <- (4 * max_grwth_rt)/ (max_S**2)
times <- c(0)
mean_coef_compt <-c(NA)
percet_competing_plants <-c(NA)
coef_vars <- c(sd(plantcomm$sz)/mean(plantcomm$sz))


saveGIF({
  for (t in 1:timesteps){
    layout(mat = matrix(c(1,1,3,2), 
                        nrow = 2, 
                        ncol = 2),
           heights = c(2, 1, 1),
           widths = c(2, 1, 1))
    
    plot_plantcomm(plantcomm, circle = T, numbers = 1:nrow(plantcomm), 
                   ws = ws, main=bquote(t==.(t)))
    
    sign_ovrcrwd <- if (n_overcrowding_plants > 0) "+" else if(n_overcrowding_plants < 0) "-" else " "
    
    mtext(text = bquote("Overcrowding "== ~ .(paste(sign_ovrcrwd, abs(n_overcrowding_plants), sep="")) ~","~ 
                        kappa==.(round(world_reachablity, 2))~","~
                        s[0]==.(round(max_initial_size, 2))~","),
          side = 1, line = 2)
    mtext(text = bquote(theta==.(round(theta, 2))~","~
                        S[max]==.(round(max_S, 3))~","~
                        M[grt]==.(round(max_grwth_rt,2))),
          side = 1, line = 4)
    par(mar=c(5, 4, 2, 2) + 0.1)
    hist(plantcomm$sz, xlab="Plant size", freq = T, main="",
         ylim = c(0, initial.n ), 
         breaks = seq(from =0, to = max_S * 1.2, length.out = 10  ))
    
    print(paste( "t ", t))
    par(mar=c(5, 4, 4, 4) + 0.1)
    plot(times, coef_vars, xlim = c(0, timesteps), ylim =c(0, y_lim_coefvar), xlab="t",
         ylab="Coefficient of Variation", lwd=2, type="l", main="")
    
    par(new = TRUE)
    plot(times, mean_coef_compt, type = "l", axes = FALSE, bty = "n", lwd=2,
         xlab = "", ylab = "", col="darkblue", ylim=c(0,1), xlim=c(0,timesteps))

    axis(side=4, at = seq(0, 1, length = 5 ) )
    mtext( bquote("Competition ("~phantom(bar(c))~","~phantom(" %")~ ")" ), 
           side=4, line=3, cex=1.2)
    mtext( bquote(phantom("Competition (")~bar(c)~phantom(", "~"%"~")")),
           side=4, line=3, col="darkblue", cex=1.2)
    mtext( bquote(phantom("Competition ("~bar(c))~phantom(", ")~"%"~ phantom(")")),
           side=4, line=3, col="magenta", cex=1.2)
    
    lines(times, percet_competing_plants, lwd=2, col="magenta")
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    
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
    
    # this optimization trick is to compute first intersections
    # with most circles, then subsets of those groups will be extracted
    # from memory instead of computed again
    most_to_least_crowded <- order(n_in_groups, decreasing = T)
    
    intersect_to_compute <- (1:n)[most_to_least_crowded]
    intersect_to_compute <- unname(intersect_to_compute[n_in_groups[most_to_least_crowded] > 1])
    # compute unbound growth rate, if this plant can grow further, compute
    # intersections areas. Otherwise, it will stay the same size, computing areas
    # will be irrelevant and thus it is omitted
    
    #cant_grow_more <- which(plantcomm$sz >= (max_S * ind_s_var))
    
    # initiate vector of resources obtained with the area of effect of each plant
    resources_obtained_p_ind <- plantcomm$sz**2 * pi
    
    intersections_list_memory <-list()
    
    for(j in intersect_to_compute){
      
      # if (j %in% cant_grow_more){
      #   if((a[j] * plantcomm$sz[j]) - (b[j] * plantcomm$sz[j]**2) > 1e-6){
      #     stop("Algo mal con optimizacion de crecimiento")
      #   }
      #   next
      # }
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
          
          intersections_list_memory[[length(intersections_list_memory) + 1]]<- list(
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

    times <- c(times, t)
    coef_vars <- c(coef_vars, sd(plantcomm$sz)/mean(plantcomm$sz))
    mean_coef_compt <- c(mean_coef_compt, mean(competition_effect))
    percet_competing_plants <- c(percet_competing_plants, sum(n_in_groups>1)/length(n_in_groups))
    
    
  }
  
}, movie.name = "example_develop_symmetryc.gif", interval = 1, 
ani.width = 1500, ani.height = 1000, clean = TRUE,  ani.res = 150)
    