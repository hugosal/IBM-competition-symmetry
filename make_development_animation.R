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
       yaxt="n", xaxt="n", yaxs="i", xaxs="i", 
       xlim = c(0,ws),
       ylim = c(0,ws), asp=1)
  
  
  
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
    polygon(c(-1, ws+5, ws+5,-1), c(0,0,-10,-10), col="white", border=NA)
    
    polygon(c(0, 0, -10,-10), c(0,ws+5,-2,-2), col="white" , border=NA)
    
    polygon(c(ws, ws, ws+5,ws+5), c(-6,ws+5,ws+5,-6), col="white", border=NA)
    
    polygon(c(-6, ws+5, ws+5,-1), c(ws,ws,ws+5,ws+5), col="white", border=NA)
    
    lines(x = c(0, ws), y = c(0, 0)  )
    lines(x = c(0, 0), y = c(ws, 0)  )
    lines(x = c(ws, 0), y = c(ws, ws)  )
    lines(x = c(ws, ws), y = c(ws, 0)  )
  }
}

plot_plantcomm_deads <- function(com, numbers = FALSE, main="", ws = NULL, circle = TRUE){
  
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
  
  colores <- col2rgb("gray")
  colores <- apply(colores, 2, FUN = function(x)rgb(x[1]/255,
                                                    x[2]/255,
                                                    x[3]/255,
                                                    alpha = 0.5)  )
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
    polygon(c(-1, ws+5, ws+5,-1), c(0,0,-10,-10), col="white", border=NA)
    
    polygon(c(0, 0, -10,-10), c(0,ws+5,-2,-2), col="white" , border=NA)
    
    polygon(c(ws, ws, ws+5,ws+5), c(-6,ws+5,ws+5,-6), col="white", border=NA)
    
    polygon(c(-6, ws+5, ws+5,-1), c(ws,ws,ws+5,ws+5), col="white", border=NA)
    
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


ws <- 20 # world size 

timesteps <- 30   # length of each run

slf_thinning_limit <- -1 # if competition effect is stronger or equal to this
#                           plants will die


# The  number of plants of a population that would have just enough resources is
intermediate_pop <- 16

# The maximum size of plants such that they have just enough resources is
max_S <- ws/(sqrt(intermediate_pop) * 2)  # ¡¡¡¡¡¡¡¡¡asumiendo que es un numero cuadrado!!!!!!!!!!!!!!!!!!!!!!!


world_reachablity <- 10
max_initial_size <- 0.5
n_overcrowding_plants <- 10
theta <-  40


seed_value <-32
set.seed(seed_value)

max_grwth_rt <- 0.1
indiviudal_var_growth_rate <- 0

y_lim_coefvar <- 0.5


initial.n <- intermediate_pop + n_overcrowding_plants
config_found <- c(4,5,9,10,13,16,17,20,25,29,36,49,64)
mean_coef_compt <-c(NA)
percet_competing_plants <-c(NA)

if (! initial.n %in% config_found){
  print("adj")
  initial.n <- config_found[which.min(abs(initial.n - config_found))]
}
print(unname(initial.n))

# set plants initial coordinates
coordinates <- generate_initial_points(N = initial.n, ws = ws)

# set plants initial coordinates
coordinates <- generate_initial_points(N = initial.n, ws = ws)

plantcomm = data.frame(
  x = coordinates[, 1],
  y = coordinates[, 2],
  ft = 1,
  sz = runif(initial.n, min = 0.5, max = max_initial_size))

dead_plants <- data.frame(x = numeric(),
                          y = numeric(),
                          ft = numeric(),
                          sz = numeric())
# modify plant coordinates accoding to spatial disarrangement parameter

random_angles <- runif(nrow(plantcomm), 
                       min = 0, max = 2*pi)

plantcomm$x <- (plantcomm$x + cos(random_angles) * world_reachablity) %% ws
plantcomm$y <- (plantcomm$y + sin(random_angles) * world_reachablity) %% ws

# plot_plantcomm(plantcomm, numbers = TRUE, ws = ws,
#                main=paste("density", round(pop_density,2),
#                           "reachabil", round(world_reachablity,2)))

ind_s_var <- runif(min = 1, 
                   max = 1 + indiviudal_var_growth_rate,
                   n = nrow(plantcomm))
ind_grwth_rt <- runif(min = 1,
                      max = 1 + indiviudal_var_growth_rate, 
                      n = nrow(plantcomm))

a <- (4 * (max_grwth_rt * ind_grwth_rt))/(max_S * ind_s_var)
b <- (4 * (max_grwth_rt * ind_grwth_rt))/((max_S * ind_s_var)**2)

times <- c(0)
coef_vars <- c(sd(plantcomm$sz)/mean(plantcomm$sz))



saveGIF({
  for (t in 1:timesteps){
    layout(mat = matrix(c(1,1,3,2), 
                        nrow = 2, 
                        ncol = 2),
           heights = c(2, 1, 1),
           widths = c(2, 1, 1))
    plot_plantcomm(plantcomm, numbers = TRUE, ws = ws, main=bquote(t==.(t)))
    if (nrow(dead_plants)>0){
    plot_plantcomm_deads(dead_plants, numbers = TRUE, ws = ws, main=bquote(t==.(t)))
    }  
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
    plot(times, coef_vars, xlim = c(0,30), ylim =c(0, y_lim_coefvar), xlab="t",
         ylab="Coefficient of Variation", lwd=2, type="l", main="")
    
    par(new = TRUE)
    plot(times, mean_coef_compt, type = "l", axes = FALSE, bty = "n", lwd=2,
         xlab = "", ylab = "", col="darkblue", ylim=c(0,1), xlim=c(0,30))

    axis(side=4, at = seq(0, 1, length = 5 ) )
    mtext( bquote("Competition ("~phantom(bar(c))~","~phantom(" %")~ ")" ), 
           side=4, line=3, cex=1.2)
    mtext( bquote(phantom("Competition (")~bar(c)~phantom(", "~"%"~")")),
           side=4, line=3, col="darkblue", cex=1.2)
    mtext( bquote(phantom("Competition ("~bar(c))~phantom(", ")~"%"~ phantom(")")),
           side=4, line=3, col="magenta", cex=1.2)
    
    
    lines(times, percet_competing_plants, lwd=2, col="magenta")
    
    
    
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
    
    ################### Mortality by self thinning. Plants will die if the
    # effect of competition is stronger than a thresholf
    
    dead_by_self_thinning <- competition_effect <= slf_thinning_limit
    
    if (any(dead_by_self_thinning)) {
      
      dead_plants <- rbind(dead_plants, plantcomm[dead_by_self_thinning, ])
      plantcomm <- plantcomm[!dead_by_self_thinning, ]
      rownames(plantcomm) = 1:nrow(plantcomm)
      a <- a[!dead_by_self_thinning]
      b <- b[!dead_by_self_thinning]
    }
    times <- c(times, t)
    coef_vars <- c(coef_vars, sd(plantcomm$sz)/mean(plantcomm$sz))
    mean_coef_compt <- c(mean_coef_compt, mean(competition_effect))
    percet_competing_plants <- c(percet_competing_plants, sum(n_in_groups>1)/length(n_in_groups))
    
    
  }
  
}, movie.name = "example_develop_1.gif", interval = 1, 
ani.width = 1500, ani.height = 1000, clean = TRUE,  ani.res = 150)
    