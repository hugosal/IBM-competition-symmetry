# script to make a gif animation of the development of a plant pop given a 
# set of population attribute parameters
source("circles_area_of_overlap.R")
seed_value <-21
library(animation)
set.seed(seed_value)
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
      return(coords / 20 * ws) # Scaled to fit ws
    }else{
      stop("Configuration not found yet")}
  }
}


ws <- 40 # world size 

timesteps <- 30   # length of each run

world_reachablity <- 3
max_initial_size <- 3
pop_density <- 30/(ws**2)
theta <- 5
max_S <- 5
max_grwth_rt <- 0.5

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

a <- ( 4 * max_grwth_rt)/max_S
b <- ( 4 * max_grwth_rt)/(max_S**2)


saveGIF({
  for (t in 1:timesteps){
    par(mfrow=c(1,2))
    plot_plantcomm(plantcomm, numbers = TRUE, ws = ws, main=bquote(t==.(t)))
    
    mtext(text = bquote("D"==over(.(round(pop_density*ws**2)) , .(ws)^2) ~","~ 
                        kappa==.(round(world_reachablity),2)~","~
                        s[0]==.(round(max_initial_size, 2))~","),
          side = 1, line = 2)
    mtext(text = bquote(theta==.(round(theta),2)~","~
                        S[max]==.(round(max_S),2)~","~
                        M[grt]==.(round(max_grwth_rt,2))),
          side = 1, line = 4)
    hist(plantcomm$sz, xlab="Plant size", freq = T, main="",
         ylim = c(0, initial.n ),  
         breaks = seq(from =0, to = max_S, length.out = 10  ))
    
    
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
        
        if (length(names_original_plantcomm) > 9){
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
          #print("usando hack de memoria")
          
          intersections_list_memory[[length(intersections_list_memory) + 1]]<- list(
            "names_plant" = names_original_plantcomm,
            "intersections" = intersections)
        }
      }
    }
    # plants increase in size asymptotically
    
    
    if(!all(resources_obtained_p_ind <= pi * plantcomm$sz**2)){
      stop("Something wrong with resource partitions")
    }
    
    # 1 in competition_effect indicates a plants without competition
    competition_effect <- resources_obtained_p_ind / (pi * plantcomm$sz**2)
    
    # logistic growth rate according to 10.1016/S0304-3800(98)00182-3
    
    indgr <-  (competition_effect * a * plantcomm$sz) - (b * plantcomm$sz**2)
    
    indgr[indgr < 0] <- 0
    
    plantcomm$sz = plantcomm$sz + indgr
    
  }
  
}, movie.name = "example_develop_1.gif", interval = 1, 
ani.width = 960, ani.height = 480)
    