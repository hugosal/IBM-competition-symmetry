library(spatstat)

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

ws <- 20 # world size 
N <- 25
points_initial <- generate_initial_points(N = N, ws = ws)

spatstats_initial_points <- ppp(x = points_initial[,1], 
                                y = points_initial[,2], c(0, ws), c(0, ws))

for( world_reachablity in c(0, 0.5,  1, 2, 5, 10, 15)){
  points_after_moving <- points_initial
  
  set.seed(26)
  random_angles <- runif(nrow(points_after_moving), 
                         min = 0, max = 2*pi)

  points_after_moving[, 1] <- (points_after_moving[, 1] + (cos(random_angles) * world_reachablity)) %% ws
  points_after_moving[, 2] <- (points_after_moving[, 2] + (sin(random_angles) * world_reachablity)) %% ws
  
  
  spatstats_second_points <- ppp(x = points_after_moving[,1], 
                                 y = points_after_moving[,2], c(0, ws), c(0, ws))
  
  png(paste("randomness_example_", world_reachablity, ".png", sep = ""),   
            width = 24, height = 12, units = "cm",
      res = 300)
  par(mfrow=c(1,2))
  par(mar=c(3, 3, 3, 3))
  plot(points_after_moving[,1],points_after_moving[,2],
       main = "", xlab="",ylab="", 
       yaxt="n", xaxt="n" , col ="black", pch=19)
  mtext(text = bquote(N==.(N)~","~ kappa==.(world_reachablity)), side = 3,
        line = - 2,
        outer = TRUE)
  par(mar=c(5, 4, 4, 2) + 0.1)
  Kest_out <-  Kest(spatstats_second_points, correction = "none" )
  plot(Kest(spatstats_second_points, correction = "none" ), 
       main = "")
  
  dev.off()
  
}

points_after_moving <- points_initial

  points_after_moving[, 1] <- c(rnorm(n = 20, mean = 15, sd = 2 ),
                                rnorm(n = 20, mean = 0, sd = 1 )) %% ws
  points_after_moving[, 2] <- c(rnorm(n = 20, mean = 15, sd = 3 ),
                                rnorm(n = 20, mean = 0, sd = 2 )) %% ws
  
  
  spatstats_second_points <- ppp(x = points_after_moving[,1], 
                                 y = points_after_moving[,2], c(0, ws), c(0, ws))
  

  par(mfrow=c(1,2))
  par(mar=c(3, 3, 3, 3))
  plot(points_after_moving[,1],points_after_moving[,2],
       main = "", xlab="",ylab="", 
       yaxt="n", xaxt="n" , col ="black", pch=19, xlim = c(0, ws), ylim = c(0, ws))
  par(mar=c(5, 4, 4, 2) + 0.1)
  plot(Kest(spatstats_second_points, correction = "none" ), 
       main = "")

  


