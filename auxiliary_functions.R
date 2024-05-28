

# This function centers xy coordinates such that they are 
# centered around center_x, center_y. This is used to 
# avoid the wrapping effect of the torus afecting computation of distances

center_world_arround <- function(center_xy, ws, coords){
  row_number <- which(rownames(coords)==as.character(center_xy))
  x_offset <- (ws/2) - coords$x[row_number]
  y_offset <- (ws/2) - coords$y[row_number]
  coords$x <- (coords$x + x_offset) %% ws
  coords$y <- (coords$y + y_offset) %% ws
  coords
}

# Function to generate the coordinates of N arranged points in a torus world of size ws
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

# Function to plot a plant population, number indicate if the number of each plant 
# is to be shown, circle indicates if the circles indicating the size of each 
# plants is to be shown.

plot_plantcomm <- function(com, numbers = NULL, main="", ws = NULL, circle = TRUE, 
                           col_circles = "#228B2280"){
  
  text_box <- function(x, y, labels, cex){
    sw   <- strwidth(labels)
    sh   <- strheight(labels)
    frsz <- 0
    
    if (strwidth(labels)>strwidth("1")){
      text(x, y, labels, font=2, col="white", cex=cex)
      text(x, y, labels, cex=cex*0.98, font=1)
    }else{
      text(x, y, labels, font=2, col="white", cex=cex)
      text(x, y, labels, cex=cex*0.98, font=1) 
    }
  }
  
  plot(com$x, com$y, col=com$ft, type="n", main=main, xlab="", ylab="",
       yaxt="n", xaxt="n", yaxs="i", xaxs="i",
       xlim = c(0,ws),
       ylim = c(0,ws), asp=1, bty="n")
  
  
  
  for (j in 1:nrow(com)){
    
    this_names <- if (!is.null(numbers)) numbers[j] else j
    
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
                               col = col_circles, border = NA)
        }else{
          points(x = com$x[j]+k[1], y = com$y[j]+k[2], col = com$ft, pch=19)
        }
        
      }}
    else{
      if(circle){
        plotrix::draw.circle(x = com$x[j], y = com$y[j], radius = com$sz[j], 
                             col = col_circles, border = NA)
      }else{
        points(com$x, com$y, col=com$ft, pch=19)
      }
    }
    offst <- ws*0.97
    letter_sz <- 0.9
    
    
    if (!is.null(numbers)){
      
      if (sum(at_border) == 0){
        if (circle){
          text_box(x = com$x[j], y = com$y[j], labels = this_names, cex = letter_sz)
        }else{
          text_box(x = com$x[j]+(ws-offst), y = com$y[j]+(ws-offst), 
                   labels = this_names, cex = letter_sz)
        }
        
      }else if (sum(at_border)==2){
        if (circle){text_box(x = ws-offst, y = ws-offst, labels = this_names, cex = letter_sz)
          text_box(x = ws-offst, y = offst, labels = this_names, cex = letter_sz)
          text_box(x = offst, y = offst, labels = this_names, cex = letter_sz)
          text_box(x = offst, y = ws-offst, labels = this_names, cex = letter_sz)
        }else{
          text_box(x = ws-offst, y = ws-offst, labels = this_names, cex = letter_sz)
          text_box(x = ws-offst, y = offst, labels = this_names, cex = letter_sz)
          text_box(x = offst, y = offst, labels = this_names, cex = letter_sz)
          text_box(x = offst, y = ws-offst, labels = this_names, cex = letter_sz)
          
        }
        
      }else if (at_border[2]){
        if (circle){
          text_box(x = com$x[j], ws-offst, labels = this_names, cex = letter_sz)
          text_box(x = com$x[j], offst, labels = this_names, cex = letter_sz)
        }else{
          text_box(x = com$x[j]+ws-offst, ws-offst, labels = this_names, cex = letter_sz)
          text_box(x = com$x[j]+ws-offst, offst, labels = this_names, cex = letter_sz)
        }
        
        
      }else if (at_border[1]){
        
        if (circle){
          text_box(x = ws-offst , y = com$y[j], labels = this_names, cex = letter_sz)
          text_box(x = offst , y = com$y[j], labels = this_names, cex = letter_sz)
          
        }else{
          text_box(x = ws-offst , y = com$y[j]+ws-offst, labels = this_names, cex = letter_sz)
          text_box(x = offst , y = com$y[j]+ws-offst, labels = this_names, cex = letter_sz)
          
        }
        
      }
    }}
  
  if (! is.null(ws)){
    polygon(c(-1000, ws+1000, ws+1000,-1000), c(0,0,-1000,-1000), col="white", border=NA)
    
    polygon(c(0, 0, -1000,-1000), c(0,ws+1000,-1000,-1000), col="white" , border=NA)
    
    polygon(c(ws, ws, ws+1000,ws+1000), c(-1000,ws+1000,ws+1000,-1000), col="white", border=NA)
    
    polygon(c(-1000, ws+1000, ws+1000,-1000), c(ws,ws,ws+1000,ws+1000), col="white", border=NA)
    
    lines(x = c(0, ws), y = c(0, 0)  )
    lines(x = c(0, 0), y = c(ws, 0)  )
    lines(x = c(ws, 0), y = c(ws, ws)  )
    lines(x = c(ws, ws), y = c(ws, 0)  )
  }
}
