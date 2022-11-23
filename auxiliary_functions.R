# intersection_three_circles.R
# Auxiliary functions to circles_area_of_overlap.R


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

# Function to generate the coordinates of N arragned points in a torus world of size ws
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
# plants is to be shown,.
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
      text(x, y, labels, cex=cex*0.98, font=1) # con  el problema del que no se veia bien
      # aca cambia el 0.98
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
                               col = col_circles)
        }else{
          points(x = com$x[j]+k[1], y = com$y[j]+k[2], col = com$ft, pch=19)
        }
        
      }}
    else{
      if(circle){
        plotrix::draw.circle(x = com$x[j], y = com$y[j], radius = com$sz[j], 
                             col = col_circles)
      }else{
        points(com$x, com$y, col=com$ft, pch=19)
      }
    }
    offst <- ws*0.97
    letter_sz <- 0.9
    
    
    if (!is.null(numbers)){
      
      if (sum(at_border) == 0){
        if (circle){
          text_box(x = com$x[j], y = com$y[j], labels = this_names, cex=letter_sz)
        }else{
          text_box(x = com$x[j]+(ws-offst), y = com$y[j]+(ws-offst), 
                   labels = this_names, cex=letter_sz)
        }
        
      }else if (sum(at_border)==2){
        if (circle){text_box(x = ws-offst, y = ws-offst, labels = this_names, cex=letter_sz)
          text_box(x = ws-offst, y = offst, labels = this_names, cex=letter_sz)
          text_box(x = offst, y = offst, labels = this_names, cex=letter_sz)
          text_box(x = offst, y = ws-offst, labels = this_names, cex=letter_sz)
        }else{
          text_box(x = ws-offst, y = ws-offst, labels = this_names, cex=letter_sz)
          text_box(x = ws-offst, y = offst, labels = this_names, cex=letter_sz)
          text_box(x = offst, y = offst, labels = this_names, cex=letter_sz)
          text_box(x = offst, y = ws-offst, labels = this_names, cex=letter_sz)
          
        }
        
      }else if (at_border[2]){
        if (circle){
          text_box(x = com$x[j], ws-offst, labels = this_names, cex=letter_sz)
          text_box(x = com$x[j], offst, labels = this_names, cex=letter_sz)
        }else{
          text_box(x = com$x[j]+ws-offst, ws-offst, labels = this_names, cex=letter_sz)
          text_box(x = com$x[j]+ws-offst, offst, labels = this_names, cex=letter_sz)
        }
        
        
      }else if (at_border[1]){
        
        if (circle){
          text_box(x = ws-offst , y = com$y[j], labels = this_names, cex=letter_sz)
          text_box(x = offst , y = com$y[j], labels = this_names, cex=letter_sz)
          
        }else{
          text_box(x = ws-offst , y = com$y[j]+ws-offst, labels = this_names, cex=letter_sz)
          text_box(x = offst , y = com$y[j]+ws-offst, labels = this_names, cex=letter_sz)
          
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


# # Next are several test to validate the functions
# # 
# # para probar esto, voy a usar el triangulo de reileux
# # dado que los tres centros estan a la distancia igual a sus radios (que
# # tambien son iguales) el area de interseccion debe ser
# # igual a 0.70477092301046  (https://rechneronline.de/pi/three-circles.php)
# #
# xes <- c(0, 1, 0.5)
# yes <-c(0, 0, sqrt(1-0.5**2))
# ra <- c(1, 1, 1)
# dist(matrix(c(xes,yes), ncol = 2))#triangulo reilauexx
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(are)))
# 
# 
# # prueba con triangulo relieaux grande
# # debe dar 281.908
# 
# xes <- c(0, 20, 20/2)
# yes <-c(0, 0, sqrt((20**2)-(10**2)))
# ra <- c(20, 20, 20)
# 
# dist(matrix(c(xes,yes), ncol = 2))#triangulo reilauexx
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(are)))
# 
# 
# 
# ## un segmento reflex
# 
# xes <- c(0, 0.5, 0.25)
# yes <-c(0, 0, 1)
# ra <- c(1, 1, 0.5)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(are)))
# 
# 
# #
# # # #uno no junto
# # #
# #
# #
# xes <- c(0, 2, 3)
# yes <-c(1, 0, -1)
# ra <- c(2, 1, 1)
# 
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#       bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# #
# #
# #
# # #
# # ### disjuntos
# #
# xes <- c(0, 5, 3)
# yes <-c(1, 0, -1)
# ra <- c(2, 1, 1)
# 
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# #
# #
# 
# uno dentro de los dos, pero no enteramente, caso dificil
# 
# xes <- c(0, 2, 1)
# yes <-c(0, 0, 0)
# ra <- c(2, 2, 1.2)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# #similar al anterios
# 
# xes <- c(0, 2, 0.5)
# yes <-c(0, 0, 0.1)
# ra <- c(2, 2, 1.6)
# 
# plot_circles_simple(xes, yes, ra)
# text(xes, yes, 1:length(xes))
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# intersection_two_circles(xes[1:2], yes[1:2], ra[1:2])
# 
# 
# #uno dentro de los dos, pero no enteramente, caso dificil, pero
# # debe ser un area muy cerca de la de c3
# 
# xes <- c(0, 2, 1)
# yes <-c(0, 0, 0)
# ra <- c(2, 2, 1.04)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#     bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
#
# 
# #uno dentro de los dos, pero no enteramente, caso  que parece dificil pero no
# #
# xes <- c(0, 2, 1)
# yes <-c(0, 0, 0)
# ra <- c(2, 2, 1.9)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# intersection_two_circles(xes[1:2], yes[1:2], ra[1:2])
# 
# # uno envolviendo a los  otros dos
# 
# xes <- c(1, 0, 2)
# yes <-c(0, 0, 0)
# ra <- c(4, 2, 2 )
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# intersection_two_circles(xes[2:3], yes[2:3], ra[2:3])
# 
# #
# #
# # # uno tocando a los otros dos, pero los otros dos no circles_intersecting caso 4
# # #
# xes <- c(0, 0.5, 3)
# yes <-c(0, 0, 0)
# ra <- c(4, 3, 2)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# intersection_two_circles(centers_x = xes[2:3],
#                          centers_y = yes[2:3],
#                          radii = ra[2:3])
# # #
# # #
# #
# # #
# # # # uno dentro de otro dentro de otro caso 6
# # #
# #
# xes <- c(0, 0, 0)
# yes <-c(0, 0, 0)
# ra <- c(4, 3, 2)
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 2**2*pi
# 
# #
# #
# #uno dentro de los otros dos caso 5
# 
# xes <- c(0, 2, 1)
# yes <-c(0, 0, 0)
# ra <- c(2, 2, 0.5)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# # # se tocan todos pero sin intersecccion
# 
# xes <- c(0, 0.5, 1)
# yes <-c(0, 1, 0)
# ra <- c(0.6, 0.6, 0.6)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# 
# ##pruebas aleatorias
# 
# 
# xes <- runif(3, -2, 2)
# yes <-runif(3, -2, 2)
# ra <- runif(3, 1, 4)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# 
# text(x = max(xes)+max(ra), y = min(yes)-max(ra),
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# fr <- c(1,2)
# intersection_two_circles(centers_x = xes[fr],
#                          centers_y = yes[fr],
#                          radii = ra[fr])
# 
# ra[3]**2*pi
# 
# 
# 
# 
# 
# #pruebas que fallaron antes
# 
# xes <- c(0.8204722, -1.3509740,  1.6051183)
# yes <- c(1.2302704, -0.3386334, -0.6597031)
# ra <- c(2.636195, 2.886777, 3.035400)
# #fallo porque no estaba un parentesis bien
# # y no  se hacia la multiplicacion cuando debia ser
# 
# plot_circles_simple(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# 
# plot_circles_simple(c(0,3,1), c(0,0,2), c(2,2,1))
# 
# intersection_three_circles(c(0,3,1), c(0,0,2), c(2,2,1))
# 
# # esta interseccion daba cero, pero no debe
# # tenia mal una distancia en la linea
# #if (! (x12 - (Dist[1,3]*costheta))**2 + (y12-(Dist[1,2]*sintheta))**2 < r[3]**2 ){return(0)}
# # ya se arreglo y se puso la segunda condicion que necesitaba
# 
# # este tenia en la seccion deel return referencia a center_x cuando ya debia ser cx
# 
# xes <- c(6, 1.5,2)
# yes <-c( 2,1.6,2)
# ra<- c( 4, 2, 1)
# plot_circles_simple(xes, yes, ra)
# intersection_two_circles(centers_x = xes[c(1,3)],
#                          centers_y = yes[c(1,3)],
#                          radii = ra[c(1,3)])
# 
# intersection_three_circles(xes, yes, ra)
# 
# #en la siguiente estaba mal que hacia referencia a radii en vez de r
# xes <- c(0, 3, 1.5, 2, 6)
# yes <-c(0, 0, 1.6, 2, 2)
# ra<- c(2, 2, 2, 1, 4)
# 
# intersection_two_circles(centers_x = xes[c(2,4)],
#                          centers_y = yes[c(2,4)],
#                          radii = ra[c(2,4)])
# 
# intersection_three_circles(xes[c(2,4,5)], yes[c(2,4,5)], ra[c(2,4,5)])
# 
# plot_circles_simple(xes[c(2,4,5)], yes[c(2,4,5)], ra[c(2,4,5)])
# 
# 
# 
# xes <- c(0.5, 0, 1)
# yes <-c( sqrt(1-0.5**2), 0, 0)
# ra <- c(2, 1, 1)
# 
# plot_circles_simple(xes, yes, ra)
# 
# intersection_two_circles(centers_x = xes[c(3,2)],
#                          centers_y = yes[c(3,2)],
#                          radii = ra[c(3,2)])
# intersection_three_circles(xes, yes, ra)
# 
# 
# # el siguiente fallaba porque hacia referencias a centers y en ves de a cy en
# # y2 <- cy[c1] + a * (cy[c2]-cy[c1])/D_c1_c2
# 
# xes <- c(1.3576603,  0.3199869, -0.2730404)
# yes<-c( 1.1513890, -0.1994831,  1.8560016)
# ra <- c(1.857299, 3.313472, 1.691115)
# 
# plot_circles_simple(xes, yes,ra)
# 
# intersection_two_circles(centers_x = xes[c(1,3)],
#                          centers_y = yes[c(1,3)],
#                          radii = ra[c(1,3)])
# 
# intersection_three_circles(xes, yes, ra)
# 
# 
# # esta tenia un error en intersection_two_circles, porque el arco era
# # muy grande, entonces estaba calculando un area de interseccion muy chica en
# # vez de muy grande, se cambio por completo la funcion
# 
# xes <- c(0,-0.6)
# yes <-  c(0,0)
# ra <- c(1, 0.5)
# plot_circles_simple(xes, yes,ra)
# points(xes[2], yes[2])
# abline(v=-0.925)
# 
# intersection_two_circles(xes, yes,ra)
# 
# (ra[2]**2*pi)-intersection_two_circles(xes, yes,ra)
# 
# 
# 
# xes <-c( 6.845009,  11.107935,10.000000 )
# yes<-c(  9.320638  ,  9.976079   , 10.000000  )
# ra <- c(3.329812, 2.047095, 1.346420)
# 
# plot_circles_simple(xes, yes,ra)
# 
# intersection_two_circles(centers_x = xes[c(2,1)],
#                          centers_y = yes[c(2,1)],
#                          radii = ra[c(2,1)])
# 
# intersection_three_circles(xes, yes, ra)
# 
# # esta daba un error, decia que si habia interzeccion entre las tres,
# # porque estaba usando el metodo de calcular la interseccion entre
# # el de adentro y que lo contiene, pero debia ser el de adento
# # y el que no es el lo contiene
# xes <- c(0, 0.5, 3)#, 0)
# yes <-c(0, 0.8660254, -2)#,-1)
# ra <- c(1,  2, 2)#, 1)
# 
# plot_circles_simple(xes, yes, ra)
# 
# intersection_three_circles(xes, yes, ra)
# 
# 
# # esta fallaba porque regresaba un area, cuando obviamente la interseccion
# # entre los tres es 0
# xes <- c( 8.897652, 10.596518,  8.858422)
# yes <- c(10.72434, 11.10207, 10.74569)
# ra <- c(0.8083996, 0.7807484, 0.7322310  )
# 
# plot_circles_simple(xes, yes, ra)
# intersection_three_circles(xes,yes, ra)
# # ya quedo, el problema era que se deteactaba como
# # un caso en el que el circulo estaba adentro de otro o algo mal hecho
# # simplififque totalmente la parte que calcula intersecciones entre dos
# # si al menos uno esta dentro de otro
# 
# 
# # en la siguiente obtenia un area negativa, obetnia esto porque en
# #    circle_minimum_intersection_dist <- which.min(lapply(intermediate_points_and_distance, #
# #function(x) x[["Distance"]])) estaba angle en vez de Distance
# xes <- c(10.87894,10.00000,  10.09625)
# yes <- c( 12.02520,10.00000,10.81140)
# ra <-c( 1.736221,1.646274,1.040372)
# #plot_circles_simple(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# # esta  parecia que no funcionaba , pero ya
# xes <- c(13.46621, 10.00000, 10.87894, 12.37296 ,10.09625)[c(1,2,4)]
# yes <- c(10.20452, 10.00000, 12.02520, 11.61642, 10.81140)[c(1,2,4)]
# ra <- c(1.670403, 1.824236, 1.884513, 1.838703, 1.083562)[c(1,2,4)]
# plot_circles_simple(xes, yes,ra)
# intersection_three_circles(xes, yes, ra)
# 
# # en la que sigue fallaba que era un caso dificil, pero al momento
# # determinar si se debia usar el arco mayor o menor, se obtenian angulos
# # negativos en atan2, osea no estaban de 0 a 2*pi entonces fallaba
# xes <- c(10.00000, 10.87894, 10.09625)
# yes <- c(10.0000, 12.0252, 10.8114)
# ra <- c(1.824236, 1.884513 ,1.083562)
# 
# 
# plot_circles_simple(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# 
# # aqui puede haber un error en el documento inicial de donde saque el algoritmo de calcular areas de tres
# # resulta que este es un caso con eso de reflex, entonces, se estaba calculando el arco menor
# # pero, segui inicialmente el algorimto al pie de la letra, consultando
# # otras formas ya vi que habia que restar r[2]**2*pi-etc, pero en el original
# # solo ponia -etc, agregando
# xes <- c(10.00000,10.15376, 10.72018 )
# yes <- c(10.000000,    9.684573, 9.555746)
# ra <- c(1.408165,  1.335831, 0.954278)
# # debe ser 2.28
# 
# plot_circles_simple(xes, yes, ra)
# 
# intersection_three_circles(xes, yes, ra)
# 
# intersection_two_circles(xes[c(2,3)], yes[c(2,3)], ra[c(2,3)])
# 
# # la que sigue tenia mal que al momento de calcular la diferencia
# # entre los angulos de interseccion dificil, la diferencia era negativa, entonces
# # parecia que el arco que debia tomarse era el mayor, pero no
# xes <- c(12.392308, 13.208398 ,10.000000 , 6.293676  ,7.250503,  7.707795 , 7.578273)
# yes <- c(6.593339 , 8.814251 ,10.000000 ,10.687837 , 8.256571  ,6.510737,  7.741898)
# ra <- c(2.143131, 2.098075, 2.229723, 2.347259, 1.714867, 2.066434, 1.439183)
# 
# settt<-c(5,6,7)
# plot_circles_simple(xes[settt], yes[settt], ra[settt])
# intersection_three_circles(xes[settt], yes[settt], ra[settt])
#
# 
# xes <- c(11.6672521762811,11.5962204686382,10,11.0147284614604,11.7018404358604,6.6888417678325,10.6840695374855,7.03345296153769,7.14089138405157,12.3342674745378,4.96631325245612)
# yes <- c(12.3230618524857,14.0627312838369,10,13.1352026226611,11.9761433395843,12.4910456440841,12.1515686531268,12.882647812017,12.3657794991253,14.839448134711,8.82829888531704)
# ra <- c(1.4689341610822,1.72477976773981,2.75153435963293,1.94509470421372,1.78272950009198,1.44599726243327,1.84572515365831,1.84519828236943,1.90716745542423,2.65134362052366,2.82634790078531)
# 
# settt<-c(3,6,8,9)
# 
# plot_circles_simple(xes[c(6,8,9)], yes[c(6,8,9)], ra[c(6,8,9)])
# intersection_three_circles(xes[c(6,8,9)], yes[c(6,8,9)], ra[c(6,8,9)])
# 
# xes <-c(0.9480961 ,1.7342106 ,0.8072776 ,1.6636906)
# yes <-c(1.2825796, 0.9270165, 1.4500958, 0.3893828)
# ra <- c(1.301322 ,1.414781, 1.496260 ,1.138502)
# 
# sett<- c(1,2,3)
# 
# plot_circles_simple(xes[sett], yes[sett], ra[sett])
# intersection_two_circles(xes[sett[c(1,2)]], yes[sett[c(1,2)]], ra[sett][c(1,2)])
# intersection_three_circles(xes[sett], yes[sett], ra[sett])
# 
# xes <- c(0.07244497, 1.09237271 ,0.91348859 ,1.24595679)[]
# yes <- c(0.1336482 ,1.5746221 ,1.4636216, 1.7886682)
# ra <-c(1.493899, 1.171559, 1.045887 ,1.424308)
# plot_circles_simple(xes[c(4,2,3)], yes[c(4,2,3)], ra[c(4,2,3)])
# 
# are<-intersection_three_circles(xes[c(4,2,3)], yes[c(4,2,3)], ra[c(4,2,3)])
# 
# text(x = 3, y = 0.5,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# xes <- c(0.57895659, 1.74898852, 1.59984376)
# yes <- c( 0.9427072, 1.5913909, 1.4699771)
# ra <- c(1.334639, 1.279570 ,1.116100)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# intersection_two_circles(xes[c(1,3)],yes[c(1,3)],ra[c(1,3)])
# text(x = 3, y = 0.5,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# xes <- c(1.7312687, 0.2162469, 1.2598769)
# yes <- c(0.9435113, 0.7290017 ,0.8874463)
# ra <- c(1.398440 ,1.132017, 1.048008)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# intersection_two_circles(xes[c(2,3)],yes[c(2,3)],ra[c(2,3)])
# 
# text(x = 3, y = 0,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# xes <- c(0.60173904, 0.54918222 , 0.63180893)
# yes <- c(0.92748815 ,0.01710736 , 0.30829904)
# ra <- c(1.361168 ,1.454600 , 1.256604)
# 
# plot_circles_simple(xes, yes, ra)
# 
# are<-intersection_three_circles(xes, yes, ra)
# intersection_two_circles(xes[c(2,3)],yes[c(2,3)],ra[c(2,3)])
# 
# text(x = 3, y = 0,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# 
# xes <- c(0.03318468, 0.57895659, 1.74898852 ,1.59984376)
# yes <-c(0.624486434273422,0.942707156762481,1.59139089891687,1.46997711015865)
# ra <- c(1.17046336608473,1.33463932166342,1.27956974727567,1.11610015179031)
# plot_circles_simple(xes[2:4], yes[2:4], ra[2:4])
# 
# are<-intersection_three_circles(xes[2:4], yes[2:4], ra[2:4])
# intersection_two_circles(xes[c(2,4)], yes[c(2,4)], ra[c(2,4)])
# text(x = 3, y = 0,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
########################################
###### End of tests



