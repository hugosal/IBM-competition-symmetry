
# centers:x, centers_y, radii, are two element vectors, and distance is one element
two_circles_inters_points <- function(centers_x, centers_y, radii, 
                                      circle_numbers = NULL,
                                      distance){
  
  a <- (radii[1]**2 -radii[2]**2+distance**2)/(2*distance)
  h<- sqrt(radii[1]**2 - a**2)
  x2 <- centers_x[1] + a * (centers_x[2]-centers_x[1])/distance
  y2 <- centers_y[1] + a * (centers_y[2]-centers_y[1])/distance
  
  dx <- h*(centers_y[2]-centers_y[1])/distance
  dy <-h*(centers_x[2]-centers_x[1])/distance
  
  x3_1 <- x2+dx
  x3_2 <- x2-dx
  
  y3_1 <- y2-dy
  y3_2 <- y2+dy
  if (! is.null(circle_numbers)){ 
    names1 <- paste(c(circle_numbers, "A"), collapse =  ":" )
    names2 <- paste(c(circle_numbers, "B"), collapse =  ":" )
  }else{
    names1 <- NULL
    names2 <- NULL}
  matrix(c(x3_1, x3_2, y3_1, y3_2), nrow = 2, 
         dimnames = list(c(names1, names2), 
                         c(NULL)))
}


intersection_two_circles <- function(centers_x, centers_y, radii){
  r1 <- radii[1]
  r2 <- radii[2]
  D <- sqrt( (centers_x[1] - centers_x[2])**2 + (centers_y[1] - centers_y[2])**2 )
  
  if (D <= abs(r1 - r2)){ # one circle contains the other
    return (pi * min(r1, r2)**2)
    
  }else if (D > r1 + r2){ # no intersection
    return (0) 
  }
  
  aa <- (r1**2 - r2**2 + D**2)/(2 * D)
  bb <- (r2**2 - r1**2 + D**2)/(2 * D)
  th1 <- 2 * acos(aa / r1)                          
  th2 <- 2 * acos(bb / r2)

  if (th1 > pi){
    a1 <- (r2**2 * (th2 + sin((2 * pi) - th2)))/2
    } else {
      a1 <- (r2**2 * (th2 - sin(th2)))/2}
  
  if (th2 > pi){
    a2 <- (r1**2 * (th1 + sin((2 * pi) - th1)))/2
    } else {
      a2 <- (r1**2 * (th1 - sin(th1)))/2
    }
  a1 + a2
  }

graficame_esta <- function(centers_x, centers_y, radii, npoints=500){
  
  plot(1, xlab="x", ylab="y", type="n", 
       xlim=c(min(centers_x)-max(radii), max(centers_x)+max(radii)*2),  
       ylim=c(min(centers_y)-max(radii), max(centers_y)+max(radii)*2))
  grid()
  colores <- col2rgb(rainbow(length(centers_y)), alpha = 0.2)

  for (i in 1:length(centers_y)){
    xes <- seq(from=centers_x[i]-radii[i], to = centers_x[i]+radii[i], 
               length.out=npoints )
    
    yes1 <- centers_y[i]+ sqrt(radii[i]**2-(xes-centers_x[i])**2)
    yes2 <- centers_y[i]- sqrt(radii[i]**2-(xes-centers_x[i])**2)
    
    xes <- c(xes[!is.nan(yes1)], rev(xes[!is.nan(yes2)]))
    yes <- c(yes1[!is.nan(yes1)], rev(yes2[!is.nan(yes2)]))
    polygon( xes, yes,  
             col=  rgb(colores[1, i], colores[2, i], 
                       colores[3, i], 255/2, maxColorValue=255))
    text(centers_x[i], centers_y[i] + radii[i], i)
    text(x =  max(centers_x)+max(radii ) , y= max(centers_y)+max(radii)*2-i+1, 
          bquote("C"~ .(i)~  A == .(round(radii[i]**2*pi ,3)) ))
    
     }
}


# This function computes the intersection area of three circles.
# several circle configurations have different overlap regions

intersection_three_circles <- function(centers_x, centers_y, radii){
  
  # circles must be sorted according to their radius
  order_rad <- order(radii, decreasing = TRUE)
  cx <- centers_x[order_rad]
  cy <- centers_y[order_rad]
  r <- radii[order_rad]
  Dist <- as.matrix(dist(matrix(c(cx, cy), ncol = 2)))
  
  if (any(combn(1:3, 2, 
                FUN = function(x){Dist[x[1], x[2] ] > r[x[1]] + r[x[2]] } ))){return(0)}
  
  # count which circles are inside of another circle.
  # Circles are ordered according to their radii, and since
  # only a larger cicle may contain a shorter one, so only 
  # combinations cmbdn(1:3) are needed
  
  circles_inside <-  combn(1:3, 2, 
                     FUN = function(x){
                       r[x[1]] >=  Dist[x[1], x[2]] + r[x[2]]} )

  # count which circles intersect
  circles_intersecting <- combn(1:3, 2, 
                 FUN = function(x){
                   if ((r[x[1]] - r[x[2]]) < Dist[x[1], x[2]] &  
                       Dist[x[1], x[2]] < (r[x[1]] + r[x[2]])){
                     x
                     }else{
                       c(-1,-1)}
                   })
  
  counts_intersect <- circles_intersecting[, circles_intersecting[1,] != -1]

  # If any circle is contained inside another, the intersection of the 
  # three circles is the intersection of the circle inside and the other one.
  if (any(circles_inside)){
    
    container_name <- c(1, 1, 2)
    contained_name <- c(2, 3, 3)
    
    maximum_container <- min(container_name[circles_inside])
    circle_most_inside <- max(contained_name[circles_inside])
    other_circle <- (1:3)[-c(maximum_container, circle_most_inside)]
    return(intersection_two_circles(centers_x = cx[c(circle_most_inside, other_circle)], 
                                    centers_y = cy[c(circle_most_inside, other_circle)],
                                    radii =  r[c(circle_most_inside, other_circle)]))
    
    } 
  
  # compute points of overlap between  a pair of circles, and see if 
  # the third other circle (g) contain both of them

  intersec_points_contained <- c(FALSE, FALSE , FALSE) 
   
  for (g in 1:3){
    c1 <- (1:3)[-g][1]
    c2 <- (1:3)[-g][2]
    D_c1_c2 <- Dist[c1, c2]
    
    intsrs_points <- two_circles_inters_points(centers_x = cx[c(c1, c2)],
                              centers_y =  cy[c(c1, c2)],
                              radii = r[c(c1, c2)],
                              distance = D_c1_c2)

    dist_matrix <- as.matrix(dist(matrix(c(cx[g], intsrs_points[, 1],
                                           cy[g], intsrs_points[, 2]), 
                                         ncol = 2  )))
    if (all(dist_matrix[1, 2:3] < r[g])){
      intersec_points_contained[g] <- TRUE
      }
    }
  
  if (sum(intersec_points_contained) == 2){

    # if one circle contains the intersection of the other two, the 
    # area is the area of the circle that contains the intersections, minus
    # the overlarp of this circle and the other two.
    
    containing_circle <- which(!intersec_points_contained)
    containing_circle_area <- r[containing_circle]**2 * pi
    
    area_remove <- 0
    
    for (f in (1:3)[-containing_circle]){
      area_inters <- intersection_two_circles(cx[c(f, containing_circle )], 
                               cy[c(f, containing_circle )], 
                                r[c(f, containing_circle )])

      area_remove <- area_remove + (containing_circle_area - area_inters) }
    
    return(containing_circle_area - area_remove)
    
  }else if (sum(intersec_points_contained) == 1){
    
    circs <- which(!intersec_points_contained)
    
    return(intersection_two_circles(centers_x = cx[circs], 
                                    centers_y = cy[circs],
                                    radii =  r[circs]))
  }
  
  # Next, is the case where the intersection area is a circular triangle
  
  if( !(r[1] - r[2]) < Dist[1, 2] &  Dist[1, 2]  < (r[1] + r[2])){return(0)}
  
  rs_1 <- r[1]**2
  rs_2 <- r[2]**2
  rs_3 <- r[3]**2
  
  d1_2 <- Dist[1, 2]
  d2_3 <- Dist[2, 3]
  d1_3 <- Dist[1, 3]
  
  x12 <- (rs_1 - rs_2 + d1_2**2) / (2 * d1_2)
  y12 <- sqrt(((2 * d1_2**2) * (rs_1 + rs_2)) - (rs_1 - rs_2)**2 - d1_2**4 )/(2 * d1_2)
  
  costheta <- (d1_2**2 + d1_3**2 - d2_3**2)/(2 * d1_2* d1_3)
  sintheta <- sqrt(1 - costheta**2)

  costhetap <- -(d1_2**2 + d2_3**2 - d1_3**2)/(2 * d1_2 * d2_3)
  sinthetap <- sqrt(1 - costhetap**2)
  
  if (! (x12 - (d1_3 * costheta))**2 + (y12 - (d1_3 * sintheta))**2 < rs_3 ){return(0)}
  
  if (! (x12 - (d1_3 * costheta))**2 + (y12 + (d1_3 * sintheta))**2 > rs_3 ){return(0)}
  
  x13p <- (rs_1 - rs_3 + d1_3**2) / (2 * d1_3)
  y13p <- -sqrt((2 * d1_3**2 * (rs_1 + rs_3)) - (rs_1 - rs_3)**2 - d1_3**4)/(2 * d1_3)
  
  x13 <- (x13p * costheta) - (y13p * sintheta)
  y13 <- (x13p * sintheta) + (y13p * costheta)

  x23pp <- (rs_2 - rs_3 + d2_3**2)/(2 * d2_3)
  y23pp <- sqrt((2 * d2_3**2 * (rs_2 + rs_3)) - (rs_2 - rs_3)**2 - d2_3**4)/(2 * d2_3)
  
  x23 <- (x23pp * costhetap) - (y23pp * sinthetap) + d1_2
  y23 <- (x23pp * sinthetap) + (y23pp * costhetap)

  c1 <- sqrt((x13 - x12)**2 + (y13 - y12)**2) # i3 j2 k1 
  c2 <- sqrt((x12 - x23)**2 + (y12 - y23)**2) # i1 j3 k2
  c3 <- sqrt((x23 - x13)**2 + (y23 - y13)**2) # i2 j1 k3

  A <- (sqrt((c1 + c2 + c3) * (c2 + c3 - c1) * (c1 + c3 - c2) * (c1 + c2 - c3))/4) +
    ((rs_1 * asin(c1 / (2 * r[1]))) - ( (c1/ 4) * sqrt((4 * rs_1) - c1**2))) +
    ((rs_2 * asin(c2 / (2 * r[2]))) - ( (c2/ 4) * sqrt((4 * rs_2) - c2**2)))

  if ((d1_3 * sintheta) < (y13 + ((y23 - y13) / (x23 - x13)) * ((d1_3 * costheta) - x13) )){
      arc_3 <- (r[3]**2*pi) - ((rs_3 * asin(c3 / (2 * r[3])))  - ( ((c3/ 4) * sqrt((4 * rs_3) - c3**2) )))
      }else{
    arc_3 <-((rs_3 * asin(c3 / (2 * r[3])))  - ( ((c3/ 4) * sqrt((4 * rs_3) - c3**2) ) ))
    }
  A + arc_3 
  }

# ### QUITAR COMENTARIOS AQUI
# 
# #pruebas
# 
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# 
# graficame_esta(c(0,3,1), c(0,0,2), c(2,2,1))
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes[c(2,4,5)], yes[c(2,4,5)], ra[c(2,4,5)])
# 
# 
# 
# xes <- c(0.5, 0, 1)
# yes <-c( sqrt(1-0.5**2), 0, 0)
# ra <- c(2, 1, 1)
# 
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes,ra)
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
# graficame_esta(xes, yes,ra)
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
# graficame_esta(xes, yes,ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# #graficame_esta(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# # esta  parecia que no funcionaba , pero ya
# xes <- c(13.46621, 10.00000, 10.87894, 12.37296 ,10.09625)[c(1,2,4)]
# yes <- c(10.20452, 10.00000, 12.02520, 11.61642, 10.81140)[c(1,2,4)]
# ra <- c(1.670403, 1.824236, 1.884513, 1.838703, 1.083562)[c(1,2,4)]
# graficame_esta(xes, yes,ra)
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
# graficame_esta(xes, yes, ra)
# intersection_three_circles(xes, yes, ra)
# 
# 
# # aqui puede haber un error en el documento inicial https://www.youtube.com/watch?v=Zr4jwxXP04w&list=RDGMEMQ1dJ7wXfLlqCjwV0xfSNbA
# # resulta que este es un caso con eso de reflex, entonces, se estaba calculando el arco menor
# # pero, segui inicialmente el algorimto al pie de la letra, consultando
# # otras formas ya vi que habia que restar r[2]**2*pi-etc, pero en el original
# # solo ponia -etc, agregando
# xes <- c(10.00000,10.15376, 10.72018 )
# yes <- c(10.000000,    9.684573, 9.555746)
# ra <- c(1.408165,  1.335831, 0.954278)
# # debe ser 2.28
# 
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes[settt], yes[settt], ra[settt])
# intersection_three_circles(xes[settt], yes[settt], ra[settt])
#
# 
# xes <- c(11.6672521762811,11.5962204686382,10,11.0147284614604,11.7018404358604,6.6888417678325,10.6840695374855,7.03345296153769,7.14089138405157,12.3342674745378,4.96631325245612)
# yes <- c(12.3230618524857,14.0627312838369,10,13.1352026226611,11.9761433395843,12.4910456440841,12.1515686531268,12.882647812017,12.3657794991253,14.839448134711,8.82829888531704)
# ra <- c(1.4689341610822,1.72477976773981,2.75153435963293,1.94509470421372,1.78272950009198,1.44599726243327,1.84572515365831,1.84519828236943,1.90716745542423,2.65134362052366,2.82634790078531)
# 
# settt<-c(3,6,8,9)
# 
# graficame_esta(xes[c(6,8,9)], yes[c(6,8,9)], ra[c(6,8,9)])
# intersection_three_circles(xes[c(6,8,9)], yes[c(6,8,9)], ra[c(6,8,9)])
# 
# xes <-c(0.9480961 ,1.7342106 ,0.8072776 ,1.6636906)
# yes <-c(1.2825796, 0.9270165, 1.4500958, 0.3893828)
# ra <- c(1.301322 ,1.414781, 1.496260 ,1.138502)
# 
# sett<- c(1,2,3)
# 
# graficame_esta(xes[sett], yes[sett], ra[sett])
# intersection_two_circles(xes[sett[c(1,2)]], yes[sett[c(1,2)]], ra[sett][c(1,2)])
# intersection_three_circles(xes[sett], yes[sett], ra[sett])
# 
# xes <- c(0.07244497, 1.09237271 ,0.91348859 ,1.24595679)[]
# yes <- c(0.1336482 ,1.5746221 ,1.4636216, 1.7886682)
# ra <-c(1.493899, 1.171559, 1.045887 ,1.424308)
# graficame_esta(xes[c(4,2,3)], yes[c(4,2,3)], ra[c(4,2,3)])
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes, yes, ra)
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
# graficame_esta(xes[2:4], yes[2:4], ra[2:4])
# 
# are<-intersection_three_circles(xes[2:4], yes[2:4], ra[2:4])
# intersection_two_circles(xes[c(2,4)], yes[c(2,4)], ra[c(2,4)])
# text(x = 3, y = 0,
#      bquote(A[1]*intersect(A[2])*intersect(A[3])==.(round(are,3))))
# 
# # #  # # ### QUITAR COMENTARIOS AQUI



