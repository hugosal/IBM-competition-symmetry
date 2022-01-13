# Librino_N
# Hugo Salinas, 2021

# The Librino_N function computes the areas of exclusive overlap 
# of N circles with x and y center coordinates  and  a positive radius.
# This is an implementation of the Librino-Levorato-Zorzi1 
# algorithm (10.1002/wcm.2305)

# Use the functions that compute the area of overlap of 2 and 3 circles

#Rprof(tmp <- tempfile(), line.profiling=TRUE)

source("intersection_three_circles.R")




# centers_x, centers_y, radii are numeric vectors of the same size
Librino_N <- function(centers_x, centers_y, radii){
  
  # This function returns the decimal number of circles in a binary code 
  
  circle_number_from_binary <- function(binary){
    which(strsplit(binary, "")[[1]]==1)
  }
  
  # This function return the transition matrix from n to n+t
  
  transition_n_to_k <- function(n, k, transitions){
    name_trans <- paste(n + k - 1, ":", n + k, sep = "")
    a_bar <- transitions[[name_trans]]
    if (k > 1){
      for (i in 1:(k - 1)){
        name_trans <- paste(n + k - i - 1, ":", n + k - i, sep = "")
        a_bar <- a_bar %*% transitions[[name_trans]]}
    }
    a_bar/factorial(k)
  }
  
  # This function returns a subtrellis given the nelement terminal node 
  # in the nvert subset 
  
  sub_trellis <- function(original, nvert, nelement, start){
    final_node <- names(original[[nvert]][nelement])
    reduced_trelis <- original[1:(nvert - 1)]
    n_final_node <- strsplit(final_node, "")[[1]]
    for (back in (nvert - 1):start){
      included_nodes <- sapply(names(reduced_trelis[[back]]), function(x){
        sum(!(strsplit(x, "")[[1]] == n_final_node)) == (nvert - back)
      })
      reduced_trelis[[back]] <- reduced_trelis[[back]][included_nodes] 
    }
    reduced_trelis
  }
  
  # This function returns the transition matrix of a corresponding 
  # list of areas of  overlap a
  
  transition_from_a_i <- function(a, start){
    transition_matrices_a <- list()
    # 1 is 0:1
    # 2 is 1:2 ...
    for (r in start:(length(a)- 1 )){
      next_lab <- paste(r, ":", r + 1, sep = "")
      names_rows <- lapply(names(a[[r + 1]]), function(x) strsplit(x, split = "")[[1]])
      names_cols <- lapply(names(a[[r]]),  function(x) strsplit(x, split = "")[[1]])
      
      this_matrix <- outer(X = 1:length(names_rows), Y =  1:length(names_cols), 
                  FUN = Vectorize(function(x, y){
                    sum(!(names_rows[[x]] ==  names_cols[[y]])) == 1}))
      
      rownames(this_matrix) <- names(a[[r+1]])
      colnames(this_matrix) <- names(a[[r]])

      transition_matrices_a[[next_lab]] <- this_matrix
    }
    transition_matrices_a
  }
  
  N <- length(centers_x)
  
  areas <- pi * radii**2
  
  if (N==2){
    intersection <- intersection_two_circles(centers_x, centers_y, radii)
    return (list("1:2" = intersection, 
                 "1" = areas[1] - intersection, 
                 "2" = areas[2] - intersection ))
  }else{
    
    d <- as.matrix(dist(matrix(c(centers_x, centers_y), ncol = 2)))
    vector_dec_from_bin <- numeric(strtoi(paste(rep("1", N), collapse = ""), 
                                          base = 2))
    
    names(vector_dec_from_bin) <- sapply(X = 1:length(vector_dec_from_bin),
                                         FUN = function(x){
                              bin <- paste(rev(as.integer(intToBits(x))), collapse="") 
                              substr(bin, start = nchar(bin) - N + 1 , stop = nchar(bin))
                                         } )
    
    a_i <- list()
    
    sum_binary <- sapply(names(vector_dec_from_bin), 
                         FUN = function(x){sum(as.numeric(strsplit(x, "")[[1]]))})
    
    for (j in 1:N){
      a_i[[j]] <- sort(which(sum_binary == j), decreasing = TRUE)
    }
    
    # # Transition matrices of the complete trellis 
    
    transition_matrices <- transition_from_a_i(a_i, start = 1)
    
    # Areas of a_1 are the areas of each circle
    
    a_i[[1]] <- sapply(names(a_i[[1]]), function(x){
      areas[circle_number_from_binary(x)]}) # nueva 
    
    # Areas of a_2 are the areas of each pairs of circles intersections
    
    a_i[[2]] <- sapply(names(a_i[[2]]), function(x){
      this_pair <- circle_number_from_binary(x)
      intersection_two_circles(centers_x = centers_x[this_pair],
                               centers_y = centers_y[this_pair],
                               radii = radii[this_pair])
    }) 
    
    # Areas of a_3 are the areas of each threecircle
    a_i[[3]] <- sapply(names(a_i[[3]]), function(x){
      this_three <- circle_number_from_binary(x)
      intersection_three_circles(centers_x = centers_x[this_three],
                                 centers_y = centers_y[this_three],
                                 radii = radii[this_three])
    })
    
    # Areas of intersection of a_n n>= 4 are computed using a_bar vectors
    # 
    if (N >= 4){
      for (u in 4:length(a_i)){
        
        for (v in seq_along(a_i[[u]])){
          Abar<-list()
          another_one <- u != 4 # this is becasue if  u=4, a_1 may be necessary
          minimum_depth <- if (u == 4) 1 else u - 2 # 1 - another_one
          subtrelis <- sub_trellis(nvert = u, nelement =  v, original = a_i, start = minimum_depth )
          transicion_sub <- transition_from_a_i(subtrelis, start = minimum_depth)

          for (e in minimum_depth:(u - 1 - another_one )){
            product_sum <- 0
            if (e < u - 1){

              for (j in (e + 1):(u - 1)){

                tnk <- transition_n_to_k(n = e,
                                         k = j - e,
                                         transitions = transicion_sub)
                this_rep <- (-1)**(j - e + 1) * (t(tnk) %*% as.matrix(subtrelis[[j]]))
                product_sum <- product_sum + this_rep
              }
            }
            Abar[[e]] <- subtrelis[[e]] - product_sum
          }

          if (u == 4){
            agam <- min(Abar[[1]])
            bgam <- max(-Abar[[2]])
            cgam <- min(Abar[[3]])

            if (any(subtrelis[[3]] <1e-6) ){
              a_i[[u]][v] <-  min(Abar[[u-1]])  # m = 0
            }else{
              if (cgam > bgam & # c > b
                  abs(cgam - cgam) < 1e-5){ # a = c

                four_circles <- circle_number_from_binary(names(a_i[[u]][v]))
                combinations <- combn(four_circles, m = 2)

                # get all circles intersection points
                intersc_pts <- lapply(1:ncol(combinations), FUN = function(x){
                  c1 <- combinations[1, x]
                  c2 <- combinations[2, x]
                  D_c1_c2 <- d[c1, c2]

                  two_circles_inters_points(centers_x = centers_x[c(c1, c2)],
                                            centers_y = centers_y[c(c1, c2)],
                                            radii = radii[c(c1, c2)],
                                            distance = D_c1_c2,
                                            circle_numbers = c(c1, c2))
                })

                intersections_mat <- do.call(rbind, intersc_pts)

                # get which intersection points are inside each circle,
                # not considering the intersections with that circle

                inters_pts_inside_circle <- lapply(four_circles, function(x){
                  distances_to_inters <- as.matrix(dist(rbind(
                    matrix(c(centers_x[x], centers_y[x]), nrow = 1),
                    intersections_mat)))[, 1]

                  not_including <-sapply(rownames(intersections_mat), function(y){
                    ! as.character(x) %in% strsplit(y, split = ":")[[1]]
                  })
                  unname(which(distances_to_inters[-1] < radii[x] & not_including))
                })

                the_one <- logical(4)

                # Find if the intersection of the points
                # inside one of the circles with the other circles is
                # just one element

                for (j in seq_along(inters_pts_inside_circle)){
                  the_one[j] <- all(sapply(seq_along(inters_pts_inside_circle)[-j],
                                           function(x){
                                             length(intersect(inters_pts_inside_circle[[j]],
                                                              inters_pts_inside_circle[[x]]))==1
                                           }))
                }

                if (all(unlist(lapply(inters_pts_inside_circle,
                                      function(x) length(x) == 3)))){
                  if (sum(the_one)==1){
                    a_i[[u]][v] <-   min(Abar[[u-1]])
                  }else{
                    a_i[[u]][v] <- max(-Abar[[u-2]])
                  }
                }else{
                  a_i[[u]][v] <-  min(Abar[[u-1]])}
              }else{
                a_i[[u]][v] <-  max(-Abar[[u-2]])
              }}
          }else{
            a_i[[u]][v] <- max(-Abar[[u-2]])
            }}
      }
    }
    
    # Compute exclusive intersection areas using a_i vectors
    Intersections_final <- list()
    
    for (i in 1:N){
      This_intersection_area <- a_i[[i]]
      product_sum <- 0
      if (i < N){ 
        for (j in (i+1):N){
          
          tnk <- transition_n_to_k(n = i, 
                                   k = j - i, 
                                   transitions = transition_matrices)
          
          this_rep <- ((-1)**(j - i + 1)) * (t(tnk) %*% as.matrix(a_i[[j]]))
          
          product_sum <- product_sum + this_rep
          
        }
      }
      
      Intersections_final[[i]] <- as.matrix( This_intersection_area - product_sum )
    }
    
    # customize labels and remove zero elements
    
    for (l in length(Intersections_final):1 ){
      this_vect <- Intersections_final[[l]]
      rownames(this_vect)<- sapply(rownames(this_vect), function(x){
        paste(circle_number_from_binary(x), collapse = ":")})
      #this_vect <- this_vect[this_vect > 1e-5, 1]
      if (length(this_vect)==0){
        Intersections_final[[l]] <- NULL
      }else{
        Intersections_final[[l]] <- this_vect[, 1]
      }
    }
    #print(a_i)
    Intersections_final
  }
} # end function


# Function to validate the Librino function.
#the area of a circle must be the sum of its partitions (intersections) 

validate_Librino <- function(librino, radii){
  passed_tests <- logical(length(radii))
  
  for(n in seq_along(radii)){
    summ <- 0
    if (class(librino)=="list"){
      for (i in seq_along(librino)){
        
        for (j in 1:nrow(librino[[i]])){
          nombres <- as.numeric(strsplit(rownames(librino[[i]]), split = ":")[[j]])
          if (n %in% nombres){
            
            summ <- summ + abs(librino[[i]][j])
          }
        } 
      }
      differencia <- (radii[n]**2*pi) - summ
      passed_tests[n] <- abs(differencia) < 1e-8 # arbitrary small number
      #print(paste(n, " diff Area original y reconstruida: ", differencia))
    }else{
      
      for (i in seq_along(librino)){
        nombres <- as.numeric(strsplit(names(librino)[i], ":")[[1]])
        if (n %in% nombres){
          summ <- summ + abs(librino[i])
        }
        differencia <- (radii[n]**2*pi) - summ}
      passed_tests[n] <- abs(differencia) < 1e-5 # arbitrary small number
      #print(paste(n, " diff Area original y reconstruida: ", differencia))
    }
  }
  if (all(passed_tests)){return(TRUE)}else{return( which(!passed_tests))}
}


###############################################
########## Simple test, passsed
# 
#  xes <- c(0, 1, 0.5)
#  yes <-c(0, 0, sqrt(1-0.5**2))
#  ra <- c(1, 1, 1)
# 
#  graficame_esta(xes, yes, ra)
# 
#  intersections <- Librino_N(centers_x = xes, centers_y = yes, radii = ra)
# 
#  intersection_two_circles(centers_x = xes[c(1,2)],
#                           centers_y = yes[c(1,2)],
#                           radii = ra[c(1,2)])
#  intersection_three_circles(xes, yes, ra)
# 
#  validate_Librino(unlist(intersections), ra)
# 
#  #
#  xes <- c(0, 0, 1, 1)
#  yes <- c(0, 1, 0, 1)
#  ra <- c(1, 1, 1, 1)
# 
#  graficame_esta(xes, yes, ra)
# 
#  intersections <- Librino_N(centers_x = xes, centers_y = yes, radii = ra)
# 
#  validate_Librino(unlist(intersections), ra)
# 
# # simple, y empiezo a complicarlo
# xes <- c(0, 1, 0.5, 3, 0)
# yes <-c(0, 0, 0.8660254, -2,-1)
# ra <- c(1, 1.3, 2, 2, 1)
# 
# graficame_esta(xes, yes, ra)
# 
# intersections <- Librino_N(centers_x = xes, centers_y = yes, radii = ra)
# intersections
# validate_Librino(unlist(intersections), ra)
# 
# # 
# # intersection_two_circles(centers_x = xes[c(2,4)],
# #                          centers_y = yes[c(2,4)],
# #                          radii = ra[c(2,4)])
# # intersection_three_circles(xes[c(1,2, 4)], yes[c(1,2,4)], ra[c(1,2,4 )])
# 
# xes <-c( 6.845009,  11.107935,10.000000 )
# yes<-c(  9.320638  ,  9.976079   , 10.000000  )
# ra <- c(3.329812, 2.047095, 1.346420)
# 
# graficame_esta(xes, yes,ra)
# 
# intersection_two_circles(centers_x = xes[c(1,2)],
#                          centers_y = yes[c(1,2)],
#                          radii = ra[c(1,2)])
# 
# intersection_three_circles(xes, yes, ra)
# 
# intersections <- Librino_N(centers_x = xes,
#           centers_y = yes,
#           radii = ra)
# 
# validate_Librino(unlist(intersections), radii = ra)
# 
# set.seed(21)
# xes <-runif(n = 8, 0, 10)#[c(1,3,6,7)]
# yes<- runif(n = 8,0, 10)#[c(1, 3,6,7)]
# ra <- runif(n = 8,1, 10)#[c(1,3,6,7)]
# 
# subset <- c(3,4,6,7)
# graficame_esta(xes[subset], yes[subset],ra[subset])
# 
# intersections <- Librino_N(centers_x = xes[subset],
#                            centers_y = yes[subset],
#                            radii = ra[subset])
# 
# validate_Librino(unlist(intersections), radii = ra[subset])
# 
# Validation using an ¿easy? setup of N circles
# 
# 
# 
# 
# 
# ene <- 10
# set.seed(26)
# 
# angles <- runif(n = ene, min = 0, max = 6*pi)
# 
# xes <- sin(angles)
# yes <- cos(angles)
# ra <- rep(2, ene)
# 
# de <- Librino_N(centers_x = xes, centers_y = yes, radii = ra)
# 
# validate_Librino(unlist(de), radii = ra)
# Rprof()
# summaryRprof(tmp, lines = "show")
# 
# 
#
# el siguiente fallaba porque estaba mal la comprobacion de si las intersecciones estan dentro de un criclo
# en librino


# ##########3 un comment aqui
# xes <- c(13.46621 ,10.00000 ,10.87894, 12.37296, 10.09625)
# yes <- c(10.20452, 10.00000 ,12.02520, 11.61642, 10.81140)
# ra<- c(1.670403, 1.824236 ,1.884513, 1.838703, 1.083562)
# 
# combn(2:5, m = 2, FUN = function(x){intersection_two_circles(xes[x], yes[x], ra[x])})
# 
# graficame_esta(xes[c(2:5)], yes[2:5], ra[2:5])
# intersss<- Librino_N(xes[2:5], yes[2:5], ra[2:5])
# validate_Librino(unlist(intersss), radii = ra[2:5])
# intersss
# kijuhg
# 
# xes <-c(10.000000,  9.562275 , 9.685901 ,11.710477 , 7.373770 , 7.513670 ,11.236690)
# yes <-c(10.000000 ,12.206473 , 9.425814 , 9.633818 ,11.798342,  8.383458,  6.138875)
# ra<-c( 1.676014, 1.919516, 1.756851, 2.224165, 2.228055 ,1.827564 ,2.465481)
# settt <- c(1, 2, 3, 4, 5)
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# 
# xes <-c(10.000000,  9.144737, 10.295128 ,13.022884  ,9.325139 , 9.520639)
# yes <- c( 10.000000 , 9.454299, 12.432144 ,10.003839 , 8.625922 ,11.333171)
# ra<-c(1.7571354 ,0.9525943 ,2.3654009 ,2.3134757, 2.3137928, 1.6474938)
# 
# settt<-c(1,3,4,5,6)
# graficame_esta(xes[c(1,3,5,6)], yes[c(1,3,5,6)], ra[c(1,3,5,6)])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# 
# # esta y la que sigue eran casos donde estaba mal determinar si era el caso espcial o no
# xes <-c(10.000000,  8.123759 , 9.619918,  8.835121, 10.765696 ,10.114364, 14.032726)
# yes <- c(10.000000 , 7.828132, 11.120877,  9.449705,  6.657127 , 8.213457, 10.154439)
# ra <-   c(1.634065, 2.025496 ,2.339616 ,2.205299, 1.869221, 2.100727 ,2.533315)
# 
# settt<-c(1, 3, 6, 7)
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# xes <- c(10.00000, 10.81787, 11.78556, 10.37228,  8.81388, 12.49571)
# yes <- c(10.000000 ,10.029176, 11.095502,  7.668324, 11.183218, 10.154315)
# ra <- c(1.730755, 1.108511, 1.908007, 1.728276, 1.276515, 1.807555)
# 
# settt<-c(1,2, 3, 6)
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# 
# 
# xes <- c(10.000000 , 8.123759 , 9.619918 , 8.835121, 10.765696, 10.114364, 14.032726)
# yes <- c(10.000000,  7.828132, 11.120877 , 9.449705  ,6.657127 , 8.213457, 10.154439)
# ra <-c(1.695014, 2.108357, 2.428915, 2.267899, 1.959408, 2.174822, 2.643399)
# 
# 
# settt<-c(1,2, 3, 6)
# graficame_esta(xes[settt], yes[settt], ra[settt])
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# # # la que sigue otra vez fallaba porque debe ser caso especial, agregue mas condiciones
# 
# xes <-c( 8.868264, 13.423852, 10.000000 ,12.619694, 11.648914 ,13.253822, 10.831798, 13.068225)
# yes <- c(12.280374, 10.771036 ,10.000000 ,10.983637 ,10.033991  ,9.995148,  9.272678,5.961888)
# ra <- c(2.693991, 1.382712 ,2.366476, 2.267781, 1.797292, 2.420490 ,2.098670, 2.878467)
# 
# settt <-c(3,4,5,8)#1:8#
# 
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# 
# xes <- c(5.928890, 10.000000  ,8.123759 , 9.619918,  8.835121, 10.765696, 10.114364,14.032726)
# yes <- c( 12.364275 ,10.000000,  7.828132, 11.120877 , 9.449705 , 6.657127 , 8.213457,10.154439)
# ra <-c(2.927484 ,1.808755 ,2.253090, 2.572428, 2.374518, 2.118149, 2.298077, 2.794330)
# settt <-c(2,4,5,7)
# #
# graficame_esta(xes[settt], yes[settt], ra[settt])
# #
# des <-unlist(Librino_N(xes[settt], yes[settt], ra[settt]))
# validate_Librino(des, radii = ra[settt])
# graficame_esta(xes[c(2,4,7)], yes[c(2,4,7)], ra[c(2,4,7)])
# intersection_three_circles(xes[c(2,4,7)], yes[c(2,4,7)], ra[c(2,4,7)])
# 
# 
# xes <- c(11.876241 ,10.000000 ,11.496159 ,10.711362, 12.641938, 11.990605,  7.137971)
# yes <-c(12.171868 ,10.000000, 13.292745, 11.621573 , 8.828995, 10.385324  ,9.115377)
# ra <-c(1.570166, 1.934042 ,2.236651, 2.134641, 1.770881, 2.016495, 2.681927)
# settt <-c(1,4, 3,5, 6)
# 
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# xes <- c(10.00000, 10.81787, 11.78556, 10.37228,  8.81388, 12.49571)
# yes <- c(10.000000 ,10.029176, 11.095502,  7.668324, 11.183218, 10.154315)
# ra <- c(2.123066 ,1.290069 ,2.342230 ,2.243891, 1.901301, 2.273014)
# 
# settt<-c(2, 3,4, 5)
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# # xcon la que sigue habia un error en el validador, error de redondeo, quitar
# # que librino automaticamente quite intersecciones pqeuenias hasta despues de validar
# xes <- c(12.392308, 13.208398 ,10.000000 , 6.293676  ,7.250503,  7.707795 , 7.578273)
# yes <- c(6.593339 , 8.814251 ,10.000000 ,10.687837 , 8.256571  ,6.510737,  7.741898)
# ra <- c(2.143131, 2.098075, 2.229723, 2.347259, 1.714867, 2.066434, 1.439183)
# 
# settt<-c(3,5,6,7)
# 
# graficame_esta(xes[settt], yes[settt], ra[settt])
# 
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
# xes <- c(11.6672521762811,11.5962204686382,10,11.0147284614604,11.7018404358604,6.6888417678325,10.6840695374855,7.03345296153769,7.14089138405157,12.3342674745378,4.96631325245612)
# yes <- c(12.3230618524857,14.0627312838369,10,13.1352026226611,11.9761433395843,12.4910456440841,12.1515686531268,12.882647812017,12.3657794991253,14.839448134711,8.82829888531704)
# ra <- c(1.4689341610822,1.72477976773981,2.75153435963293,1.94509470421372,1.78272950009198,1.44599726243327,1.84572515365831,1.84519828236943,1.90716745542423,2.65134362052366,2.82634790078531)
# 
# settt<-c(3,6,8,9)
# 
# graficame_esta(xes[settt], yes[settt], ra[settt])
# intersection_two_circles(xes[c(6,8)], yes[c(6,8)], ra[c(6,8)])
# validate_Librino(unlist(Librino_N(xes[settt], yes[settt], ra[settt])), radii = ra[settt])
# 
#
# random tests
# set.seed(26)
# enes <-1000
# for (i in 1:enes){
   xes <- runif(6, min = 0, max = 2)
   yes <- runif(6, min = 0, max = 2)
   ra <- runif(6, min = 1, max = 1.5)                               #                                 radii = ra
   validate_Librino(unlist(Librino_N(xes, yes, ra)),  radii = ra)
                    #                                 radii = ra)
   # if (is.numeric(validate_Librino(unlist(Librino_N(xes, yes, ra)), 
#                                 radii = ra))){
#  
#   stop( paste("falla en ", i))
# }
#  if (i==enes) print("exito")}
# 
# #graficame_esta(xes, yes, ra)
#   
# library(rbenchmark)
# benchmark("2ndopt" = {
#   set.seed(26)
#   enes <-1
#   for (i in 1:enes){
#     xes <- runif(11, min = 0, max = 4)
#     yes <- runif(11, min = 0, max = 4)
#     ra <- runif(11, min = 1, max = 1.5)
#     if (is.numeric(validate_Librino(unlist(Librino_N(xes, yes, ra)),
#                                     radii = ra))){
# 
#       stop( paste("falla en ", i))
#     }
#     #if (i==enes) print("exito")
#   }
# },
# replications = 2,
# columns = c("test", "replications", "elapsed",
#             "relative", "user.self", "sys.self"))
# 
# 
