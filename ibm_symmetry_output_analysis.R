#### Analysis of simulated populations over time. 


# install.packages("Bolstad")
# install.packages("extraDistr")

library(Bolstad)
library(extraDistr)

# This file contains the results of a simulation 
# performed with "ibm_symmetry_experiment_main.R"

results_over_time <- read.csv("IBM_res_world_reachability_max_initial_sz_n_overcrowding_plants_comp_symmetry_20000_reps_20ws_50tsteps_26_seed.csv",
                              stringsAsFactors = F)


# Time 0 is the population at initialization, before any competition occurs, so it
# is omitted from this analysis
results_over_time <- results_over_time[results_over_time$t > 0, ]

set.seed(25)
time_limits <- c(1, 10, 20, 30, 40, 50)
accuracy_meassures <- data.frame(acc_BF_sym_asym = numeric(length(time_limits)),
                           mean_BF_sym = numeric(length(time_limits)),
                           miss_clas_sym = numeric(length(time_limits)),
                           meanBF_asym = numeric(length(time_limits)),
                           miss_clas_asym = numeric(length(time_limits))
)


rownames(accuracy_meassures) <- time_limits

results_over_time$comp_symmetry_factor <- factor(
  sapply(results_over_time$comp_symmetry, 
         function(x) if (x>1) "Assymetric" else "Symmetric"))

counter_days <- 1


train_index <- sample(unique(results_over_time$re), 
                      size = length(unique(results_over_time$re)) * 2/3, 
                      replace = F)

### revisa si esto hace falta 
  
# use 2/3 of simulation as training set, and 1/3 of simulations as testing set
trainin_repetitions <- results_over_time[results_over_time$re %in% train_index, ]
testing_repetitions <- results_over_time[!results_over_time$re %in% train_index, ]

training_repetitions_at_end <- trainin_repetitions[trainin_repetitions$t == 50, ]



reg_asym_mult <- bayes.lm(formula = coef_var ~ mean_no_competitors * sd_no_competitors,
                          data = training_repetitions_at_end[training_repetitions_at_end$comp_symmetry_factor=="Assymetric", ], 
                          center = F, x = T, y = T, sigma = F, model = T)

reg_sym_mult <- bayes.lm(formula = coef_var ~ mean_no_competitors*sd_no_competitors,
                         data = training_repetitions_at_end[training_repetitions_at_end$comp_symmetry_factor=="Symmetric", ], 
                         center = F,
                         x = T, y = T, sigma = F, model = T)

table_BFlm <- data.frame(cbind(
  c(reg_sym_mult$coefficients),
  c(summary(reg_sym_mult)$std.err),
  c(reg_asym_mult$coefficients),
  c(summary(reg_asym_mult)$std.err)
))

colnames(table_BFlm) <- c("Symetric Mean", "Symetric Std. Error",
                          "Asymetric Mean", "Asymetric Std. Error")
row.names(table_BFlm) <- c("Intercept", "Mean N", "SD N", "Mean N:SD N" )

# print(xtable::xtable(table_BFlm,digits=-2, display= c("s","g", "g", "g", "g")), 
#       math.style.exponents = TRUE)

summary(reg_sym_mult)

vbetainv <- reg_asym_mult$post.var[2:4, 2:4] # needs to be full rank, check with the determinant != 0

bbeta <- reg_asym_mult$coefficients[2:4]

(t(bbeta) %*% vbetainv) %*% bbeta

# This value must be lower than

qchisq(p = 0.05, lower.tail = F, df = ncol(vbetainv))

color_symmetric <- "#00BFC466"
color_asymmetric <- "#F8766D66" 

for (t_lim in time_limits){
  print(t_lim)
  
  testing_repetitions_at_time_t <- testing_repetitions[testing_repetitions$t == t_lim, ]
  
  # To get likelihood of the data under each model
  
  # First for asymmetric
  x_mean <- colMeans(reg_asym_mult$x[, c(2, 3) ]) 
  y_mean <- mean(reg_asym_mult$y)
  x_new <- as.matrix(cbind(rep(1, nrow(testing_repetitions_at_time_t)), 
                           testing_repetitions_at_time_t[, c("mean_no_competitors",
                                                             "sd_no_competitors") ]))
  df_asym <- reg_asym_mult$df.residual
  mean_asym <- apply(x_new, MARGIN = 1, FUN = function(xxx){
    t(as.matrix(reg_asym_mult$coefficients)) %*% c(xxx, xxx[2]*xxx[3] ) }) 
  XX <- as.matrix(reg_asym_mult$x)
  sigmaasym2 <- 1/df_asym *( t( reg_asym_mult$y -  XX %*% reg_asym_mult$coefficients) %*%
                               (reg_asym_mult$y - XX %*% reg_asym_mult$coefficients) )
  
  sigmaasym2 <- c(sigmaasym2) + apply(x_new, MARGIN = 1,
                                      function(xxx){
                                        t(c(xxx, xxx[2]*xxx[3] )) %*% solve(t(XX) %*% XX) %*% c(xxx, xxx[2]*xxx[3] )
                                      } )
  
  likelihood_asymm <- extraDistr::dlst(x = testing_repetitions_at_time_t$coef_var, 
                                       df = df_asym, 
                                       mu = mean_asym, sigma = sqrt(sigmaasym2) )
  
  # for symmetric
  x_mean_sym <- colMeans(reg_sym_mult$x[, c(2, 3) ]) 
  
  y_mean <- mean(reg_asym_mult$y)
  
  x_new <- as.matrix(cbind(rep(1, nrow(testing_repetitions_at_time_t)), 
                           testing_repetitions_at_time_t[, c("mean_no_competitors",
                                                             "sd_no_competitors") ] ))
  
  df_sym <- reg_sym_mult$df.residual
  
  mean_sym <- apply(x_new, MARGIN = 1, FUN = function(xxx){
    t(as.matrix(reg_sym_mult$coefficients)) %*% c(xxx, xxx[2]*xxx[3]) }) 
  
  XX <- as.matrix(reg_sym_mult$x)
  
  sigmasym2 <- 1/df_asym *( t( reg_sym_mult$y -  XX %*% reg_sym_mult$coefficients) %*%
                              (reg_sym_mult$y - XX %*% reg_sym_mult$coefficients) )
  
  sigmasym2 <- c(sigmasym2) + apply(x_new, MARGIN = 1,
                                    function(xxx){
                                      t(c(xxx, xxx[2]*xxx[3] )) %*% solve(t(XX) %*% XX) %*% c(xxx, xxx[2]*xxx[3] )
                                    } )
  
  likelihood_symm <- extraDistr::dlst(x = testing_repetitions_at_time_t$coef_var, 
                                      df = df_sym, 
                                      mu = mean_sym, sigma = sqrt(sigmasym2) )
  
  # BF_sym_asym > 1 will indicate that it is more likely that the population experienced symmetryc competiton, and viceversa
  BF_sym_asym <- likelihood_symm/likelihood_asymm
  
  # Distribution of parameters of missclasified populations under symmetric and asymmetric compeititon
  png(paste(t_lim ,"_sym_missclasified.png", sep = ""), width = 190, height = 160,
      units = 'mm', res = 600)

  par(mfrow = c(2,2))

  fancy_names <- c(bquote("Spatial randomness ("~kappa~")"),
                   bquote("Initial size variation ("~S[0]~")"),
                   "Overcrowding",
                   bquote("Competition symmetry ("~theta~")"))

  ranges <- list(world_reachability = 0:10,
                 max_initial_sz = round(seq(from = 0.5, to = 1.25, length.out = 12),2),
                 n_overcrowding_plants = seq(from = -10, to = 16, by = 2),
                 comp_symmetry = round(seq(from = 0, to = 1, length.out = 10), 2) )
  count <- 1

  for (var_name in c("world_reachability", "max_initial_sz",
                  "n_overcrowding_plants","comp_symmetry")){
    
        par(mar = c(5, if (count %% 2 != 0) 4 else 0.5, 2, 0.1) + 0.1)
    
    # The simulations that where symmetric and where incorrectly classified as asymmetric
    table_counts_missclassified_symm <- table(cut(testing_repetitions_at_time_t[testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" & BF_sym_asym < 1, var_name], 
                              breaks = ranges[[var_name]]))
    
    
    to_histogram <- list(breaks = ranges[[var_name]],
                         counts = table_counts_missclassified_symm / sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric") * 100,
                         density = table_counts_missclassified_symm / sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric") * 100,  
                         xname = var_name
                         )

        
    class(to_histogram) <- "histogram"

    plot(to_histogram, 
                freq = FALSE, col = color_symmetric, 
                xaxt = "n", yaxt = "n",
                main ="", ylab = if (count %% 2 != 0)  "Misclassified simulations (%)" else "",
                ylim = c(0, 5),
                xlab = fancy_names[count])
    
    axis(side = 1, at = to_histogram$breaks)
    if (count %% 2 != 0) axis(side = 2, at = c(0, 2.5,  5))
    count <- count + 1

  }

  dev.off()
 
  # Plots for missclatied asymmetric simulations

  png(paste(t_lim ,"_asym_missclasified.png", sep = ""),width = 190, height = 160,
      units = 'mm', res = 600)
  par(mfrow = c(2,2))
  
  fancy_names <- c(bquote("Spatial randomness ("~kappa~")"),
                   bquote("Initial size variation ("~S[0]~")"),
                   "Overcrowding",
                   bquote("Competition symmetry ("~theta~")"))
  
  ranges <- list(world_reachability = 0:10,
                 max_initial_sz = round(seq(from = 0.5, to = 1.25, length.out = 12),2),
                 n_overcrowding_plants = seq(from = -10, to = 16, by = 2),
                 comp_symmetry = seq(from = 1, to = 100, length.out = 12) )
  count <- 1
  
  for (var_name in c("world_reachability", "max_initial_sz",
                     "n_overcrowding_plants","comp_symmetry")){
    par(mar = c(5, if (count %% 2 != 0) 4 else 0.5, 2, 0.1) + 0.1)
    
    # The simulations that where symmetric and where incorrectly classified as asymmetric
    table_counts_missclassified_asymm <- table(cut(testing_repetitions_at_time_t[testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" &
                                                                                  BF_sym_asym > 1, var_name], 
                                                  breaks = ranges[[var_name]]))
    
    
    to_histogram <- list(breaks = ranges[[var_name]],
                         counts = (table_counts_missclassified_asymm /sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric")) * 100,
                         density = (table_counts_missclassified_asymm /sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric")) * 100,  
                         xname = var_name
    )
    class(to_histogram) <- "histogram"
    
    plot(to_histogram, 
         freq = FALSE, col = color_asymmetric, 
         xaxt = "n", yaxt = "n",
         main ="", ylab = if (count %% 2 != 0)  "Misclassified simulations (%)" else "",
         ylim = c(0, 10),
         xlab = fancy_names[count])
    
    axis(side = 1, at = to_histogram$breaks)
    
    if (count %% 2 != 0) axis(side = 2, at = c(0, 5,  10))
    count <- count + 1
    
  }
  
  dev.off()
  
  missc_sym <-sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" &
                    BF_sym_asym < 1)/sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric")
  
  missc_asym <- sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" &
                    BF_sym_asym > 1)/sum(testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric")
  

  accuracy_meassures[counter_days, 1] <- sum( (testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" & BF_sym_asym > 1) | 
                                                (testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" & BF_sym_asym < 1) ) / length(BF_sym_asym)*100
  accuracy_meassures[counter_days, 2] <- mean(BF_sym_asym[testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" & BF_sym_asym > 1 ])
  accuracy_meassures[counter_days, 3] <-  missc_sym * 100
  accuracy_meassures[counter_days, 4] <- mean(BF_sym_asym[testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" & BF_sym_asym < 1])
  accuracy_meassures[counter_days, 5] <- missc_asym * 100
  
  counter_days <- counter_days + 1
  
  
  
}

# print(xtable::xtable(accuracy_meassures, digits = c(1, 2, 2, 2, 2, 2)))

postscript_file <- FALSE

if(postscript_file){grDevices::cairo_ps("competiton_per_neighbours.eps",
                                        height = 5, width = 10,
                                        bg = F,fallback_resolution = 600)
}else{
  png("competiton_per_neighbours.png",
      width = 24, height = 12, units = "cm", res = 600)}
par(mfrow=c(1, 2), mar = c(5, 4, 4, 0) )
plot(y = results_over_time$mean_competiton[results_over_time$t==50],
     x = results_over_time$mean_no_competitors[results_over_time$t==50],
     col = c(color_asymmetric, 
             color_symmetric)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     pch = c(16,17)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     cex=0.4,  xlab = "Mean Neighbours", ylab = "Mean competition coefficient",
     yaxt="n",  xaxt ='n',yaxs="i")


axis(side = 1, at = c(0, 2, 4, 6))
axis(side = 2, at =c(0, 0.5, 0.75, 1), las =1)

par(mar = c(5, 0, 4, 2) )

plot(y = results_over_time$mean_competiton[results_over_time$t==50],
     x = results_over_time$sd_no_competitors[results_over_time$t==50],
     col = c(color_asymmetric, 
             color_symmetric)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     pch = c(16,17)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     cex=0.4,  xlab = "SD Neighbours",
     xlim = c(0, 3),
     ylab="", yaxt="n",  xaxt ='n',yaxs="i")

legend("topright", legend = c("Asymmetric", "Symmetric"),
       col = c(color_asymmetric, color_symmetric), cex = 1, pch = c(16,17), 
       title = "Symmetry of competition", bg = NA, box.col = NA)

axis(side = 1, at = c(0,1, 2, 3, 4))


dev.off()



postscript_file <- FALSE
if(postscript_file){grDevices::cairo_ps("final_coefvar_per_mean_sd_comp_symmetry.eps",
                                        height = 5, width = 10,
                                        bg = F,fallback_resolution = 600)
}else{
  png("final_coefvar_per_mean_sd_comp_symmetry.png",
      width = 24, height = 12, units = "cm", res = 600)}
par(mfrow=c(1, 2), mar = c(5, 4, 4, 0) )
plot(y = results_over_time$coef_var[results_over_time$t==50],
     x = results_over_time$mean_no_competitors[results_over_time$t==50],
     col = c(color_asymmetric, color_symmetric)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     pch = c(16,17)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     cex = 0.4,  xlab = "Mean Neighbours", 
     ylab = "Size coefficient of variation",
     yaxt="n",  xaxt ='n')

axis(side = 1, at = c(0, 2, 4, 6))
axis(side = 2, at =c(0, 0.25, 0.5), las = 1)

par(mar = c(5, 0, 4, 2) )

plot(y = results_over_time$coef_var[results_over_time$t==50],
     x = results_over_time$sd_no_competitors[results_over_time$t==50],
     col = c(color_asymmetric, color_symmetric)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     pch = c(16,17)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
     cex=0.4,  xlab = "SD Neighbours",
     xlim = c(0,3),
     ylab="", yaxt="n",  xaxt ='n',yaxs="i")

legend("topright", legend = c("Asymmetric", "Symmetric"),
       col = c(color_asymmetric, color_symmetric), cex = 1, pch = c(16,17), 
       title = "Symmetry of competition", bg = NA, box.col = NA)

axis(side = 1, at = c(0,1, 2, 3, 4))

dev.off()



# library(plotly)
# 
# plot_ly(results_over_time[results_over_time$t==50,], size = 1, type = "scatter3d",
#         x = ~mean_no_competitors, y = ~sd_no_competitors, z = ~coef_var, 
#         color = ~comp_symmetry_factor, colors =c("red", "blue"),
#         symbols = c('circle','x'), marker = list(size = 1))
# 
# plot3D::scatter3D(z = results_over_time$coef_var[results_over_time$t==50],
#                   x = results_over_time$sd_no_competitors[results_over_time$t==50],
#                   y = results_over_time$mean_no_competitors[results_over_time$t==50],
#                   col = c("red", "blue")[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
#                   pch = c(16,17)[results_over_time$comp_symmetry_factor[results_over_time$t==50]],
#                   cex=0.4,
#                   phi = 0, theta = 0, type = "p")




###########3 borrar

# plot(y = results_over_time$sd_no_competitors/results_over_time$mean_no_competitors,
#      x= results_over_time$mean_no_competitors, 
#      col = results_over_time$t)
# 
# 
# ranges_sym <- unique(results_over_time$t)
# colores_tiempo <- colorRampPalette(colors = c("darkblue", "lightblue"))
# 
# plot(y = results_over_time$sd_no_competitors,
#      x= results_over_time$mean_no_competitors,
#      pch=20, cex=0.6,
#      #xlim=c(1,3), ylim=c(1,1.7)
#      
# )
# 
# 
# plot(y = results_over_time$sd_no_competitors,
#      x= results_over_time$mean_no_competitors,
#      col=colores_tiempo(length(unique(results_over_time$t)))[
#        findInterval(x = results_over_time$t, vec = ranges_sym )],
#      pch=20, cex=0.6,
#      #xlim=c(1,3), ylim=c(1,1.7)
#           )
# 
# plot(y = results_over_time$sd_no_competitors[results_over_time$t==0],
#      x= results_over_time$mean_no_competitors[results_over_time$t==0],
#      pch=20, cex=0.6,
#      #xlim=c(1,3), ylim=c(1,1.7)
# )
# 
# plot(y = results_over_time$sd_no_competitors,
#      x= results_over_time$mean_no_competitors,
#      col=colores_tiempo(length(unique(results_over_time$t)))[
#        findInterval(x = results_over_time$t, vec = ranges_sym )],
#      pch=20, cex=0.6,
#      #xlim=c(1,3), ylim=c(1,1.7)
#      
# )
# 
# 
# 
# for(j in c( 2187
#             , 10586, 17580
#             )){
# lines(y = results_over_time$sd_no_competitors[results_over_time$re==j],
#        x= results_over_time$mean_no_competitors[results_over_time$re==j],
#        col="red", pch = 19, cex =0.5)
# }
# 
# results_over_time$re[results_over_time$sd_no_competitors <0.01  & 
#                        results_over_time$mean_no_competitors < 4.01 & 
#                        results_over_time$mean_no_competitors > 3.9]
# 
# abline(v=c(1,2,3,4,5), col="red")
