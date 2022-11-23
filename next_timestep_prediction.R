# Script to model and predict variation tren only for the next timestep
# of simluations
library(tree)
library(Bolstad)
results_over_time <- read.csv("IBM_res_world_reachability_max_initial_sz_n_overcrowding_plants_comp_symmetry_8400_reps_20ws_50tsteps_23_seed.csv",
                              stringsAsFactors = F)

set.seed(26)
time_limits <- c(1, 10, 20, 30, 40, 50)
res_acc_days <- data.frame(acc_BF = numeric(length(time_limits)),
                           mean_bf_sym = numeric(length(time_limits)),
                           meanbf_asym = numeric(length(time_limits)),
                           reg_tree_deltacoef_acc = numeric(length(time_limits)),
                           mean_deltacoef_sym  = numeric(length(time_limits)),
                           mean_deltacoef_asym = numeric(length(time_limits))           
)


rownames(res_acc_days) <- time_limits
results_over_time$comp_symmetry_factor <- factor(
  sapply(results_over_time$comp_symmetry, 
         function(x) if (x>1) "Assymetric" else "Symmetric"))



counter_days <- 1
for (t_lim in time_limits){
  
  train_index <- sample(unique(results_over_time$re), 
                        size = length(unique(results_over_time$re))*2/3, 
                        replace = F)
  
  # use 2/3 of simulation as training set, and 1/3 of simulations as testing set
  trainin_repetitions <- results_over_time[results_over_time$re %in% train_index, ]
  testing_repetitions <- results_over_time[!results_over_time$re %in% train_index, ]
  
  time_of_interest <- t_lim # max(trainin_repetitions$t)
  
  training_repetitions_at_time_t <- trainin_repetitions[
    trainin_repetitions$t == time_of_interest, ]
  
  testing_repetitions_at_time_t <- testing_repetitions[
    testing_repetitions$t == time_of_interest, ]
  
  ###### Compare three methods for prediction thew symmetry of competition
  ## 1) regression tree using only increase of decrease from the previous time to this one
  # for a certain time, compute the change in varaition with the previous step
  training_deltacoef <-trainin_repetitions[trainin_repetitions$t != max(trainin_repetitions$t), ]
  training_deltacoef$dcoef_Var <- trainin_repetitions$coef_var[trainin_repetitions$t != min(results_over_time$t)] - training_deltacoef$coef_var
  training_deltacoef <- training_deltacoef[training_deltacoef$t == time_of_interest - 1, ]
  
  testing_deltacoef <-testing_repetitions[testing_repetitions$t != max(testing_repetitions$t), ]
  testing_deltacoef$dcoef_Var <- testing_repetitions$coef_var[testing_repetitions$t != min(results_over_time$t)] - testing_deltacoef$coef_var
  testing_deltacoef <- testing_deltacoef[testing_deltacoef$t == time_of_interest - 1, ]
  
  # x_lim <- max(abs(training_deltacoef$dcoef_Var))
  # plot(training_deltacoef$dcoef_Var[training_deltacoef$comp_symmetry > 1], 
  #      rep(0.3, sum(training_deltacoef$comp_symmetry > 1)), pch='|', xlim = c(-x_lim, x_lim),
  #      main=time_of_interest, col="red", ylim = c(0,0.4), yaxt = "n", 
  #      ylab="", xlab="d Coef. Var.")
  # 
  # points(training_deltacoef$dcoef_Var[training_deltacoef$comp_symmetry < 1], 
  #        rep(0.1, sum(training_deltacoef$comp_symmetry < 1)), pch='|', col = "blue")
  # axis(side = 2, at = c(0.1, 0.3), labels = c("Symmetric", "Asymmetric"))
  # abline(v=0)
  
  delta_coef_tree <- tree(comp_symmetry_factor ~ dcoef_Var,
                          data = training_deltacoef)
  # Prediction of symmetry
  tree_accuracy_this_rep_dcoef <- sum(testing_deltacoef$comp_symmetry_factor == predict(object = delta_coef_tree, 
                                                                                        newdata =  testing_deltacoef,
                                                                                        type="class"))/nrow(testing_deltacoef)
  
  mean_dcoef_asym <- mean(training_deltacoef$dcoef_Var[training_deltacoef$comp_symmetry_factor=="Assymetric"])
  mean_dcoef_sym <- mean(training_deltacoef$dcoef_Var[training_deltacoef$comp_symmetry_factor=="Symmetric"])
  # multiple linear regressio bayesian model
  
  reg_asym_mult <- bayes.lm(formula = coef_var ~ mean_no_competitors * sd_no_competitors,
                            data=training_repetitions_at_time_t[training_repetitions_at_time_t$comp_symmetry_factor=="Assymetric", ], center = F,
                            x = T, y = T, sigma = F, model = T)
  
  
  x_mean <- colMeans(reg_asym_mult$x[, c(2, 3) ]) 
  
  # xes <- seq(from=0,to=7, length.out=20)
  # yes <- seq(from=0,to=3, length.out=20)
  # z <-outer(X = xes, Y = yes,FUN = Vectorize(function(x, y){
  #         t(as.matrix(reg_asym_mult$coefficients)) %*% matrix(c(1, x, y, x*y), nrow = 4)
  # 
  #       }  ))
  # saveGIF({
  # for(ki in seq(0, 360, length.out=30)){
  # pespective <- persp(xes, yes, z, theta = ki, phi = -2, col = "red", shade = 0.5, 
  #                     zlim = c(0, 0.5), xlab = "Mean. comp",
  #                     ylab = "SD. comp", zlab = "Final Coef. Var.")
  # 
  # points( trans3d(x = testing_repetitions_at_time_t$mean_no_competitors,
  #                 y = testing_repetitions_at_time_t$sd_no_competitors,
  #                 z = testing_repetitions_at_time_t$coef_var, 
  #                 pespective),
  #         pch=18,  
  #         #col=red2blue(n = 8)[color_BF], 
  #         col =  c("red", "blue") [testing_repetitions_at_time_t$comp_symmetry_factor], 
  #         cex=0.7)}
  # }, movie.name = "posteriormean_asym.gif", interval = 1, 
  # ani.width = 1500, ani.height = 1000, clean = TRUE,  ani.res = 150)
  
  x_mean <- colMeans(reg_asym_mult$x[, c(2, 3) ]) 
  
  y_mean <- mean(reg_asym_mult$y)
  
  x_new <- as.matrix(cbind(rep(1, nrow(testing_repetitions_at_time_t)), 
                           testing_repetitions_at_time_t[, c(8, 9) ]))
  
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
  
  ############### for symmetry model 
  reg_sym_mult <- bayes.lm(formula = coef_var ~ mean_no_competitors*sd_no_competitors,
                           data=training_repetitions_at_time_t[training_repetitions_at_time_t$comp_symmetry_factor=="Symmetric", ], 
                           center = F,
                           x = T, y = T, sigma = F, model = T)
  
  # zasym <-outer(X = xes, Y = yes, FUN = Vectorize(function(x, y){
  #   t(as.matrix(reg_sym_mult$coefficients)) %*% matrix(c(1, x, y, x*y), nrow = 4)
  #   } ))
  # saveGIF({
  # for(ki in seq(0, 360, length.out=30)){
  #   pespective <- persp(xes, yes, zasym, theta = ki, phi = -2, col = "blue", shade = 0.5, 
  #                       zlim = c(0, 0.5), xlab = "Mean. comp",
  #                       ylab = "SD. comp", zlab = "Final Coef. Var.")
  #   points( trans3d(x = testing_repetitions_at_time_t$mean_no_competitors,
  #                   y = testing_repetitions_at_time_t$sd_no_competitors,
  #                   z = testing_repetitions_at_time_t$coef_var, 
  #                   pespective),
  #           pch=18,  col =  c("red", "blue") [testing_repetitions_at_time_t$comp_symmetry_factor], 
  #           cex=0.7)}
  # }, movie.name = "posteriormean_sym.gif", interval = 1, 
  # ani.width = 1500, ani.height = 1000, clean = TRUE,  ani.res = 150)
  
  
  x_mean_sym <- colMeans(reg_sym_mult$x[, c(2, 3) ]) 
  
  y_mean <- mean(reg_asym_mult$y)
  
  x_new <- as.matrix(cbind(rep(1, nrow(testing_repetitions_at_time_t)), 
                           testing_repetitions_at_time_t[, c(8, 9) ] ))
  
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
  
  BF <- likelihood_symm/likelihood_asymm
  
  res_acc_days[counter_days, 1] <- sum( (testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" & BF > 1) | 
                                          (testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" & BF < 1) ) / length(BF)*100
  
  res_acc_days[counter_days, 2] <- mean(BF[testing_repetitions_at_time_t$comp_symmetry_factor == "Symmetric" & BF > 1 ])
  
  res_acc_days[counter_days, 3] <- mean(BF[testing_repetitions_at_time_t$comp_symmetry_factor == "Assymetric" & BF < 1])
  
  res_acc_days[counter_days, 4] <- tree_accuracy_this_rep_dcoef  * 100
  
  res_acc_days[counter_days, 5] <- mean_dcoef_sym
  
  res_acc_days[counter_days, 6] <- mean_dcoef_asym
  
  counter_days <- counter_days + 1
  
  
  postscript_file <- FALSE
  if(postscript_file){grDevices::cairo_ps("final_coefvar_per_mean_sd_comp_symmetry.eps",
                                          height = 5, width = 10,
                                          bg = F,fallback_resolution = 500)
  }else{
    png(paste(t_lim ,"final_coefvar_per_mean_sd_comp_symmetry.png"),width = 190, height = 160, 
        units = 'mm', res = 600)}
  only_last <-  results_over_time[results_over_time$t==t_lim, ]
  # tiff("final_coefvar_per_mean_sd_comp_symmetry.tiff",
  #     width = 16, height = 8, units = "cm", res = 500)
  par(mfrow=c(1, 2), mar = c(5, 4, 4, 0) )
  plot(x = only_last$mean_no_competitors, y = only_last$coef_var,
       col = c("red", "blue")[only_last$comp_symmetry_factor],
       pch=20, cex=0.4, xlab= "Mean Neighbours",
       ylab="Size coefficient of variation", xaxs="i", yaxs="i",
       xaxt = "n",  yaxt = "n", ylim=c(0,0.6))
  
  axis(side = 1, at = c(0,2,4,6))
  axis(side = 2, at = c(0,0.2,0.4), las=1)
  
  
  
  legend("topleft", legend = c("Asymmetric", "Symmetric"),
         fill = c("red","blue"), cex=1.2)
  par(mar = c(5, 0, 4, 2) )
  plot(x = only_last$sd_no_competitors, y = only_last$coef_var,
       col = c("red", "blue")[only_last$comp_symmetry_factor],
       pch=20, cex=0.4,  xlab = "SD Neighbours",
       ylab="", yaxt="n",  xaxt ='n',yaxs="i",
       ylim=c(0,0.6))
  
  axis(side = 1, at = c(0,2,4,6))
  
  
  dev.off()
  
  
  
}

print(xtable::xtable(res_acc_days, digits =c(1,2,2,2,2,-2,-2)))

colnames(testing_repetitions_at_time_t)
#hist(testing_repetitions_at_time_t$mean_no_competitors[false_asym])

# 
# ######
# 
# 
# model_final_coefvar <- lm(coef_var~mean_no_competitors*sd_no_competitors-1,
#                           data=only_last)
# summary(model_final_coefvar)
# 
# 
# png("fimnal_coefvar_per_meancomp_symmetry.png",
#     width = 24, height = 12, units = "cm",
#     res = 300)
# par(mfrow=c(1,2))
# plot(x = only_last$mean_no_competitors, y = only_last$coef_var,
#      #col= rgb(blue=only_last$sd_no_competitors/max(only_last$sd_no_competitors), green = 0,red = 0),
#      #col= rgb(blue=only_last$max_initial_sz-0.5/max(only_last$max_initial_sz), green = 0,red = 0),
#      col = c("red", "blue")[only_last$comp_symmetry_factor],
#      #col=topo.colors(25)[as.factor(only_last$n_overcrowding_plants)],
#      pch=20, cex=0.8, xlab="Mean Comp.",
#      ylab="Final coefficient of variation" )
# 
# legend("topleft", legend = c("Asymmetric", "Symmetric"),
#        fill = c("red","blue"), cex=0.8)
# 
# plot(x = only_last$sd_no_competitors, y = only_last$coef_var,
#      #col= rgb(blue=only_last$sd_no_competitors/max(only_last$sd_no_competitors), green = 0,red = 0),
#      #col= rgb(blue=only_last$max_initial_sz-0.5/max(only_last$max_initial_sz), green = 0,red = 0),
#      col = c("red", "blue")[only_last$comp_symmetry_factor],
#      #col=topo.colors(25)[as.factor(only_last$n_overcrowding_plants)],
#      pch=20, cex=0.8, xlab="SD Comp.",
#      ylab="Final coefficient of variation" )
# legend("topleft", legend = c("Asymmetric", "Symmetric"),
#        fill = c("red","blue"), cex=0.8)
# #abline(model_final_coefvar, lwd=2, lty=1, col="red")
# 
# dev.off()
# 
# only_last <- results_over_time[results_over_time$t == max(results_over_time$t),]
# 
# 
# reg_asym <- bayes.lin.reg(y = only_last$coef_var[only_last$comp_symmetry_factor=="Assymetric"],
#                           x = only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Assymetric"],
#                           plot.data = F, slope.prior = "flat", intcpt.prior = "flat", 
#                           pred.x=testing_repetitions_at_time_t$mean_no_competitors,
#                           drawPlot=F, plot = FALSE, quiet=T )
# 
# reg_symm <- bayes.lin.reg(y = only_last$coef_var[only_last$comp_symmetry_factor=="Symmetric"],
#                           x = only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Symmetric"],
#                           plot.data = F, slope.prior = "flat", intcpt.prior = "flat",
#                           pred.x=testing_repetitions_at_time_t$mean_no_competitors, drawPlot=F,
#                           plot = FALSE, quiet=T)
# 
# 
# 
# x_range <-seq(from=0, to=, 10, length.out=100)
# 
# interval_assym <- bayes.lin.reg(y = only_last$coef_var[only_last$comp_symmetry_factor=="Assymetric"],
#                                 x = only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Assymetric"],
#                                 plot.data = F, slope.prior = "flat", intcpt.prior = "flat", 
#                                 pred.x=x_range,  plot = FALSE, quiet=T)
# 
# mean_asym <-  interval_assym$pred.y
# #reg_asym$post.coef[1] +
# #(x_range-mean(only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Assymetric"]))*reg_asym$post.coef[2]
# 
# sd_asym <- interval_assym$pred.se
# #(reg_asym$post.coef.sd[1]**2) +
# #   (((x_range-mean(only_last$mean_no_competitors[
# #     only_last$comp_symmetry_factor=="Assymetric"]))**2) * (reg_asym$post.coef.sd[2]**2))+
# #   var(only_last$coef_var[only_last$comp_symmetry_factor=="Assymetric"])
# 
# upr95_asym <- qnorm(p = 0.975, mean = mean_asym, sd= sd_asym)
# low95_asym <- qnorm(p = 0.025, mean = mean_asym, sd= sd_asym)
# 
# 
# interval_sym <- bayes.lin.reg(y = only_last$coef_var[only_last$comp_symmetry_factor=="Symmetric"],
#                               x = only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Symmetric"],
#                               plot.data = F, slope.prior = "flat", intcpt.prior = "flat", 
#                               pred.x=x_range,  plot = FALSE, quiet=T)
# 
# mean_sym <-  interval_sym$pred.y
# #reg_asym$post.coef[1] +
# #(x_range-mean(only_last$mean_no_competitors[only_last$comp_symmetry_factor=="Assymetric"]))*reg_asym$post.coef[2]
# 
# sd_sym <- interval_sym$pred.se
# 
# #(reg_asym$post.coef.sd[1]**2) +
# #   (((x_range-mean(only_last$mean_no_competitors[
# #     only_last$comp_symmetry_factor=="Assymetric"]))**2) * (reg_asym$post.coef.sd[2]**2))+
# #   var(only_last$coef_var[only_last$comp_symmetry_factor=="Assymetric"])
# 
# upr95_sym <- qnorm(p = 0.975, mean = mean_sym, sd= sd_sym)
# low95_sym <- qnorm(p = 0.025, mean = mean_sym, sd= sd_sym)
# 
# 
# 
# png("BF_symm_prediction_coef_95.png",
#     width = 12, height = 12, units = "cm",
#     res = 300)
# 
# plot(1, type="n", xlim=c(0, 9), ylim=c(0, 0.5), xlab="Mean competitors", 
#      ylab="Coefficient of variation", xaxs="i" )
# lines(x_range, mean_asym, col="red", lwd=2, lty = 1)
# polygon(y = c(upr95_asym, rev(low95_asym)), x= c(x_range, rev(x_range)), 
#         col = rgb(1, 0,0, alpha = 0.08), border = "red", lty = 2)
# lines(x_range, mean_sym, col="blue", lwd=2, lty = 1)
# 
# polygon(y = c(upr95_sym, rev(low95_sym)), x= c(x_range, rev(x_range)), 
#         col = rgb(0, 0,1, alpha = 0.08), border = "blue", lty = 2)
# 
# points(x = testing_repetitions_at_time_t$mean_no_competitors, 
#        y = testing_repetitions_at_time_t$coef_var, pch=19,cex=0.6, 
#        col = c(rgb(as.data.frame(t(col2rgb( "red" ))), alpha=90, maxColorValue=255 ), 
#                rgb(as.data.frame(t(col2rgb( "blue" ))), alpha=90, maxColorValue=255 ))[testing_repetitions_at_time_t$comp_symmetry_factor ])
# 
# legend("topleft", legend = c("Asymmetric", "Symmetric"), 
#        fill = c("red","blue"), cex=0.8)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# likelihood_asymm <- dnorm(x = testing_repetitions_at_time_t$coef_var , 
#                           mean =  reg_asym$pred.y,
#                           sd = reg_asym$pred.se)
# 
# likelihood_symm <- dnorm(x = testing_repetitions_at_time_t$coef_var, 
#                          mean =  reg_symm$pred.y, sd = reg_symm$pred.se)
# 
# BF <- likelihood_symm/likelihood_asymm
# 
# 
# sum(testing_repetitions_at_time_t$comp_symmetry_factor[BF > 1]=="Symmetric")/sum(BF > 1)*100
# 
# sum(testing_repetitions_at_time_t$comp_symmetry_factor[BF < 1]=="Assymetric")/sum(BF < 1)*100
# 
# 
# predicted_final <- predict(newdata = testing_repetitions_at_time_t,
#                            object = model_final_coefvar)
# 
# 
# 
# color_BF <- sapply(BF, function(x){
#   if (x < 1){
#     ranges_menores <- c(-Inf, seq(to = 1, from = 0.1, length.out = 5))
#     which((x > ranges_menores[1:5]) & (x < ranges_menores[2:6]))
#   }else{
#     ranges_mayores <- c(seq(from = 1, to = 5, length.out = 5), Inf)
#     which((x > ranges_mayores[1:5]) & (x < ranges_mayores[2:6]))+4
#   }
# } )
# 
# red2blue <- colorRampPalette(c("red","blue"))
# 
# # png("BF_symm_prediction.png",
# #     width = 12, height = 12, units = "cm",
# #     res = 300)
# 
# layout(mat = matrix(c(1,  1, 1,0,2,0 ), 
#                     nrow = 3, 
#                     ncol = 2),
#        heights = c(1, 2, 1),
#        widths = c(5, 1))
# 
# plot(x = testing_repetitions_at_time_t$mean_no_competitors, 
#      y = testing_repetitions_at_time_t$coef_var, cex=1.3,
#      pch = c(17, 19)[testing_repetitions_at_time_t$comp_symmetry_factor ],
#      col = red2blue(n = 8)[color_BF], xlab="Mean Comp.", 
#      ylab="Final coefficient of variation")
# lines(x_range, mean_sym, col="blue", lwd=2, lty = 1)
# lines(x_range, mean_asym, col="red", lwd=2, lty = 1)
# 
# legend("topleft", legend = c("Asymmetric", "Symmetric"), 
#        pch= c(17,19), cex=0.8)
# 
# par(mar=c(2, 1, 1, 4))
# 
# pts_colorbar <- c(0, seq(to = 1, from = 0.1, length.out = 5),
#                   seq(from = 1, to = 5, length.out = 5)[-1], 100)
# 
# image(matrix(seq_along(pts_colorbar), nrow = 1),  axes=F,
#       col = red2blue(9), zlim=c(1, length(pts_colorbar)))
# 
# box(which = "plot", lty = "solid")
# 
# axis(4, at = seq(0, to =1, length=length(pts_colorbar)), 
#      labels = pts_colorbar, 
#      lwd = 0, lwd.ticks = 1, cex.axis=0.8,  las = 2)
# 
# mtext(text = "BF", side = 2, cex = 0.8, line = 1)
# 
# dev.off()
# 
# png("obs_pred_assymetry.png",
#     width = 12, height = 12, units = "cm",
#     res = 300)
# plot(y = predicted_final, 
#      x = testing_repetitions_at_time_t$coef_var, 
#      xlab = "Observed", ylab = "Predicted", pch=20, asp=1, cex=0.8,ylim=c(0,0.45),
#      col=c("red","blue")[testing_repetitions$comp_symmetry_factor[testing_repetitions$t==max(trainin_repetitions$t)]])
# 
# legend("topleft", legend = c("Asymmetric", "Symmetric"), 
#        fill = c("red","blue"), cex=0.8)
# abline(0,1, col="gray", lwd=2, lty=1)
# points(y = predicted_final, 
#        x = testing_repetitions_at_time_t$coef_var, 
#        pch=20, asp=1, cex=0.8, 
#        col=c("red","blue")[testing_repetitions_at_time_t$comp_symmetry_factor])
# 
# dev.off()
# 
# 
# error <- predicted_final-testing_repetitions_at_time_t$coef_var
# 
# symmetry_lab <- testing_repetitions_at_time_t$comp_symmetry_factor
# 
# sum(symmetry_lab[error<0]=="Assymetric")/sum(error<0)
# 
# sum(symmetry_lab[error>0]=="Symmetric")/sum(error>0)
# 
# ### hay erroe negaativo en asimetria, osea, la variacion final es mayor
# # que la predicha. 
# 
# 
# plot(variable[symmetry_lab=="Assymetric"], 
#      error[symmetry_lab=="Assymetric"], ylim=c(-0.3,0.3))
# abline(h=0)
# 
# plot(variable[symmetry_lab=="Symmetric"], 
#      error[symmetry_lab=="Symmetric"], ylim=c(-0.3,0.3))
# abline(h=0)
# 
# ##################################################
# 
# # create regression tree that predicts if on the next timestep the variation
# # will increase or decrease
# 
# difference_critical <- 1e-4
# t0 <- trainin_repetitions[trainin_repetitions$t != max(trainin_repetitions$t), ]
# 
# # coef_vart1 es el coef de var en el siguiente tiempo
# t0$coef_vart1<-trainin_repetitions$coef_var[trainin_repetitions$t != min(results_over_time$t)]
# t0$dcoef_Var <- t0$coef_vart1 - t0$coef_var
# 
# # plot(t0$coef_var, t0$dcoef_Var, pch=19, cex=0.8, 
# #      col=heat.colors(8)[as.factor(round(t0$mean_no_competitors))])
# # 
# # abline(h=0, col="red", lwd=2)
# 
# t0$var_will_increase <- as.factor(sapply(1:nrow(t0), function(x){
#   if (t0$dcoef_Var[x] > difference_critical) 
#     "inc" else if (t0$dcoef_Var[x] < -difference_critical)
#       "dec" else
#         "flat"
# }))
# 
# 
# #### repito lo mismo pero para el grupo de prueba
# 
# test_t0 <- testing_repetitions[testing_repetitions$t != max(testing_repetitions$t), ]
# 
# # coef_vart1 es el coef de var en el siguiente tiempo
# test_t0$coef_vart1<-testing_repetitions$coef_var[testing_repetitions$t != min(results_over_time$t)]
# test_t0$dcoef_Var <- test_t0$coef_vart1 - test_t0$coef_var
# 
# 
# test_t0$var_will_increase <- as.factor(sapply(1:nrow(test_t0), function(x){
#   if (test_t0$dcoef_Var[x] > difference_critical) 
#     "inc" else if (test_t0$dcoef_Var[x] < -difference_critical)
#       "dec" else
#         "flat"
# }))
# 
# 
# for (time_of_interest_delta in time_limits-1){
#   training_deltacoef_model <- t0[t0$t == time_of_interest_delta, ]
#   testing_deltacoef_model <- test_t0[test_t0$t == time_of_interest_delta,]
#   
#   training_deltacoef_model$comp_symmetry_factor <- factor(sapply(training_deltacoef_model$comp_symmetry, 
#                                                                  function(x) if (x>1) "Assymetric" else "Symmetric"))
#   
#   testing_deltacoef_model$comp_symmetry_factor <- factor(sapply(testing_deltacoef_model$comp_symmetry, 
#                                                                 function(x) if (x>1) "Assymetric" else "Symmetric"))
#   delta_coef_tree <- tree(comp_symmetry_factor ~ dcoef_Var,
#                           data=training_deltacoef_model)
#   summary(delta_coef_tree) 
#   
#   
#   tree_accuracy_this_rep <- sum(testing_deltacoef_model$comp_symmetry_factor == predict(object = delta_coef_tree, 
#                                                                                         newdata =  testing_deltacoef_model,
#                                                                                         type="class"))/nrow(testing_deltacoef_model)
#   print(round(c(time_of_interest_delta, tree_accuracy_this_rep),2))
# }
# 
# png("coefvar_stepwise_change_longer_simuls.png",
#     width = 12, height = 12, units = "cm",
#     res = 300)
# 
# plot(1, xlim=c(0, max(t0$t)), ylim=c(0,0.4), xlab="t", ylab="Coef. Var")
# hu<-1
# simuls <- c(2 ,  3 ,  4,   6,   9,  10 , 13 , 16)
# for (fre in simuls){
#   cvar <- t0$coef_var[t0$re==fre]
#   t <- t0$t[t0$re==fre]
#   type_inc <- t0$var_will_increase[t0$re==fre]
#   
#   for (se in 1:(length(t)-1)){
#     if (type_inc[se] =="inc"){ colo <- "red"} else if(type_inc[se] =="dec"){
#       colo <- "blue" } else { colo <- "green"}
#     lines(c(t[se], t[se+1]), c(cvar[se], cvar[se+1]),
#           col = colo, lwd=2)
#   }
#   
#   # predicted_max_coefvar <- predict(newdata = t0[t0$re==fre,][1,],
#   #                                  object = model_final_coefvar)
#   #abline(h=predicted_max_coefvar, lwd=2, col=rainbow(length(simuls))[hu])
#   #lines(t, cvar, col = "black", lwd=1)
#   
#   #text(t[length(t)], cvar[length(cvar)], labels = fre, col=rainbow(length(simuls))[hu])
#   hu<-hu+1
# }
# 
# dev.off()
# 
# 
# #### ya puedo predecir el final. ahora usar el final para predecir el siguiente timestep
# 
# 
# 
# colnames(t0)
# will_increase_tree <- tree(var_will_increase~coef_var+n_overcrowding_plants+
#                              mean_no_competitors+total_biomass+
#                              comp_symmetry+max_initial_sz, 
#                            data=t0)
# 
# summary(will_increase_tree)
# sum(test$var_will_increase == predict(will_increase_tree, 
#                                       test, type="class"))/nrow(test)*100
# 
# 
# plot(cv.tree(will_increase_tree))
# 
# will_increase_tree_prune <- prune.tree(will_increase_tree, best = 10)
# summary(will_increase_tree_prune)
# 
# sum(test$var_will_increase == predict(will_increase_tree_prune, 
#                                       test, type="class"))/nrow(test)*100
# 
# 
# plot(will_increase_tree_prune)
# text(will_increase_tree_prune)
# 
# 
# 
# plot(1, xlim=c(0, max(t0$t)), ylim=c(0,0.5), xlab="t", ylab="Coef. Var")
# 
# for (i_recons in 1:20){
#   
#   cvar <- t0$coef_var[t0$re==i_recons]
#   t <- t0$t[t0$re==i_recons]
#   type_inc <- t0$var_will_increase[t0$re==i_recons]
#   
#   for (se in 1:(length(t)-1)){
#     predc <- predict(will_increase_tree_prune, t0[t0$re==i_recons, ][se,], type="class")
#     actual <- type_inc[se]
#     if (actual=="inc"){
#       color_pred <- if (predc==actual) "red" else "pink"
#     } else if (actual=="dec"){
#       color_pred <- if (predc==actual) "blue" else 	"lightblue"
#     } else {
#       color_pred <- if (predc==actual) "green" else "lightgreen"
#     }
#     
#     linea <- if (predc == actual) 1 else 3
#     
#     lines(c(t[se], t[se+1]), c(cvar[se], cvar[se+1]), 
#           col = color_pred , lwd=2, lty = linea)
#   }
#   
# }
# 
# 
# # # plot of regression tree 
# # coef_var_tree_discrete$frame$var <- gsub(x = coef_var_tree_discrete$frame$var, pattern = "coef_var", replacement = "Coefficient of variation")
# # coef_var_tree_discrete$frame$var <- gsub(x = coef_var_tree_discrete$frame$var, pattern = "mean_no_competitors", replacement = "Mean No. Comp.")
# # coef_var_tree_discrete$frame$var <- gsub(x = coef_var_tree_discrete$frame$var, pattern = "sd_no_competitors", replacement = "SD No. Comp. growth rate")
# # levels(coef_var_tree_discrete$frame$yval) <- c("Asymmetric", "Symmetric")
# # 
# # png(paste("classifciation_tree_", t_lim, "_tlim.png", sep = ""),
# #     width = 27, height = 18, units = "cm",
# #     res = 300)
# # 
# # plot(coef_var_tree_discrete)
# # text(coef_var_tree_discrete, pretty = 0)
# # mtext(text =   bquote(.(round(tree_accuracy_this_rep*100, 2))~"% Cross-validation accuracy "), side = 1, line = 3)
# # 
# # dev.off()
# 
# 
# # # to  plot the predictive distribution of a specific case
# # for (case_want in c(6)){
# # curve(from= 0, to = 0.5, extraDistr::dlst(x, df = df_sym, mu = mean_sym[case_want],
# #                                           sigma = sqrt(sigmasym2[case_want])), col="blue",
# #       main=case_want, ylab = "Probability density", xlab = "Final Coef. Var.", ylim=c(0, 5))
# # # curve(from= 0, to = 2, dnorm(x, mean = mean_sym[case_want], 
# # #                              sd = sqrt(sigmasym2[case_want])),col = 'blue'  , add=T, lty=2)
# # curve(from= 0, to = 0.5, extraDistr::dlst(x, df = df_asym, mu = mean_asym[case_want],
# #                                           sigma = sqrt(sigmaasym2[case_want])), add=T, col="red" )
# # # curve(from= 0, to = 2, dnorm(x, mean = mean_asym[5], sd = sqrt(sigmaasym2[case_want])),
# # #       col = 'red'  , add=T, lty=2)
# # abline(v = testing_repetitions_at_time_t$coef_var[case_want], 
# #        col = c("red","blue")[testing_repetitions_at_time_t$comp_symmetry_factor[case_want]], lwd=2)
# # 
# # points(testing_repetitions_at_time_t$coef_var[case_want], 0,
# #        pch = if ((testing_repetitions_at_time_t$comp_symmetry_factor[case_want] == "Symmetric" &
# #                  BF[case_want] > 1) |
# #                  (testing_repetitions_at_time_t$comp_symmetry_factor[case_want] == "Assymetric" &
# #                                       BF[case_want] < 1 ) )  19 else 4, col = "black", cex=4)
# # 
# # legend("topright", legend = c("Asymmetric", "Symmetric"), 
# #        fill = c("red","blue"), cex=0.8)
# # }


