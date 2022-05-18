#### Analysis of simulated populations over time

results_over_time <- read.csv("IBM_res_world_reachability_max_initial_sz_n_overcrowding_plants_comp_symmetry_2100_reps_20ws_50tsteps_23_seed.csv",
                              stringsAsFactors = F)

# Time 0 is the population at initialization, before any competition occurs, so it
# is omitted from this analysis
results_over_time <- results_over_time[results_over_time$t > 0, ]

## Asses the relationship between competition and variation

lm_model <- lm(mean_competiton~coef_var,  data = results_over_time)
summary(lm_model)


my_colors <- colorRampPalette(c("#440154FF", "#482878FF", "#3E4A89FF", 
                                "#31688EFF", "#26828EFF", "#1F9E89FF",
                                "#35B779FF", "#6DCD59FF", "#B4DE2CFF",
                                "#FDE725FF"))
  
  
png("coefvar_competiton_per_time.png",
    width = 12, height = 12, units = "cm", res = 300)

layout(mat = matrix(c(1, 1, 1, 0, 2, 0),  nrow = 3, ncol = 2),
       heights = c(2,3, 2), widths = c(7,1))

plot(x = results_over_time$coef_var, 
     y = results_over_time$mean_competiton,
     pch=20, cex = 0.05,  yaxs="i", cex.lab = 1.3,
     ylab = "Mean competition coefficient", 
     xlab = "Size coefficient of variation",
     col = my_colors(length(unique(results_over_time$t)))[results_over_time$t])

# Colorbar
pts_colorbar <- 0:max(results_over_time$t)
par(mar = c(2, 0, 1, 4))
image(matrix(pts_colorbar, nrow = 1), axes = F,
      col = my_colors(length(unique(results_over_time$t))), 
      zlim = c(0, max(pts_colorbar)))
box(which = "plot", lty = "solid")

ticks_at <- seq(from = min(results_over_time$t), to = max(results_over_time$t), 
                length.out = 6)
axis(4, at = seq(0, 1, length.out = length(ticks_at)), 
     labels = pts_colorbar[ticks_at + 1], 
     lwd = 0, lwd.ticks = 1, las = 1)

mtext(text = "Time-step", side = 3, cex = 0.8, line = 1)

dev.off()

cor.test(x = results_over_time$coef_var, 
         y = results_over_time$mean_competiton)

### plots f examples of curves 
###### esta fallando lo de abajo, los ejes del paneld erecho no quedan bieeeeeeeeen
png("examples_curves_variation.png",
    width = 14, height = 12, units = "cm", res = 300)
set.seed(26)
par(mar = c(5, 4, 4, 1) + 0.1)
layout(mat = matrix(c(1, 2),  nrow = 1, ncol = 2),
       heights = c(2), widths = c(4, 1))
plot(1, xlim = c(1, max(results_over_time$t)), 
     ylim=c(0, 0.4), yaxs = "i", xaxs = "i",
     xlab="Time-step", ylab="Size Coefficient of Variation", )
hu<-1
simuls <- sample(1:max(results_over_time$re), size = 9)

for (fre in simuls){
  cvar <- results_over_time$coef_var[results_over_time$re==fre]
  t <- results_over_time$t[results_over_time$re==fre]
  lines(t, cvar, col = RColorBrewer::brewer.pal(length(simuls), "Set1")[hu], 
        lwd=2)
  hu<-hu+1
}
par(mar = c(10, 0 , 4, 3))
plot(1, xlim = c(0, 0.5), ylim=c(0, 1), 
     xlab="", ylab="", xaxt = "n", yaxt = "n", bty="n")

axis_positions <- seq(from = 0, to = 0.5, length.out = 4)
axis(side = 4, at = c(0, 1), labels = c("Lower limit", "Upper limit"), 
     las = 1, cex.axis = 0.5, tick = FALSE)
axis(side = 1, at = axis_positions, labels = c("Spatial disarrangement", 
                                               "Initial size variation",
                                               "Overcrowding",
                                               "Competition symmetry"), 
     las = 2, cex.axis = 0.8, tick = FALSE)

for (xcord in axis_positions){
  lines(x = c(xcord, xcord), y = c(0, 1), 
        lty = 1, col = "grey" )
  lines(x = c(xcord * 0.9, xcord * 1.1), y = c(1, 1), 
        lty = 1, col = "grey" )
  lines(x = c(xcord * 0.9, xcord * 1.1), y = c(0, 0), 
        lty = 1, col = "grey" )
  }
hu <- 1
for (fre in simuls){
    points(axis_positions[1], results_over_time$world_reachability[results_over_time$re==fre][1]/max(results_over_time$world_reachability),
           col = RColorBrewer::brewer.pal(length(simuls), "Set1")[hu],
           pch = "-", cex= 2)
    points(axis_positions[2], results_over_time$max_initial_sz[results_over_time$re==fre][1]/max(results_over_time$max_initial_sz),
           col = RColorBrewer::brewer.pal(length(simuls), "Set1")[hu],
           pch = "-", cex= 2)
    points(axis_positions[3], results_over_time$n_overcrowding_plants[results_over_time$re==fre][1]/max(results_over_time$world_reachability),
           col = RColorBrewer::brewer.pal(length(simuls), "Set1")[hu],
           pch = "-", cex= 2)
    points(axis_positions[4], results_over_time$comp_symmetry[results_over_time$re==fre][1]/max(results_over_time$comp_symmetry),
           col = RColorBrewer::brewer.pal(length(simuls), "Set1")[hu],
           pch = "-", cex= 2)
    hu<-hu+1
  }

dev.off()




### bayesian posterior distribution of parameters

bayesian_posterior_parameter <- function(values, 
                                         total_values,
                                         values_for_joint = NULL,
                                         total_values_for_joint = NULL,
                                         ranges_init_joint =NULL, 
                                         ranges_init = NULL, 
                                         prior = NULL, 
                                         n_categ = 10,
                                         n_categ_joint = 10){
  N <- length(total_values)
  
  if (is.null(ranges_init)){
    ranges_init <- seq(from = min(total_values),
                       to = max(total_values),
                       length.out = n_categ +1 )
    }else{
      n_categ <- length(ranges_init)
      ranges_init <- c(ranges_init, Inf)}
  
  
  
  if ( is.null(values_for_joint)) {
    
   prior <- sapply(1:n_categ, function(x){ sum(total_values >= ranges_init[x] & total_values < ranges_init[x+1])  } )
   prior <- prior/N

   likelyhood <- sapply(1:n_categ, function(x){ sum(values >= ranges_init[x] & values < ranges_init[x+1])  } )
   likelyhood <- likelyhood/N
  
  return(list(posterior = (prior*likelyhood)/sum(prior*likelyhood),
       prior = prior,
       param_vals = ranges_init[1:n_categ]))

  }else{
    
    if (is.null(ranges_init_joint)){
    ranges_init_joint <- seq(from = min(total_values_for_joint),
                         to = max(total_values_for_joint),
                         length.out = n_categ_joint + 1)
    }else{
      n_categ_joint <- length(ranges_init_joint)
      ranges_init_joint <- c(ranges_init_joint, Inf)}
    
  
    prior <- outer(X = 1:n_categ,  
                   Y = 1:n_categ_joint,
                   FUN = Vectorize(function(x, y){
                     is_x <- (total_values >= ranges_init[x]) & (total_values < ranges_init[x+1])
                     is_y <- (total_values_for_joint >= ranges_init_joint[y]) & (total_values_for_joint < ranges_init_joint[y+1]) 
                     
                     return(sum(is_x & is_y)) }))
    prior <- prior/N
    
   
    
    likelyhood <- outer(X = 1:n_categ,  
                        Y = 1:n_categ_joint,
                        FUN = Vectorize(function(x, y){
                          is_x <- (values >= ranges_init[x]) & (values < ranges_init[x+1])
                          is_y <- (values_for_joint >= ranges_init_joint[y]) & (values_for_joint < ranges_init_joint[y+1])
                          return(sum(is_x & is_y)) }))
    
    likelyhood <- likelyhood/N
    
    posterior <- (prior*likelyhood)/sum(prior*likelyhood)

    # print(n_categ_joint)
    # print(ranges_init_joint)
    # 
    # print(n_categ)
    # print(ranges_init)
    # 
    # 
    # print(dim(posterior))
    
    colnames(posterior) <- ranges_init_joint[1:n_categ_joint]
    rownames(posterior) <- ranges_init[1:n_categ]
    
    
    return(list(posterior = posterior,
                prior = prior,
                param_vals = ranges_init[1:n_categ],
                param_vals_joint = ranges_init_joint[1:n_categ_joint]))
    
  }  
}






variables <- c('Spatial disarrangement' = "world_reachability",
               'Initial size variation' = "max_initial_sz",
               'Overcrowding' = "overcrowding",
               'Competition symmetry' = "comp_symmetry"
               #'Maximum growth rate' = "max_growth_rate",
               #'Ind. variation in growth rate' = "indiviudal_var_growth_rate"
               )

table(curves_description$type)

type_simul <- "inc_c"

fr <- which(curves_description$type==type_simul)
subset <- curves_description[1:nrow(curves_description) %in% fr ,]


combinations <- combn(1:4, 2)

range_for_assymetry <- c(seq(0, 1, length.out=11), seq(1, 100, length.out=11)[-c(1, 11)])

range_for_overcrowding <- seq(from=-10, to=14, by=1)

 png(paste("parameter_posterior_joint_symm_ovcr_", type_simul, ".png", sep="") ,
     width = 18, height = 12, units = "cm",
     res = 300)
{
layout(mat = matrix(c(1, 1, 2, 2, 3,  3,  4,  4, 5, 5, 6, 6, 0,  7, 7, 0), 
                    nrow = 4, 
                    ncol = 4),
       heights = c(2, 2, 2, 2),
       widths = c(3, 3, 3, 1))

for (k in 1:ncol(combinations)){
  
  variables_for_joint <- variables[combinations[, k]]
  n_groups_paramter <- 13
  ev <- (n_groups_paramter-1)/2
    {
  
  range_this_rep <- if  (variables_for_joint[1] == "comp_symmetry"
                         ) range_for_assymetry else if (variables_for_joint[1] == "overcrowding"
                                                        ) range_for_overcrowding else NULL
      
  range_this_rep_joint <- if  (variables_for_joint[2] == "comp_symmetry"
      ) range_for_assymetry else if (variables_for_joint[2] == "overcrowding"
      ) range_for_overcrowding else NULL
      
  
      prob_joint <- bayesian_posterior_parameter(values = subset[[variables_for_joint[1]]], 
                                       total_values =  curves_description[[variables_for_joint[1]]],
                                       values_for_joint = subset[[variables_for_joint[2]]], 
                                       total_values_for_joint = curves_description[[variables_for_joint[2]]],
                                       ranges_init = range_this_rep ,
                                       ranges_init_joint = range_this_rep_joint,
                                       n_categ = 10, n_categ_joint = 10)
  max_colobar <- max(prob_joint$posterior)
  palette <- colorRampPalette(c('#ffffff','#200080'))(10)
  par(c(3, 4, 3, 1) + 0.1)
  image(prob_joint$posterior, col = palette, axes=FALSE,
        main="", zlim = c(0, max_colobar))
  
  grid(nx = nrow(prob_joint$posterior), ny = ncol(prob_joint$posterior),
       col = "gray", lty = 1)
  
  box(which = "plot", lty = "solid")
  
  axis(2, at = seq(0,1, length.out = ncol(prob_joint$posterior) ), 
       labels = round(as.numeric(colnames(prob_joint$posterior)), 2),
      lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
  
  axis(1, at = seq(0, 1, length.out=nrow(prob_joint$posterior) ), 
       labels = round(as.numeric(rownames(prob_joint$posterior)), 2),
       lwd = 0, lwd.ticks = 1, cex.axis = 0.75)
  
  mtext(1, text = names(variables_for_joint[1]), line = 2, cex=0.8)
  mtext(2, text = names(variables_for_joint[2]), line = 2, cex=0.8)
  }

}
par(mar=c(2, 1, 1, 4))

pts_colorbar <- seq(0, max_colobar, length.out = 10)
image(matrix(pts_colorbar, nrow = 1), axes=F,
      col = palette, zlim=c(0, max_colobar))
box(which = "plot", lty = "solid")

axis(4, at = seq(0, 1, length.out = length(pts_colorbar)), 
     labels = round(pts_colorbar, 2), 
     lwd = 0, lwd.ticks = 1)

mtext(text = "Probability", side = 2, cex = 0.8, line = 1)



mtext(type_simul, side = 3, line = -3,
      outer = TRUE)
 }
 
dev.off()
# 
# # png(paste("parameter_posterior", type_simul, ".png", sep="") , 
# #     width = 24, height = 12, units = "cm",
# #     res = 300)
{par(mfrow=c(2,2),  mar = c(4,4,2,1),  cex.lab=1.1, oma = c(4,1,1,1))

  for (va in c(1,2,3,4)){#seq_along(variables)){

    if (variables[va] =="comp_symmetry"){
      prob <- bayesian_posterior_parameter(values = subset$comp_symmetry,
                                           total_values =  curves_description$comp_symmetry,
                                           ranges_init =  range_for_assymetry,
                                           n_categ = length(range_for_assymetry))

      }else if (variables[va] =="overcrowding"){
        prob <- bayesian_posterior_parameter(values = subset$overcrowding,
                                             total_values =  curves_description$overcrowding,
                                             ranges_init =  range_for_overcrowding,
                                             n_categ = length(range_for_overcrowding))


    }else{
      prob <- bayesian_posterior_parameter(values = subset[[variables[va]]],
                                           total_values =  curves_description[[variables[va]]],
                                           n_categ = 12)
      }
    barplot(prob$posterior, col = rgb(0, 0 , 1, alpha = 0.3),
            names.arg = round(prob$param_vals,2),xlab=names(variables)[va],
            main="", ylab="Probability")

    barplot(prob$prior, add = T, col = rgb(0, 1 , 0, alpha = 0.3))

  }

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom',legend = c("Posterior", "Prior"),
         fill=c(rgb(red = c(0, 0), green = c(0, 1), blue = c(1, 0),
                    alpha = 0.5 )), xpd = TRUE, horiz = TRUE, cex = 1)

  mtext(type_simul, side = 3, line = - 2,
        outer = TRUE)
}
# 
# #dev.off()
# 
