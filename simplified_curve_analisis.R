#library(animation)
results_over_time <- read.csv("LHS_results_world_reachability_max_initial_sz_n_overcrowding_plants_comp_symmetry_4_reps_20ws_50tsteps_23_seed.csv",
                              stringsAsFactors = F)

# results_over_time <- read.csv("tree_steps/tree4_variables_world_reachability_max_initial_sz_n_overcrowding_plants_comp_symmetry_400_reps.csv", 
#                              stringsAsFactors = F)


############## relationship between competition and variation

data_with_competiton <- results_over_time[results_over_time$mean_competiton < 1, ]


png("coefvar_competiton_per_time.png",
    width = 12, height = 12, units = "cm", res = 300)

plot(data_with_competiton$coef_var, 
     data_with_competiton$mean_competiton, pch=".", 
     ylab="Mean competition coefficient", 
     xlab=" Coefficient of variation",
     col = heat.colors(31)[as.factor(data_with_competiton$t)])

model_competiton<- lm(mean_competiton~coef_var, 
                      data=data_with_competiton)

abline(model_competiton, col="black", lwd=1.5)

dev.off()

summary(model_competiton)

plot(model_competiton)
######## Curve classification, identify what sort of temporal patterns doies
# the coef var has?


clasfify_Curve <- function(cvar, t, full_classification = TRUE){
  #these are indexes
  init_t <- which(t==0)
  final_t <- length(t)
  
  final_cvar <- cvar[final_t]
  initial_cvar <- cvar[init_t]
  
  significativos <- abs(c(NA, cvar[-length(cvar)]) - cvar) > 1e-5
  
  if (all(!(significativos[-1]))) {return("flat")}
  
  if ( final_cvar > initial_cvar) {
    type_curve <- "inc"
    aumenta_coefvar <- c(NA, cvar[-length(cvar)]) < cvar  
    
    }else{
      aumenta_coefvar <- c(NA, cvar[-length(cvar)]) > cvar 
    type_curve <- "dec"}
  state <- aumenta_coefvar[2]
  
  inflpts <- c()
  
  for (i in 3:length(aumenta_coefvar)){
    
    if (state != aumenta_coefvar[i] & significativos[i]){
      state <- aumenta_coefvar[i]
      inflpts <- c( inflpts, state )# esto va a decir si es concavo o convexo la diferencia
    }
  }
  #  print(aumenta_coefvar)
    #inflpts <- if (length(inflpts)>0) inflpts[seq(from= 1, to=length(inflpts), by=2)]

  if (full_classification) {
    for (i in inflpts){
      type_curve <- paste(type_curve, if (i) "_convex" else "_concave", sep="")}
      #type_curve <- paste(type_curve, if (i) "_c" else "_c", sep="")}
  }

  type_curve
}

quiero <- 27
clasfify_Curve(t=results_over_time$t[results_over_time$re==quiero],
                 cvar = results_over_time$coef_var[results_over_time$re==quiero],
               full_classification = T)
results_over_time$n_overcrowding_plants[results_over_time$re==67][1]
plot(x=results_over_time$t[results_over_time$re==quiero],
              y = results_over_time$coef_var[results_over_time$re==quiero], type="l")

n_sims <- max(results_over_time$re)

curves_description <- data.frame(
  type = character(n_sims),
  mean_potential_competitors  = numeric(n_sims),
  sd_potential_competitors  = numeric(n_sims),
  world_reachability= numeric(n_sims),
  max_initial_sz = numeric(n_sims),
  overcrowding = numeric(n_sims),
  comp_symmetry= numeric(n_sims),
  stringsAsFactors = F)

for (q in 1:n_sims){
  one_curve <- results_over_time[results_over_time$re==q, ]

  curves_description[q, ] <- c(clasfify_Curve(cvar = one_curve$coef_var, 
                                              t = one_curve$t, 
                                              full_classification = T)
                               , one_curve[1, 8:13])
  }

table(curves_description$type)/sum(table(curves_description$type))*100



#barplot(table(curves_description$type))

#  png("coefvar_per_biomass_dec_c.png", 
#      width = 12, height = 12, units = "cm",
#      res = 300)
 plot(1, xlim=c(0,max(results_over_time$t)), ylim=c(0.05,0.3), 
      xlab="t", ylab="Coef. Var")
 hu<-1
 simuls <- which(curves_description$type == 'inc_convex')[1:10]
 for (fre in simuls){
   cvar <- results_over_time$coef_var[results_over_time$re==fre]
   t <- results_over_time$t[results_over_time$re==fre]
   #plot(t, cvar, type = "b", main=bquote(.(fre)))
   lines(t, cvar, col = rainbow(length(simuls))[hu], lwd=2)
   text(t[31], cvar[31], labels = fre, col=rainbow(length(simuls))[hu])
   hu<-hu+1
 }
# dev.off()
# 
library(tree)

curves_description$type <- factor(curves_description$type)

curves_description$type_inc_dec <- factor(sapply(curves_description$type, function(x){
  if (x %in% c("inc_c", "inc")) "inc" else "dec"
}))

curves_description$type_c_noc <- factor(sapply(curves_description$type, function(x){
  if (x %in% c("inc_c", "dec_c")) "C" else "N"
}))


curves_description$overcrowding <- as.numeric(curves_description$overcrowding)

names_frequent <- names(table(curves_description$type)[table(curves_description$type) > 20])

curves_description_frequents <- curves_description[curves_description$type %in% names_frequent,  ]

# re train with cross validation
set.seed(26)

test_index <- sample(x = 1:nrow(curves_description_frequents), 
                     size = 1/3*nrow(curves_description_frequents), replace = F)
# test_index <-c(
#   sample(which(curves_description_frequents$type == "inc_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "dec_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "inc"), size = 80, replace = F),
#   sample(which(curves_description_frequents$type == "dec"), size = 80, replace = F))


train <- curves_description_frequents[test_index, ]
test <- curves_description_frequents[-test_index, ]

coef_var_tree_discrete <- tree(type~ mean_potential_competitors+
                                 max_initial_sz +
                                 overcrowding +
                                 comp_symmetry,
                               data=train)


plot(coef_var_tree_discrete)
text(coef_var_tree_discrete)
summary(coef_var_tree_discrete)

#accuracy of first model
sum(test$type == predict(coef_var_tree_discrete, test, type="class"))/nrow(test)

# to know what type of simulations are harder to predict

wrongly_predicted <- test[test$type_inc_dec != predict(coef_var_tree_discrete, test, type="class"),]


table(curves_description$type_inc_dec)
table(wrongly_predicted$type_inc_dec)


plot(cv.tree(coef_var_tree_discrete))

coef_var_tree_discrete_prune <- prune.tree(coef_var_tree_discrete, best = 8)
summary(coef_var_tree_discrete_prune)

# accudary of second
pcntg <- sum(test$type == predict(coef_var_tree_discrete_prune, 
                                          test, type="class"))/nrow(test)*100
pcntg <- paste(round(pcntg, 1),"%", sep="" )

coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_initial_sz", replacement = "Initial size variation")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "world_reachability", replacement = "Spatial disarrangement")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "overcrowding", replacement = "Overcrowding")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_growth_rate", replacement = "Maximum growth rate")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "comp_symmetry", replacement = "Competition symmetry")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "t_potential_competitors", replacement = "Potential competitors")


levels(coef_var_tree_discrete_prune$frame$yval) <- c("Dec.",  "Inc.")

png(paste("classifciation_tree_detail_description.png", sep = ""),
    width = 27, height = 18, units = "cm",
    res = 300)

plot(coef_var_tree_discrete_prune)
text(coef_var_tree_discrete_prune, pretty = 0)

text(7.2,3000,  bquote(.(pcntg)~" Cross-validation\n predicccion accuracy "))
dev.off()
############# repeat for c or no c

set.seed(25)

test_index <- sample(x = 1:nrow(curves_description_frequents), size = 1000, replace = F)
# test_index <-c(
#   sample(which(curves_description_frequents$type == "inc_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "dec_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "inc"), size = 80, replace = F),
#   sample(which(curves_description_frequents$type == "dec"), size = 80, replace = F))


train <- curves_description_frequents[test_index, ]
test <- curves_description_frequents[-test_index, ]

coef_var_tree_discrete <- tree(type_c_noc ~ world_reachability+
                                 max_initial_sz+
                                 overcrowding+
                                 comp_symmetry,
                               data=train)


plot(coef_var_tree_discrete)
text(coef_var_tree_discrete)
summary(coef_var_tree_discrete)

#accuracy of first model
sum(test$type_c_noc == predict(coef_var_tree_discrete, test, type="class"))/nrow(test)

plot(cv.tree(coef_var_tree_discrete))

coef_var_tree_discrete_prune <- prune.tree(coef_var_tree_discrete, best = 8)
summary(coef_var_tree_discrete_prune)

# accudary of second
pcntg <- sum(test$type_c_noc == predict(coef_var_tree_discrete_prune, 
                                          test, type="class"))/nrow(test)*100
pcntg <- paste(round(pcntg, 1),"%", sep="" )


coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_initial_sz", replacement = "Initial size variation")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "world_reachability", replacement = "Spatial disarrangement")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "overcrowding", replacement = "Overcrowding")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_growth_rate", replacement = "Maximum growth rate")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "comp_symmetry", replacement = "Competition symmetry")


levels(coef_var_tree_discrete_prune$frame$yval) <- c("Dip", "No dip")

png(paste("classifciation_tree_dip_nodip.png", sep = ""),
    width = 27, height = 18, units = "cm",
    res = 300)

plot(coef_var_tree_discrete_prune)
text(coef_var_tree_discrete_prune, pretty = 0)
text(2,935,  bquote(.(pcntg)~" Cross-validation\n predicccion accuracy "))

dev.off()


# 
# library(animation)
# set.seed(24)
# saveGIF({
# for (q in sample(which(curves_description_frequents$type=="dec_c"), 30)){
# 
#   one_curve <- results_over_time[results_over_time$re==q, ]
# 
#   final_t <-max(one_curve$t)
#   final_y <- one_curve$total_biomass[one_curve$t==final_t]
# 
#   initial_y <- one_curve$total_biomass[one_curve$t==0]
# 
#   plot(one_curve$total_biomass, one_curve$coef_var, type="l", ylim=c(0,0.4),
#       xlim=c(0, 400),
#        col="black", main=q, xlab="Total Biomass",
#        ylab = "Coefficient of variation", lwd=3)
# 
#   #points(final_t, final_y, pch=19, cex=1.5)
#   #points(0, initial_y,  pch=19, cex=1.5)
#   #points(infl_t, infl_pt, col="purple", pch=19, cex=1.5)
#   #abline(v=infl_t, lty=2, lwd=1.5)
# 
#   }
#   }, movie.name = "coef_var_dec_biomass.gif", interval = 2,
#   ani.width = 480, ani.height = 480)
# 



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
    
    # if( is.factor(values)){
    #   
    #   ranges_init <- sort(as.numeric(levels(values)))
    #   prior <- sapply(seq_along(ranges_init), function(x){ sum(total_values == ranges_init[x])} )
    #   prior <- prior/N
    #   
    #   
    #   likelyhood <- sapply(seq_along(ranges_init), function(x){ sum(values == ranges_init[x])  } )
    #   likelyhood <- likelyhood/length(values)
    #   
    #   return(list(posterior = (prior*likelyhood)/sum(prior*likelyhood),
    #               prior = prior,
    #               param_vals = ranges_init))
    # }
    
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
