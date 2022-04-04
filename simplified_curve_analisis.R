#library(animation)
results_over_time <- read.csv("LHS_sampling_results_over_time_1500res_20ws_25_seed.csv", 
                              stringsAsFactors = F)



######## Curve classification, identify what sort of temporal patterns doies
# the coef var has?

clasfify_Curve <- function(cvar, t, full_classification = TRUE){
  type_curve <- NULL
  
  #these are indexes
  init_t <- which(t==0)
  final_t <- length(t)
  
  final_cvar <- cvar[final_t]
  initial_cvar <- cvar[init_t]
  
  significativos <- abs(c(NA, cvar[-length(cvar)]) - cvar) > 1e-4
  
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
    for (i in inflpts[1]){
      #type_curve <- paste(type_curve, if (i) "_convex" else "_concave", sep="")}
      type_curve <- paste(type_curve, if (i) "_c" else "_c", sep="")}
  }

  type_curve
}

quiero <- 92
clasfify_Curve(t=results_over_time$t[results_over_time$re==quiero],
                 cvar = results_over_time$coef_var[results_over_time$re==quiero], 
               full_classification = T)

plot(x=results_over_time$t[results_over_time$re==quiero],
              y = results_over_time$coef_var[results_over_time$re==quiero], type="l")

n_sims <- max(results_over_time$re)

curves_description <- data.frame(
  type = character(n_sims),
  world_reachability= numeric(n_sims),
  max_initial_sz = numeric(n_sims),
  overcrowding = factor(numeric(n_sims), 
                        levels = c(0:max(results_over_time$n_overcrowding_plants),
                                   min(results_over_time$n_overcrowding_plants): -1)),
  comp_symmetry= numeric(n_sims),
  max_growth_rate= numeric(n_sims),
  indiviudal_var_growth_rate= numeric(n_sims), stringsAsFactors = F)

for (q in 1:n_sims){
  one_curve <- results_over_time[results_over_time$re==q, ]
  
  curves_description[q, ] <- c(clasfify_Curve(cvar = one_curve$coef_var, 
                                              t = one_curve$t, 
                                              full_classification = T)
                               , one_curve[1, 7:12])
  }

table(curves_description$type)


barplot(table(curves_description$type))

plot(1, xlim=c(0,30), ylim=c(0,0.5))
hu<-1
simuls <- which(curves_description$type == 'inc_c_c')[1:20]
for (fre in simuls){
  cvar <- results_over_time$coef_var[results_over_time$re==fre]
  t <- results_over_time$t[results_over_time$re==fre]
  #plot(t, cvar, type = "b", main=bquote(.(fre)))
  lines(t, cvar, col = rainbow(length(simuls))[hu], lwd=2)
  text(t[31], cvar[31], labels = fre, col=rainbow(length(simuls))[hu])
  hu<-hu+1
}


library(tree)

curves_description$type <- factor(curves_description$type)
curves_description$overcrowding <- as.numeric(curves_description$overcrowding)

names_frequent <- names(table(curves_description$type)[table(curves_description$type) > 20])

curves_description_frequents <- curves_description[curves_description$type %in% names_frequent,  ]

# re train with cross validation
set.seed(26)

test_index <- sample(x = 1:nrow(curves_description_frequents), size = 1000, replace = F)
# test_index <-c(
#   sample(which(curves_description_frequents$type == "inc_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "dec_c"), size = 160, replace = F),
#   sample(which(curves_description_frequents$type == "inc"), size = 80, replace = F),
#   sample(which(curves_description_frequents$type == "dec"), size = 80, replace = F))


train <- curves_description_frequents[test_index, ]
test <- curves_description_frequents[-test_index, ]

coef_var_tree_discrete <- tree(type ~ world_reachability+
                                 max_initial_sz+
                                 overcrowding+
                                 comp_symmetry+
                                 max_growth_rate+
                                 indiviudal_var_growth_rate,
                               data=train)


plot(coef_var_tree_discrete)
text(coef_var_tree_discrete)
summary(coef_var_tree_discrete)

#accuracy of first model
sum(test$type == predict(coef_var_tree_discrete, test, type="class"))/nrow(test)

plot(cv.tree(coef_var_tree_discrete))

coef_var_tree_discrete_prune <- prune.tree(coef_var_tree_discrete, best = 8)
summary(coef_var_tree_discrete_prune)

# accudary of second
sum(test$type == predict(coef_var_tree_discrete_prune, test, type="class"))/nrow(test)


coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_initial_sz", replacement = "Initial size variation")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "world_reachability", replacement = "Spatial disarrangement")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "overcrowding", replacement = "Overcrowding")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "max_growth_rate", replacement = "Maximum growth rate")
coef_var_tree_discrete_prune$frame$var <- gsub(x = coef_var_tree_discrete_prune$frame$var, pattern = "comp_symmetry", replacement = "Competition symmetry")


levels(coef_var_tree_discrete_prune$frame$yval) <- c("Dec.", "Dec. w dip",
                                                     "Inc.", "Inc. w dip")

png(paste("classifciation_tree_inc_dec_dips.png", sep = ""),   
    width = 20, height = 12, units = "cm",
    res = 300)

plot(coef_var_tree_discrete_prune)
text(coef_var_tree_discrete_prune, pretty = 0)

dev.off()



# cuales falla laclasificacion

fr <-test_index[test$type != predict(coef_var_tree_discrete_prune, test, type="class")]

fr <- fr[curves_description$type[fr]=="inc"]
fr

table(curves_description_frequents$type[fr])

# 
# curves_description <- curves_description[curves_description$simulation!=0, ]
# 
# curves_description[, -8] <- data.frame(scale(curves_description[, -8]))
# 
# #overcrowding+world_reachability+max_initial_sz +comp_symmetry+max_growth_rate+indiviudal_var_growth_rate
# 
# 
# # value of the inflection point
# lm_inf_pt <- lm(y_inflec ~ overcrowding+world_reachability+max_initial_sz +indiviudal_var_growth_rate,
#                 data = curves_description)
# 
# 
# summary(lm_inf_pt)
# 
# #position in time of inflection point
# lm_y_infl_t <- lm(y_infl_t ~ overcrowding+world_reachability+max_initial_sz +max_growth_rate
# ,
#                  data = curves_description)
# 
# summary(lm_y_infl_t)
# 
# # variation at the end
# lm_yfinal <- lm(yfinal ~ overcrowding+world_reachability+max_initial_sz +comp_symmetry
# ,
#                 data = curves_description)
# summary(lm_yfinal)
# 
# # variation at start
# lm_initial_y <- lm(yinit ~  max_initial_sz -1,
#                    data = curves_description)
# summary(lm_initial_y)
# 


library(animation)
set.seed(24)
saveGIF({
for (q in sample(which(curves_description_frequents$type=="dec_c"), 30)){

  one_curve <- results_over_time[results_over_time$re==q, ]

  final_t <-max(one_curve$t)
  final_y <- one_curve$coef_var[one_curve$t==final_t]

  initial_y <- one_curve$coef_var[one_curve$t==0]

  plot(one_curve$t, one_curve$coef_var, type="l", ylim=c(0,0.4),
       col="black", main=q, xlab="t",
       ylab = "Coefficient of variation", lwd=3)

  points(final_t, final_y, pch=19, cex=1.5)
  points(0, initial_y,  pch=19, cex=1.5)
  #points(infl_t, infl_pt, col="purple", pch=19, cex=1.5)
  #abline(v=infl_t, lty=2, lwd=1.5)

  }
  }, movie.name = "coef_var_inc_examples.gif", interval = 2,
  ani.width = 480, ani.height = 480)
# # 
# 
# ########### que esta pasando con las que no convergen
# 
# results_over_time_weirds <- results_over_time[results_over_time$re %in% not_infl, ]
# 
# #7  15  17  19  20  22  31  32  37  39  40  57  58  59  63
# 
# 
# plot(results_over_time$t,
#      results_over_time$coef_var,
#      type="n", xlab="t", ylab="Coefficient of variation")
# for (re in not_converged){
#   lines(results_over_time$t[results_over_time$re==re],
#         results_over_time$coef_var[results_over_time$re==re],
#         col = "gray", lwd=2)
# }
# 
# 
# 
# #png("parameter_space_normal_vs_conspicuous.png",
# #    width = 24, height = 12, units = "cm",
# #    res = 300)
# 

 results_over_time_weirds <- results_over_time[results_over_time$re %in% fr, ]
 
{par(mfrow=c(2,3),  mar = c(4,4,2,1),  cex.lab=1.1, oma = c(4,1,1,1))
  variables <- c('Spatial disarrangement' = "world_reachability",
                 'Initial size variation' = "max_initial_sz",
                 'Overcrowgind' = "n_overcrowding_plants",
                 'Competition symmetry' = "comp_symmetry",
                 'Maximum growth rate' = "max_growth_rate",
                 'Ind. variation in growth rate' = "indiviudal_var_growth_rate")
  for (va in seq_along(variables)){
    hist(results_over_time_weirds[[variables[va]]][results_over_time$t==1],
         col = rgb(0,0,1,0.5), xlab=names(variables)[va],main="", freq = F,
         breaks = seq(min(results_over_time[[variables[va]]]),
                      max(results_over_time[[variables[va]]]), length.out=10))
    hist(results_over_time[[variables[va]]][results_over_time_weirds$t==1],
         col = rgb(0,1,0, 0.5), add=T, freq = F,
         breaks = seq(min(results_over_time[[variables[va]]]),
                      max(results_over_time[[variables[va]]]), length.out=10))
  }

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom',legend = c("All simulations", "Simulations with increasing variation"),
         fill=c("green","blue"), xpd = TRUE, horiz = TRUE, cex = 1.5)
  # xpd = TRUE makes the legend plot to the figure
}

 

 
 table(curves_description$type)
 barplot(table(curves_description$type))
 


bayesian_posterior_parameter <- function(values, total_values, 
                                         ranges_init = NULL, prior = NULL, 
                                         n_categ = 10){
  N <- length(total_values)
  
  if(is.factor(values)){
    
    ranges_init <- sort(as.numeric(levels(values)))
    prior <- sapply(seq_along(ranges_init), function(x){ sum(total_values == ranges_init[x])} )
    prior <- prior/N
    
    
    likelyhood <- sapply(seq_along(ranges_init), function(x){ sum(values == ranges_init[x])  } )
    likelyhood <- likelyhood/length(values)
    
    return(list(posterior = (prior*likelyhood)/sum(prior*likelyhood),
         prior = prior,
         param_vals = ranges_init))
  }
  
  if (is.null(ranges_init)){
    ranges_init <- seq(from = min(total_values), to = max(total_values), length.out = n_categ)
  }

  
 prior <- sapply(1:(n_categ-1), function(x){ sum(total_values >= ranges_init[x] & total_values < ranges_init[x+1])  } )
 prior[n_categ] <- sum(total_values >=  ranges_init[n_categ] )
 prior <- prior/N


 likelyhood <- sapply(1:(n_categ-1), function(x){ sum(values >= ranges_init[x] & values < ranges_init[x+1])  } )
 likelyhood[n_categ] <- sum(values >=  ranges_init[n_categ] )
 likelyhood <- likelyhood/length(values)

list(posterior = (prior*likelyhood)/sum(prior*likelyhood),
     prior = prior,
     param_vals = ranges_init)
}


range_for_assymetry <- c(seq(from = 0, to= 1, length.out = 7)[-7], seq(from = 1, to= 100, length.out = 6)) 


type_simul <- "inc_c"

fr <- which(curves_description$type==type_simul)
subset <- curves_description[1:nrow(curves_description) %in% fr ,]
#subset$overcrowding <- as.factor(subset$overcrowding)

png(paste("parameter_posterior", type_simul, ".png", sep="") , 
    width = 24, height = 12, units = "cm",
    res = 300)
{par(mfrow=c(2,3),  mar = c(4,4,2,1),  cex.lab=1.1, oma = c(4,1,1,1))
  variables <- c('Spatial disarrangement' = "world_reachability",
                 'Initial size variation' = "max_initial_sz",
                 'Overcrowding' = "overcrowding",
                 'Competition symmetry' = "comp_symmetry",
                 'Maximum growth rate' = "max_growth_rate",
                 'Ind. variation in growth rate' = "indiviudal_var_growth_rate")
  
  for (va in seq_along(variables)){
    
    if (variables[va] =="comp_symmetry"){
      prob <- bayesian_posterior_parameter(values = subset$comp_symmetry, 
                                           total_values =  curves_description$comp_symmetry, 
                                           ranges_init =  range_for_assymetry, 
                                           n_categ = length(range_for_assymetry))
      
      
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
  
  mtext("Increasing variation with dip", side = 3, line = - 2,
        outer = TRUE)
}
dev.off()

