library(animation)
results_over_time <- read.csv("LHS_sampling_results_over_time_800res_20ws_21seed_0.1thinn.csv", 
                              stringsAsFactors = F)



######## Curve classification, identify what sort of temporal patterns doies
# the coef var has?

clasfify_Curve <- function(cvar, t){
  type_curve <- NULL
  
  #these are indexes
  init_t <- which(t==0)
  final_t <- length(t)
  
  final_cvar <- cvar[final_t]
  initial_cvar <- cvar[init_t]
  
  significativos <- abs(c(NA, cvar[-length(cvar)]) - cvar) > 1e-4

  if (all(!significativos)) {return("flat")}
  
  if ( final_cvar > initial_cvar) {
    type_curve <- "inc"
    aumenta_coefvar <- c(NA, cvar[-length(cvar)]) < cvar  
    
    }else{
      aumenta_coefvar <- c(NA, cvar[-length(cvar)]) > cvar 
    type_curve <- "dec"}

  difer <- diff(cumsum(aumenta_coefvar[2:length(aumenta_coefvar)]))

  state <- difer[1]
  
  inflpts <- c()
  for (i in 2:length(difer)){
    if (state != difer[i] & significativos[i]){
      state <- difer[i]
      inflpts <- c( inflpts, state )# esto va a decir si es concavo o convexo la diferencia
    }
  }
  inflpts <- if (length(inflpts)>0) inflpts[seq(from= 1, to=length(inflpts), by=2)]

   # for (i in inflpts){
   #   type_curve <- paste(type_curve, if (i) "_convex" else "_concave", 
   #                       sep="")}
  type_curve
}

quiero <- 15
clasfify_Curve(t=results_over_time$t[results_over_time$re==quiero],
                 cvar = results_over_time$coef_var[results_over_time$re==quiero])

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
                                              t = one_curve$t)
                               , one_curve[1, 7:12])
  }

names(table(curves_description$type))[]
barplot(table(curves_description$type))

fr <- which(curves_description$type=="dec")
plot(1, xlim=c(0,30), ylim=c(0,0.5))
hu<-1
for (fre in fr){
cvar <- results_over_time$coef_var[results_over_time$re==fre]
t<- results_over_time$t[results_over_time$re==fre]
#plot(t, cvar, type = "b", main=bquote(.(fre)))
lines(t, cvar, col = rainbow(length(fr))[hu], lwd=2)
text(t[31], cvar[31], labels = fre, col=rainbow(length(fr))[hu])
hu<-hu+1
}

curves_description <- curves_description[
  curves_description$type %in% names(table(curves_description$type))[table(curves_description$type)>10],]


library(tree)
# library(rpart)
# 
curves_description$type <- factor(curves_description$type)
curves_description$overcrowding <- factor(curves_description$overcrowding)
# 
# 
# fit <- rpart(type ~ world_reachability+
#                max_initial_sz+
#                overcrowding,
#              method="class", data=curves_description[1:100,])
# 
# printcp(fit) # display the results
# plot(fit) # visualize cross-validation results
# text(fit)
# summary(fit) # detailed summary of splits
# 
# 
coef_var_tree_discrete <- tree(type ~ world_reachability+
                                 max_initial_sz+
                                 overcrowding+
                                 comp_symmetry+
                                 max_growth_rate+
                                 indiviudal_var_growth_rate,
                               data=curves_description)

plot(coef_var_tree_discrete)
text(coef_var_tree_discrete)
summary(coef_var_tree_discrete)

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
# # set.seed(24)
# # saveGIF({
# # for (q in sample(x = 1:max(results_over_time$re), size = 20)){
# # 
# #   one_curve <- results_over_time[results_over_time$re==q, ]
# # 
# #   final_t <-max(one_curve$t)
# #   final_y <- one_curve$coef_var[one_curve$t==final_t]
# # 
# #   initial_y <- one_curve$coef_var[one_curve$t==0]
# # 
# #   min_y <- min(one_curve$coef_var)
# #   max_y <- max(one_curve$coef_var)
# # 
# #   infl_pt <- if (min_y < initial_y) min_y else max_y
# #   infl_t <- one_curve$t[abs(one_curve$coef_var-infl_pt) < 1e-12 ]
# #   if(length(infl_t)>1){next}
# # 
# #   plot(one_curve$t, one_curve$coef_var, type="l", ylim=c(0,0.4),
# #        col="red", main=q, xlab="t",
# #        ylab = "Coefficient of variation", lwd=2.5)
# # 
# #   points(final_t, final_y, col="blue", pch=19, cex=1.5)
# #   points(0, initial_y, col="green", pch=19, cex=2)
# #   points(infl_t, infl_pt, col="purple", pch=19, cex=1.5)
# #   abline(v=infl_t, lty=2, lwd=1.5)
# # 
# #   }
# #   }, movie.name = "coef_var_curves_examples.gif", interval = 2,
# #   ani.width = 480, ani.height = 480)
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



