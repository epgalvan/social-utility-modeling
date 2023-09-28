################ 1_4_0

### Create Construct Value Functions

payout_maximization = function(investment, multiplier, returned){
  return(((investment * multiplier) - returned)/(investment * multiplier))
}
inequity = function(investment, multiplier, returned, endowment){
  return(((investment * multiplier - returned)/(investment * multiplier - endowment - investment))**2)
}
guilt = function(investment, believed_multiplier, returned, multiplier){
  return((((investment * believed_multiplier)/2 - returned)/(investment * multiplier))**2)
}

################ 1_5_0

### Preallocating and Defining
trialList = data.frame(Investment = rep(seq(1, 10, 1), times = 6), Multiplier = rep(c(2, 4, 6), each = 20), Believed_Multiplier = rep(4, 60), Endowment = rep(10, 60))

utility = function(theta, phi, guilt, inequity, payout){
  return(theta*payout + (1-theta)*min(guilt + phi, inequity - phi))
}

freeParameters = data.frame(theta = rep(seq(0, 0.5, 0.005), each = 101),
                            phi = rep(seq(-0.1, 0.1, 0.002), times = 101))

predictions = data.frame()

### Define All Loops

non_specific = rep(0, length(freeParameters[,1]))

for (i in 1:length(freeParameters[,1])){
  Theta = freeParameters[i,1]
  Phi = freeParameters[i,2]
  
  for (k in 1:length(trialList[,1])){
    I = trialList[k, 1]
    M = trialList[k, 2]
    B = trialList[k, 3]
    E = trialList[k, 4]
    Choices = seq(0, (I * M), 1)
    
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(theta = Theta,
                           phi = Phi,
                           guilt = guilt(I, B, Choices[n]),
                           inequity = inequity(I, M, Choices[n], E),
                           payout = payout_maximization(I, M, R))
    }
    correct_choice = which(Utility == max(Utility))
    if (length(correct_choice) > 1){
      correct_choice = correct_choice[sample(correct_choice, 1)]
      non_specific[i] =+ 1
    }
    predictions[i,k] = Choices[correct_choice]
  }
}

################ 1_6_0

### Objective Functions

obj_function = function(params, decisions, method = "OLS") {
  Theta = params[1]
  Phi = params[2]
  
  predicted_utility = vector('numeric', length(trialList[,1]))
  chosen = decisions + 1
  for (k in 1:length(trialList[,1])){
    I = trialList[k, 1]
    M = trialList[k, 2]
    B = trialList[k, 3]
    E = trialList[k, 4]
    Choices = seq(0, (I * M), 1)
    
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(Theta, Phi, guilt(I, B, Choices[n], M), inequity(I, M, Choices[n], E), payout_maximization(I, M, Choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[chosen[k]]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
}

library(pracma)

initial_params = c(0, 0)
lower_bounds = c(0, -0.1)
upper_bounds = c(0.5, 0.1)
theta_recovered = vector('numeric', 11**2)
phi_recovered = vector('numeric', 11**2)
theta_true = rep(seq(0, 0.5, 0.05), each = 11)
phi_true = rep(seq(-0.1, 0.1, 0.02), times = 11)

for (i in 1:length(theta_true)) {
  this_idx = which(theta_true[i] == freeParameters$theta & phi_true[i] == freeParameters$phi)
  result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   decisions = as.numeric(predictions[this_idx,]))
  
  theta_recovered[i] = result$par[1]
  phi_recovered[i] = result$par[2]
}

### Verifying Free Parameter Recovery Process

library(ggplot2)

distance = (2*(theta_recovered - theta_true))**2 + (5*(phi_recovered - phi_true))**2
qplot(x = theta_true, y = phi_true, color = distance, size = distance, geom = 'point') + scale_radius(limits=c(0, sqrt(2)), range=c(0, 20))

distance_new = (2*(theta_recovered - theta_true))**2 + (5*(0.5-theta_true)*(phi_recovered - phi_true))**2
qplot(x = theta_true, y = phi_true, color = distance_new, size = distance_new, geom = 'point') + scale_radius(limits=c(0, sqrt(2)), range=c(0, 20))

