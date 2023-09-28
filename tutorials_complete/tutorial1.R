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
trialList = data.frame(Investment = rep(seq(1, 10, 1), times = 6),
                       Multiplier = rep(c(2, 4, 6), each = 20),
                       Believed_Multiplier = rep(4, 60),
                       Endowment = rep(10, 60))

utility = function(theta, phi, guilt, inequity, payout){
  return(theta*payout + (1-theta)*min(guilt + phi, inequity - phi))
}

freeParameters = data.frame(theta = rep(seq(0, 0.5, 0.005), each = 101),
                            phi = rep(seq(-0.1, 0.1, 0.002), times = 101))

predictions = data.frame()

### Define All Loops

non_specific = rep(0, length(freeParameters[,1])) # Just Added This Line

for (i in 1:length(freeParameters[,1])){
  Theta = freeParameters[i,1]
  Phi = freeParameters[i,2]
  
  for (k in 1:length(trialList[,1])){
    I = trialList[k, 1]
    M = trialList[k, 2]
    B = trialList[k, 3]
    E = trialList[k, 4]
    Choices = seq(0, (I * M), 1)
    
    # Just Added
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
