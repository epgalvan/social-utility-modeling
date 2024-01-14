################ 1_4_0

### Create Construct Value Functions

harm = function(shocksThisChoice, shocksOtherChoice){
  return(shocksThisChoice - shocksOtherChoice)
}

payout = function(moneyThisChoice, moneyOtherChoice){
  return(moneyThisChoice - moneyOtherChoice)
}

################ 1_5_0

### Preallocating and Defining Functions, TrialList, and Parameters

trialList = data.frame(Default = rep(c(1, 2), each = 100),
                       MoneyA = 10,
                       MoneyB = rep(rep(seq(11, 20), each = 10), times = 2),
                       ShocksA = rep(rep(seq(0, 9), times = 10), times = 2),
                       ShocksB = 10)

utility = function(Payout, Harm, kappa, lambda){
  if (Payout < 0) {LM = lambda} else {LM = 1}
  if (Harm > 0) {LS = lambda} else {LS = 1}
  return((Payout * kappa * LM) - (Harm * (1 - kappa) * LS))
}

freeParameters = data.frame(kappa = rep(seq(0, 1, 0.1), 11), #ranging from 0 to 1, the inverse of kappa has the same range as kappa
                            lambda = rep(seq(1, 3, 0.2), 11), #loss aversion is 2.25 according to CPT, 1 is no loss aversion
                            gamma = sample(seq(0, 2, 0.01), 121), #completely random to low stochasticity
                            epsilon = 0) #no additional, non-task related noise
predictions = data.frame()

### Determine Predictions

for (i in 1:length(freeParameters[,1])){
  Kappa = freeParameters[i,1]
  Lambda = freeParameters[i,2]
  for (k in 1:length(trialList[,1])){
    shocksThis = c(trialList$ShocksA[k], trialList$ShocksB[k])[trialList$Default[k]]
    shocksAlternative = c(trialList$ShocksB[k], trialList$ShocksA[k])[trialList$Default[k]]
    moneyThis = c(trialList$MoneyA[k], trialList$MoneyB[k])[trialList$Default[k]]
    moneyAlternative = c(trialList$MoneyB[k], trialList$MoneyA[k])[trialList$Default[k]]
    
    # Just Added
    Utility = utility(Payout = harm(shocksThis, shocksAlternative),
                      Harm = payout(moneyThis, moneyAlternative),
                      kappa = Kappa,
                      lambda = Lambda)
    if (Utility > 0) {
      predictions[i,k] = 1 #choose alternative
    } else if (Utility < 0) {
      predictions[i,k] = 0 #don't choose alternative
    } else {
      predictions[i,k] = sample(c(0, 1), 1) #random
    }
  }
}

################ 1_6_0

### Objective Functions

obj_function_oneCondition = function(params, df, method = "MLE") { #we want 2 kappa parameters (1 per condition), but we don't yet have conditions
  Kappa = params[1]
  Lambda = params[2]
  Gamma = params[3]
  Epsilon = params[4]
  
  ProbCorrect = vector('numeric', length(df[,1]))
  for (k in 1:length(df[,1])){
    moneyDefault = c(df[k, 2], df[k, 3])[df[k, 1]]
    moneyAlternative = c(df[k, 3], df[k, 2])[df[k, 1]]
    shocksDefault = c(df[k, 4], df[k, 5])[df[k, 1]]
    shocksAlternative = c(df[k, 5], df[k, 4])[df[k, 1]]
    
    DiffUtility = utility(Payout = harm(shocksAlternative, shocksDefault),
                          Harm = payout(moneyAlternative, moneyDefault),
                          kappa = Kappa,
                          lambda = Lambda)
    ProbAlternative = ((1)/(1 + exp(-1 * Gamma * DiffUtility))) * (1 - (2 * Epsilon)) + Epsilon
    ProbCorrect[k] = c((1 - ProbAlternative), ProbAlternative)[(df[k, 6] + 1)]
  }
  if (method == "OLS"){
    return(sum((ProbAlternative - df[, 6])**2))
  } else if (method == "MLE"){
    return(-1 * sum(ProbCorrect, mean = 1, log = TRUE))
  }
}

### Optimizers

library(pracma)

initial_params = c(0.5, 2, 1, 1)
lower_bounds = c(0, 1, 0, 0)
upper_bounds = c(1, 5, 10, 2)
freeParameters$kappaRecovered = 0
freeParameters$lambdaRecovered = 0
freeParameters$gammaRecovered = 0
freeParameters$epsilonRecovered = 0

for (i in 1:length(freeParameters$kappa)) {
  trialList$Predictions = as.numeric(predictions[i,])
  result = fmincon(obj_function_oneCondition, x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   df = trialList)
  
  freeParameters$kappaRecovered[i] = result$par[1]
  freeParameters$lambdaRecovered[i] = result$par[2]
  freeParameters$gammaRecovered[i] = result$par[3]
  freeParameters$epsilonRecovered[i] = result$par[4]
}

### Verifying Free Parameter Recovery Process

library(ggplot2)

#how well are we able to recover kappa and lambda?
qplot(data = freeParameters, x = kappa, y = kappaRecovered) + geom_smooth()
qplot(data = freeParameters, x = lambda, y = lambdaRecovered) + geom_smooth()

#does our inverse temperature parameter (gamma) actual capture variance in stochasticity
qplot(data = freeParameters, x = gamma, y = gammaRecovered) + geom_smooth(method = 'lm')

#is epsilon capturing loss aversion? remember, epsilon is set to 0 so it should be a nonfactor
qplot(data = freeParameters, x = lambda, y = epsilonRecovered)

