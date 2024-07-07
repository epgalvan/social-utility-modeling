################ 1_4_0

### Create Construct Value Functions

harm = function(shocksLarger, shocksSmaller){
  return(shocksLarger - shocksSmaller)
}

payout = function(moneyLarger, moneySmaller){
  return(moneyLarger - moneySmaller)
}

################ 1_5_0

### Preallocating and Defining Functions, TrialList, and Parameters

trialList = data.frame(Left = rep(c(T, F), 100),
                       MoneyA = rep(rep(seq(11, 20), each = 10), times = 2),
                       MoneyB = 10,
                       ShocksA = 10,
                       ShocksB = rep(rep(seq(0, 9), times = 10), times = 2))

utility = function(Payout, Harm, kappa){
  return((Payout * kappa) - (Harm * (1 - kappa)))
}

probability = function(beta, epsilon, gamma, utilitydiff, isLeft){
  prob = 1 / (1 + exp(-(beta * utilitydiff)))
  prob = prob * (1 - 2 * epsilon) + epsilon 
  if (isLeft) {
    prob = prob + gamma * (2 * epsilon)
  } else {
    prob = prob - gamma * (2 * epsilon)
  }
  return(prob)
}

freeParameters = data.frame(kappa = rep(seq(0, 0.95, 0.05), times = 60) + sample(seq(0, 0.05, 0.001), 20*6*10, replace = T), #ranging from 0 to 1, the inverse of kappa has the same range as kappa
                            beta = sample(seq(0,100), 20*6*10, replace = T), # stochasticity
                            epsilon = rep(rep(seq(0, 0.5, 0.1), times = 10), each = 20), # noise
                            gamma = rep(seq(-0.5, 0.4, 0.1), each = 20*6) + sample(seq(0, 0.1, 0.001), 20*6*10, replace = T)) # bias
predictions = data.frame()

### Determine Predictions

for (i in 1:length(freeParameters[,1])){
  Kappa = freeParameters[i,1]
  Beta = freeParameters[i,2]
  Epsilon = freeParameters[i,3]
  Gamma = freeParameters[i,4]
  for (k in 1:length(trialList[,1])){
    shocksMore = trialList$ShocksA[k]
    shocksLess = trialList$ShocksB[k]
    moneyMore = trialList$MoneyA[k]
    moneyLess = trialList$MoneyB[k]
    
    # Just Added
    Utility = utility(Harm = harm(shocksMore, shocksLess),
                      Payout = payout(moneyMore, moneyLess),
                      kappa = Kappa)
    predictions[i,k] = round(probability(Beta, Epsilon, Gamma, Utility, trialList$Left[k]))
  }
}

################ 1_6_0

### Objective Functions

obj_function_oneCondition = function(params, df, optimMethod = "MLE") { #we want 2 kappa and 2 beta parameters (1 per condition), but we don't yet have conditions
  Kappa = params[1]
  Beta = params[2]
  Epsilon = params[3]
  Gamma = params[4]
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 6])
  for (k in 1:length(df[,1])){
    isLeft = df[k, 1]
    moneyMore = df[k, 2]
    moneyLess = df[k, 3]
    shocksMore = df[k, 4]
    shocksLess = df[k, 5]
    
    Utility = utility(Harm = harm(shocksMore, shocksLess),
                      Payout = payout(moneyMore, moneyLess),
                      kappa = Kappa)
    ProbHarm[k] = max(min(probability(Beta, Epsilon, Gamma, Utility, isLeft), 0.9999999999), 0.00000000001)
  }
  if (optimMethod == "OLS"){
    return(sum((choices - ProbHarm)**2))
  } else if (optimMethod == "MLE"){
    return(-sum(choices * log(ProbHarm) + (1 - choices) * log(1 - ProbHarm)))
  }
}

### Optimizers

library(pracma)

initial_params = c(0.5, 4, 0.25, 0)
lower_bounds = c(0, 0, 0, -0.5)
upper_bounds = c(1, 100, 0.5, 0.5)
freeParameters$kappaRecovered = 0
freeParameters$betaRecovered = 0
freeParameters$epsilonRecovered = 0
freeParameters$gammaRecovered = 0

optimize = function(df){
  tryCatch({
    fmincon(obj_function_oneCondition, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "MLE")
  }, error = function(e){
    fmincon(obj_function_oneCondition, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "OLS")
  })
}

for (i in 1:length(freeParameters$kappa)) {
  trialList$Predictions = as.numeric(predictions[i,])
  result = optimize(df = trialList)
  
  freeParameters$kappaRecovered[i] = result$par[1]
  freeParameters$betaRecovered[i] = result$par[2]
  freeParameters$epsilonRecovered[i] = result$par[3]
  freeParameters$gammaRecovered[i] = result$par[4]
  if (mod(i, 100) == 0){message(round(100* (i/(length(freeParameters$kappa)))), '% there', sep = '')}
}

### Verifying Free Parameter Recovery Process

library(ggplot2)

freeParameters$Epsilon = factor(freeParameters$epsilon)

#how well are we able to recover parameters?
qplot(data = freeParameters, x = kappa, y = kappaRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline()
qplot(data = freeParameters, x = beta, y = betaRecovered, color = Epsilon) + geom_smooth(se = F)+ geom_abline()
qplot(data = freeParameters, x = epsilon, y = epsilonRecovered) + geom_smooth(method = 'lm') + geom_abline()
qplot(data = freeParameters, x = gamma, y = gammaRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline() 

