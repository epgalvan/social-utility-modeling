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
    predictions[i,k] = as.numeric(round(probability(Beta, Epsilon, Gamma, Utility, trialList$Left[k]) * 1000)> sample(1000,1))
  }
}

################ 1_6_0

### Objective Functions

obj_function_noConditions = function(params, df, optimMethod = "MLE") { #we want 2 kappa and 2 beta parameters (1 per condition), but we don't yet have conditions
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
    fmincon(obj_function_noConditions, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "MLE")
  }, error = function(e){
    fmincon(obj_function_noConditions, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "OLS")
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

################ 2_1_0

### Preallocating and Defining Functions

library(R.matlab) #data are stored in .mat files so we'll need this to work while in R

parentfolder = 'C:/Users/DELL/Downloads/fMRI_choice_data/run' #the parentfolder where the subject folder/file is
included_subjects = list.files(paste(parentfolder, '1', sep = ''), full.names = F) #get subject files in run 1
included_subjects = gsub('.mat', '', gsub('BC_1', '', included_subjects)) #remove run 1 identifier and .mat suffix

subjectData = data.frame()

grab_data = function(subject){
  if (sum(grepl(pattern = subject, x = list.files(paste(parentfolder, '2/', sep = ''), full.names = F))) == 0){
    return("skipping")
  }
  for (r in 1:2){
    datafile = paste(parentfolder, r, '/BC_', r, subject, '.mat', sep = '')
    temp = readMat(datafile)
    if (r == 1){
      choices = temp["choice"]
      df = data.frame(temp["trialinfo"])
      df = cbind(df[1:76, ], as.data.frame(choices)[1:76, ])
    } else {
      choices = temp["choice"]
      temp = data.frame(temp["trialinfo"])
      temp = cbind(temp[1:76, ], as.data.frame(choices)[1:76, ])
      df = rbind(df, temp) #some have duplicate values
    }
  }
  df[which(is.na(df[,11])), 11] = sample(c(1,2), size = sum(is.na(df[,11])), replace = T)
  df[, 10] = (df[, 10] * -1) + 2 # 1 is harm 0 is help now
  df[, 11] = (df[, 11] * -1) + 2 # 1 is self 0 is other now
  colnames(df) = c("trialNumber", "moneyLess", "moneyDifference", "moneyMore", "shocksLess", "shocksDifference", "shocksMore", "isLeft", "runNumber", "isSelf", "choseHarm")
  
  return(df)
}

generate_predictions = function(df, result){
  df$ProbHarm = 0
  Epsilon = result$par[5]
  Gamma = result$par[6]
  for (k in 1:length(df[,1])){
    shocksMore = df$shocksMore[k]
    shocksLess = df$shocksLess[k]
    moneyMore = df$moneyMore[k]
    moneyLess = df$moneyLess[k]
    isLeft = as.logical(df$isLeft[k])
    isSelf = as.logical(df$isSelf[k])
    if (isSelf) {Kappa = result$par[1]; Beta = result$par[3]} else {Kappa = result$par[2]; Beta = result$par[4]}
    
    # Just Added
    Utility = utility(Harm = harm(shocksMore, shocksLess),
                      Payout = payout(moneyMore, moneyLess),
                      kappa = Kappa)
    df$ProbHarm[k] = max(min(probability(Beta, Epsilon, Gamma, Utility, isLeft), 0.9999999999), 0.00000000001)
  }
  return(df)
}

obj_function = function(params, df, optimMethod = "MLE") { 
  #we want kappa and beta per condition, so we'll assign those based on the condition
  #we need to read the codebook to change the columns to reflect how the data is stored in the .mat files
  Epsilon = params[5]
  Gamma = params[6]
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 11])
  for (k in 1:length(df[,1])){
    isLeft = as.logical(df[k, 8])
    isSelf = as.logical(df[k, 10])
    if (isSelf){Kappa = params[1]; Beta = params[3]} else {Kappa = params[2]; Beta = params[4]}
    moneyMore = as.numeric(df[k, 4])
    moneyLess = as.numeric(df[k, 2])
    shocksMore = as.numeric(df[k, 7])
    shocksLess = as.numeric(df[k, 5])
    
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

#also need to include the correct number of parameter constraints
initial_params = c(0.5, 0.5, 4, 4, 0.25, 0)
lower_bounds = c(0, 0, 0, 0, 0, -0.5)
upper_bounds = c(1, 1, 100, 100, 0.5, 0.5)

#and overcome some errors that seem to be limited to using Log-Odds with fmincon
optimize = function(df){
  tryCatch({
    fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "MLE")
  }, error = function(e){
    fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "OLS")
  })
}

### Recover Free Parameters and Define Predicted Decisions for these Free Parameters

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  if (is.character(df)){next}
  
  result = optimize(df)
  
  # Just Added
  
  df = generate_predictions(df, result)
  
  model_SS = sum((df$choseHarm - df$ProbHarm)**2)
  model_NLL = -2*sum(df$choseHarm * log(df$ProbHarm) + (1 - df$choseHarm) * log(1 - df$ProbHarm))
  
  subjectData[i, 1:9] = c(included_subjects[i], result$par[1], result$par[2], result$par[3], result$par[4], result$par[5], result$par[6],
                     model_SS, model_NLL)
  
  df = cbind(data.frame(SubjectID = rep(included_subjects[i], length(df$trialNumber))), df)
  if (i == 1) {trialData = df} else {trialData = rbind(trialData, df)}
}

colnames(subjectData) = c("subjectID", "Kappa_Self", "Kappa_Other", "Beta_Self", "Beta_Other", "Epsilon", "Gamma", "SS", "Deviance")

################ 2_2_0

subjectData$BIC = as.numeric(subjectData$Deviance) + log(152) * 6

################ 2_3_0

### Preallocating and Defining Functions

#adjusting
optimize = function(of, using, df){
  tryCatch({
    fmincon(of, x0 = initial_params[using], lb = lower_bounds[using], ub = upper_bounds[using], df = df, optimMethod = "MLE")
  }, error = function(e){
    fmincon(of, x0 = initial_params[using], lb = lower_bounds[using], ub = upper_bounds[using], df = df, optimMethod = "OLS")
  })
}

obj_functionNoGamma = function(params, df, optimMethod = "MLE") { 
  #we want kappa and beta per condition, so we'll assign those based on the condition
  #we need to read the codebook to change the columns to reflect how the data is stored in the .mat files
  Epsilon = params[5]
  Gamma = 0 #0 means no bias, everything else is the same
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 11])
  for (k in 1:length(df[,1])){
    isLeft = as.logical(df[k, 8])
    isSelf = as.logical(df[k, 10])
    if (isSelf){Kappa = params[1]; Beta = params[3]} else {Kappa = params[2]; Beta = params[4]}
    moneyMore = as.numeric(df[k, 4])
    moneyLess = as.numeric(df[k, 2])
    shocksMore = as.numeric(df[k, 7])
    shocksLess = as.numeric(df[k, 5])
    
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
obj_functionNoEpsilon = function(params, df, optimMethod = "MLE") { 
  #we want kappa and beta per condition, so we'll assign those based on the condition
  #we need to read the codebook to change the columns to reflect how the data is stored in the .mat files
  Epsilon = 0 #0 means no noise
  Gamma = 0 #and 0 means no bias, everything else is the same
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 11])
  for (k in 1:length(df[,1])){
    isLeft = as.logical(df[k, 8])
    isSelf = as.logical(df[k, 10])
    if (isSelf){Kappa = params[1]; Beta = params[3]} else {Kappa = params[2]; Beta = params[4]}
    moneyMore = as.numeric(df[k, 4])
    moneyLess = as.numeric(df[k, 2])
    shocksMore = as.numeric(df[k, 7])
    shocksLess = as.numeric(df[k, 5])
    
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
obj_functionNoConditions = function(params, df, optimMethod = "MLE") { 
  #we want kappa and beta per condition, so we'll assign those based on the condition
  #we need to read the codebook to change the columns to reflect how the data is stored in the .mat files
  Kappa = params[1] #moved up here because not changing depending on trial
  Beta = params[2] #same as above
  Epsilon = params[3] #change from 5 since 2 fewer parameters
  Gamma = params[4] #change from 6 for same reason
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 11])
  for (k in 1:length(df[,1])){
    isLeft = as.logical(df[k, 8])
    isSelf = as.logical(df[k, 10])
    moneyMore = as.numeric(df[k, 4])
    moneyLess = as.numeric(df[k, 2])
    shocksMore = as.numeric(df[k, 7])
    shocksLess = as.numeric(df[k, 5])
    
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
obj_functionNoKappa = function(params, df, optimMethod = "MLE") { 
  #we want kappa and beta per condition, so we'll assign those based on the condition
  #we need to read the codebook to change the columns to reflect how the data is stored in the .mat files
  Kappa = 0 #nonfactor
  Beta = 0 #nonfactor
  Epsilon = 0.5 #the only factor is noise and bias
  Gamma = params[1] #bias can vary
  
  ProbHarm = vector('numeric', length(df[,1]))
  choices = as.numeric(df[, 11])
  for (k in 1:length(df[,1])){
    isLeft = as.logical(df[k, 8])
    isSelf = as.logical(df[k, 10])
    moneyMore = as.numeric(df[k, 4])
    moneyLess = as.numeric(df[k, 2])
    shocksMore = as.numeric(df[k, 7])
    shocksLess = as.numeric(df[k, 5])
    
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
NoGamma = seq(1, 5)
NoEpsilon = seq(1, 4)
NoConditions = c(1, 3, 5, 6)
NoKappa = 6

altSubjectData = data.frame()

### Recover Free Parameters and Determine Predicted Decisions 

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  if (is.character(df)){next}
  
  resultNoGamma = optimize(obj_functionNoGamma, NoGamma, df)
  resultNoEpsilon = optimize(obj_functionNoEpsilon, NoEpsilon, df)
  resultNoConditions = optimize(obj_functionNoConditions, NoConditions, df)
  resultNoKappa = optim(obj_functionNoKappa,par = initial_params[NoKappa], lower = lower_bounds[NoKappa], upper = upper_bounds[NoKappa], df = df, method = 'L-BFGS-B')
  
  dfNoGamma = generate_predictions(df, data.frame(par = c(resultNoGamma$par[1:5], 0)))
  dfNoEpsilon = generate_predictions(df, data.frame(par = c(resultNoEpsilon$par[1:4], 0, 0)))
  dfNoConditions = generate_predictions(df, data.frame(par = c(resultNoConditions$par[c(1, 1, 2, 2, 3, 4)])))
  dfNoKappa = generate_predictions(df, data.frame(par = c(0, 0, 0, 0, 0.5, resultNoKappa$par[1])))
  
  modelNG_SS = sum((dfNoGamma$choseHarm - dfNoGamma$ProbHarm)**2)
  modelNG_NLL = -2*sum(dfNoGamma$choseHarm * log(dfNoGamma$ProbHarm) + (1 - dfNoGamma$choseHarm) * log(1 - dfNoGamma$ProbHarm))
  modelNE_SS = sum((dfNoEpsilon$choseHarm - dfNoEpsilon$ProbHarm)**2)
  modelNE_NLL = -2*sum(dfNoEpsilon$choseHarm * log(dfNoEpsilon$ProbHarm) + (1 - dfNoEpsilon$choseHarm) * log(1 - dfNoEpsilon$ProbHarm))
  modelNC_SS = sum((dfNoConditions$choseHarm - dfNoConditions$ProbHarm)**2)
  modelNC_NLL = -2*sum(dfNoConditions$choseHarm * log(dfNoConditions$ProbHarm) + (1 - dfNoConditions$choseHarm) * log(1 - dfNoConditions$ProbHarm))
  modelNK_SS = sum((dfNoKappa$choseHarm - dfNoKappa$ProbHarm)**2)
  modelNK_NLL = -2*sum(dfNoKappa$choseHarm * log(dfNoKappa$ProbHarm) + (1 - dfNoKappa$choseHarm) * log(1 - dfNoKappa$ProbHarm))
  
  altSubjectData[i, 1:23] = c(included_subjects[i], 
                              resultNoGamma$par[1], resultNoGamma$par[2], resultNoGamma$par[3], resultNoGamma$par[4], resultNoGamma$par[5], 
                              resultNoEpsilon$par[1], resultNoEpsilon$par[2], resultNoEpsilon$par[3], resultNoEpsilon$par[4], 
                              resultNoConditions$par[1], resultNoConditions$par[2], resultNoConditions$par[3], resultNoConditions$par[4], 
                              resultNoKappa$par[1],
                              modelNG_SS, modelNG_NLL, modelNE_SS, modelNE_NLL, modelNC_SS, modelNC_NLL, modelNK_SS, modelNK_NLL)
  
  
  df = cbind(data.frame(SubjectID = rep(included_subjects[i], length(df$trialNumber))), df)
  df$ProbHarmNG = dfNoGamma$ProbHarm
  df$ProbHarmNE = dfNoEpsilon$ProbHarm
  df$ProbHarmNC = dfNoConditions$ProbHarm
  df$ProbHarmNK = dfNoKappa$ProbHarm
  if (i == 1) {altTrialData = df} else {altTrialData = rbind(altTrialData, df)}
}

colnames(altSubjectData) = c('SubjectID', 'KappaSelfM2', 'KappaOtherM2', 'BetaSelfM2', 'BetaOtherM2', 'EpsilonM2',
                             'KappaSelfM3', 'KappaOtherM3', 'BetaSelfM3', 'BetaOtherM3', 
                             'KappaM4', 'BetaM4', 'EpsilonM4', 'GammaM4', 
                             'GammaM5', 
                             'SSM2', 'DevianceM2', 'SSM3', 'DevianceM3', 'SSM4', 'DevianceM4', 'SSM5', 'DevianceM5')

for (i in 2:ncol(altSubjectData)){altSubjectData[,i] = as.numeric(altSubjectData[,i])}

### Now Compute BIC

altSubjectData$BICM4 = as.numeric(altSubjectData$DevianceM4) + log(152) * 4

### Now Compare BIC

modelBIC = c(sum(subjectData$BIC[-20]), sum(altSubjectData$BICM2[-20]), sum(altSubjectData$BICM3[-20]), sum(altSubjectData$BICM4[-20]), sum(altSubjectData$BICM5[-20]))
which(modelBIC == min(modelBIC))

################ 2_4_0

### Assessing Model Performance

sum(trialData$choseHarm == round(trialData$ProbHarm))/nrow(trialData)

qplot(subjectData$BIC[-20], geom = 'density')
worstExplained = which(subjectData$BIC[-20] > as.numeric(summary(subjectData$BIC[-20])[4]))

qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = moneyDifference, y = choseHarm - ProbHarm, group = SubjectID) + geom_smooth()
qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = shocksDifference, y = choseHarm - ProbHarm, group = SubjectID) + geom_smooth()

### Checking Assumptions

ggplot(data = trialData) + geom_smooth(aes(x = ProbHarm, y = choseHarm, group = isSelf, color = factor(isSelf))) 
ggplot(data = trialData) + geom_density(aes(x = choseHarm, fill = factor(choseHarm)))
normvals = rnorm(1000, mean = 0, sd = sd(trialData$ProbHarm - trialData$choseHarm))
qplot(x = trialData$ProbHarm - trialData$choseHarm, geom = 'density', bw = sd(trialData$ProbHarm - trialData$choseHarm), color = 'Actual') +
  geom_density(aes(x = normvals, color = 'Predicted'), bw = sd(trialData$ProbHarm - trialData$choseHarm))
qplot(x = trialData$shocksDifference, y = trialData$ProbHarm - trialData$choseHarm, geom = 'jitter') + geom_smooth()
qplot(x = trialData$moneyDifference, y = trialData$ProbHarm - trialData$choseHarm, geom = 'jitter') + geom_smooth()
qplot(x = trialData$shocksDifference, y = trialData$ProbHarm - trialData$choseHarm, color = factor(trialData$isSelf), geom = 'smooth') + lims(y = c(-1,1))
qplot(x = trialData$moneyDifference, y = trialData$ProbHarm - trialData$choseHarm, color = factor(trialData$isSelf), geom = 'smooth') + lims(y = c(-1,1))

### Assessing Independence

library(lme4)
library(MuMIn)

ris_model = glmer(data = trialData, choseHarm ~ ProbHarm + (1 + ProbHarm | SubjectID), family = "binomial")
r.squaredGLMM(ris_model)

ri_model = glmer(data = trialData, choseHarm ~ ProbHarm + (1 | SubjectID), family = "binomial")
r.squaredGLMM(ri_model)

ric_model = glmer(data = trialData, choseHarm ~ ProbHarm + isSelf + (1 | SubjectID), family = "binomial")
r.squaredGLMM(ric_model)

### Fivefold Validation

fivefold = data.frame() #preallocate for parameters and errors from the fivefold validation to go into
trialData$ProbHarm_ff = 0
optimize = function(df){
  tryCatch({
    fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "MLE")
  }, error = function(e){
    fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "OLS")
  })
}
adj = 1
for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  if (is.character(df)){adj = adj + 1; next}
  df$ProbHarm = 0
  
  order = sample(152)
  KSelf_ff = vector('numeric', length = 5)
  KOther_ff = vector('numeric', length = 5)
  BSelf_ff = vector('numeric', length = 5)
  BOther_ff = vector('numeric', length = 5)
  E_ff = vector('numeric', length = 5)
  G_ff = vector('numeric', length = 5)
  Deviance_ff = 0
  df$Pred = 0
  for (z in 1:5){
    j = round((z - 1) * (152/5) + 1)
    n = round(z * (152/5))
    withheld = order[j:n]
    
    result_ff = optimize(df[-withheld,])
    
    KSelf_ff[z] = result_ff$par[1]
    KOther_ff[z] = result_ff$par[2]
    BSelf_ff[z] = result_ff$par[3]
    BOther_ff[z] = result_ff$par[4]
    E_ff[z] = result_ff$par[5]
    G_ff[z] = result_ff$par[6]
    df[withheld, ] = generate_predictions(df[withheld, ], result_ff)
  }
  Deviance_ff = -2*sum(df$choseHarm * log(df$ProbHarm) + (1 - df$choseHarm) * log(1 - df$ProbHarm))
  fivefold[i, 1:32] = c(included_subjects[i], Deviance_ff, KSelf_ff, KOther_ff, BSelf_ff, BOther_ff, E_ff, G_ff)
  j = (i-adj)*152 + 1; n = (i-adj + 1)*152
  trialData$ProbHarm_ff[j:n] = df$ProbHarm
}
colnames(fivefold) = c('SubjectID', 'Deviance', 'KS_F1', 'KS_F2', 'KS_F3', 'KS_F4', 'KS_F5',
                       'KO_F1','KO_F2','KO_F3','KO_F4','KO_F5', 'BS_F1', 'BS_F2', 'BS_F3', 'BS_F4', 'BS_F5',
                       'BO_F1','BO_F2','BO_F3','BO_F4','BO_F5', 'E_F1', 'E_F2', 'E_F3', 'E_F4', 'E_F5', 
                       'G_F1', 'G_F2', 'G_F3', 'G_F4', 'G_F5')

sum(round(trialData$ProbHarm_ff) == trialData$choseHarm)/nrow(trialData)
fivefold$BIC = as.numeric(fivefold$Deviance) + log(152) * 6
t.test(fivefold$BIC, subjectData$BIC, paired = T) #test fivefold MFI against normal MFI for this model

library(lsa)
cosines = vector('numeric', length = 60)
for (i in 1:5){
  cosines[i] = cosine(as.numeric(subjectData$Kappa_Self[-20]), as.numeric(fivefold[-20, (i + 2)]))
  cosines[(i+5)] = cosine(as.numeric(subjectData$Kappa_Other[-20]), as.numeric(fivefold[-20, (i + 7)]))
  cosines[(i+10)] = cosine(as.numeric(subjectData$Beta_Self[-20]), as.numeric(fivefold[-20, (i + 12)]))
  cosines[(i+15)] = cosine(as.numeric(subjectData$Beta_Other[-20]), as.numeric(fivefold[-20, (i + 17)]))
  cosines[(i+20)] = cosine(as.numeric(subjectData$Epsilon[-20]), as.numeric(fivefold[-20, (i + 22)]))
  cosines[(i+25)] = cosine(as.numeric(subjectData$Gamma[-20]), as.numeric(fivefold[-20, (i + 27)]))
}

mean(cosines[1:5]) 
mean(cosines[6:10])
mean(cosines[11:15]) 
mean(cosines[16:20]) 
mean(cosines[21:25]) 
mean(cosines[26:30])  

################ 3_1_0

t.test(subjectData$BIC[-20], altSubjectData$BICM2[-20], paired = T)
t.test(subjectData$BIC[-20], altSubjectData$BICM3[-20], paired = T)
t.test(subjectData$BIC[-20], altSubjectData$BICM4[-20], paired = T)
t.test(subjectData$BIC[-20], altSubjectData$BICM5[-20], paired = T)

################ 3_2_0

t.test(as.numeric(subjectData$Kappa_Self[-20]), as.numeric(subjectData$Kappa_Other[-20]), paired = T)
mean(as.numeric(subjectData$Kappa_Self[-20]))
mean(as.numeric(subjectData$Kappa_Other[-20]))

################ 3_3_0

resultNID = fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = trialData, optimMethod = "OLS")
dfNID = generate_predictions(trialData, resultNID)
sum(dfNID$choseHarm == round(dfNID$ProbHarm))/nrow(dfNID)

altSubjectData$Deviance_NID = 0

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i]); if (is.character(df)){next} else {trials = which(included_subjects[i] == dfNID$SubjectID)}
  df = dfNID[trials, ]
  altSubjectData$Deviance_NID[i] = -2*sum(df$choseHarm * log(df$ProbHarm) + (1 - df$choseHarm) * log(1 - df$ProbHarm))
}

altSubjectData$BIC_NID = as.numeric(altSubjectData$Deviance_NID) + log(152) * 6/(nrow(dfNID))

t.test(subjectData$BIC[-20], altSubjectData$BIC_NID[-20], paired = T)
which(modelBIC > sum(altSubjectData$BIC_NID))
