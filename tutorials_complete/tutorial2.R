################ 1_4_0

### Create Construct Value Functions

new_value = function(initial_allocation, tax_rate_decimal, number_tokens_game = 100, number_players_game = 10){
  return(round(initial_allocation - (tax_rate_decimal * initial_allocation) + ((number_tokens_game * tax_rate_decimal)/(number_players_game))))
}

payout_maximization = function(chosen_outcome_self, possible_outcomes_self){
  return(chosen_outcome_self/max(possible_outcomes_self))
}

equality = function(chosen_outcomes_all, intial_allocations_all, perfect_equality = 100/10){
  return((1 - sum((chosen_outcomes_all - perfect_equality)**2)/sum((intial_allocations_all - perfect_equality)**2)))
}

equity = function(chosen_outcomes_all, intial_allocations_all, perfect_inequity = 100/10){
  return((1 - sum((chosen_outcomes_all - intial_allocations_all)**2)/sum((perfect_equality - intial_allocations_all)**2)))
}

################ 1_5_0

### Preallocating and Defining Functions, TrialList, and Parameters

#first, create a noisy resource distribution that has gini between 0.3 and 0.4 where the maximum benefit or loss is approximately going to be 10 tokens
shares = seq(0.05,0.95, 0.1)**1.25
df = data.frame()
for (k in 1:10){
  df[1:20, k] = rnorm(20, mean=shares[k], sd=0.01*sum(shares))
  df[which(df[,k] < 0),] = 0
}

#second, convert to a rounded percent
for (k in 1:length(df[,1])) {
  df[k,1:10] = round((df[k,1:10]/sum(df[k,1:10]))*100)
}

#third, ensure that there are exactly 100 tokens on each trial
for (k in 1:length(df[,1])) {
  if (sum(df[k,1:10]) < 100){
    for (j in 1:(100-sum(df[k,1:10]))){
      i = sample(1:10, 1)
      df[k, i] = df[k, i] + 1
    }
  }
  if (sum(df[k,1:10]) > 100){
    for (j in 1:(sum(df[k,1:10]) - 100)){
      i = sample(which(df[k,1:10] > 0), 1)
      df[k, i] = df[k, i] - 1
    }
  }
}

trialList = data.frame()

#fourth, ensure that our subject is in each decile the same number of times
for (k in 1:length(df[,1])){
  i = round((k/2)+0.05) #because this function rounds down on even numbers and up on odd numbers
  trialList[k, 1] = df[k, i]
  intermediate = df[k, -i]
  trialList[k, 2:10] = intermediate[1,sample(9)] #to insure that other players on screen are not all
}

#trialList above
choices = seq(0, 1, 0.1) #tax rate

utility = function(theta, phi, Equity, Equality, Payout){
  return((theta * Payout) + ((1 - theta) * ((phi * Equality) + ((1 - phi) * Equity))))
}

freeParameters = data.frame(theta = seq(0, 1, 0.01),
                            phi = seq(0, 1, 0.01))

predictions = data.frame()

### Determine Predictions

non_specific = rep(0, length(freeParameters[,1])) # Just Added This Line

for (i in 1:length(freeParameters[,1])){
  Theta = freeParameters[i,1]
  Phi = freeParameters[i,2]
  for (k in 1:length(trialList[,1])){
    
    # Just Added
    Utility = vector('numeric', length(choices))
    for (n in 1:length(choices)){
      Utility[n] = utility(theta = Theta,
                           phi = Phi,
                           Equity = equity(new_value(trialList[k, 1:10], choices[n]), trialList[k, 1:10], choices[n]),
                           Equality = equality(new_value(trialList[k, 1:10], choices[n]), trialList[k, 1:10], choices[n]),
                           Payout = payout(new_value(trialList[k, 1], choices[n]), trialList[k, 1], choices[n]))
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

obj_function = function(params, df, method = "OLS") {
  Theta = params[1]
  Phi = params[2]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  choices = seq(0, 1, 0.1)
  
  for (k in 1:length(df[,1])){
    Utility = vector('numeric', length(choices))
    for (n in 1:length(choices)){
      Utility[n] = utility(theta = Theta,
                           phi = Phi,
                           Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                           Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                           Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[which(df[k, 11] == choices)]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
}

### Optimizers

library(pracma)

initial_params = c(0.5, 0.5)
lower_bounds = c(0, 0)
upper_bounds = c(1, 1)
theta_recovered = vector('numeric', 11**2)
phi_recovered = vector('numeric', 11**2)
theta_true = rep(seq(0, 1, 0.1), length(seq(0, 1, 0.1)))
phi_true = rep(seq(0, 1, 0.1), length(seq(0, 1, 0.1)))

for (i in 1:length(theta_true)) {
  this_idx = which(round(theta_true[i], 2) == round(freeParameters$theta, 2) & round(phi_true[i], 2) == round(freeParameters$phi, 2))
  trialList_this = trialList
  trialList$Decision = as.numeric(predictions[this_idx, 1:length(trialList[1,])])
  result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   df = trialList_this)
  
  parameter1_recovered[i] = result$par[1]
  parameter2_recovered[i] = result$par[2]
}

### Verifying Free Parameter Recovery Process

library(ggplot2)

distance = (theta_recovered - theta_true)**2 + (phi_recovered - phi_true)**2
qplot(x = theta_true, y = phi_true, color = distance, size = distance, geom = 'point') + scale_radius(limits=c(0, sqrt(2)), range=c(0, 20))

distance_new = (theta_recovered - theta_true)**2 + (((phi_recovered - 0.5) * theta_recovered) - ((phi_true - 0.5) * theta_true))**2
qplot(x = theta_true, y = phi_true, color = distance_new, size = distance_new, geom = 'point') + scale_radius(limits=c(0, sqrt(2)), range=c(0, 20))

################ 1_7_0

### Cluster Your Data Using HAC

distance_mat = dist(predictions, method = 'euclidean')
set.seed(240)
hierarchical = hclust(distance_mat, method = 'average')
plot(hierarchical)
fit = cutree(hierarchical, k = 4)

### Identify Where Your Clusters Are

freeParameters$Strategy = as.character(fit)
freeParameters$Strategy[which(freeParameters$Strategy == freeParameters$Strategy[1])] = 'Equity-Seeking'
freeParameters$Strategy[which(freeParameters$Strategy == freeParameters$Strategy[10101])] = 'Equality-Seeking'
freeParameters$Strategy[which(freeParameters$Strategy != 'Equity-Seeking' & freeParameters$Strategy != 'Equality-Seeking' & freeParameters$Strategy != 'Guilt-Averse')] = 'Payout-Maximizers';
freeParameters$Strategy = as.factor(freeParameters$Strategy) #Strategy clusters
model_space = ggplot(data = freeParameters, aes(x = theta, y = phi, color = Strategy)) +
  labs(x = 'Theta', y = 'Phi', color = 'Strategy') + geom_point(size = 2.5) +
  scale_color_manual(values = c(rgb(0,130,229, maxColorValue = 255), rgb(255,25,0, maxColorValue = 255), rgb(174,0,255, maxColorValue = 255))); model_space

### Examine Model Predictions Efficiently

toPlot = data.frame()
for (i in 1:length(freeParameters[,1])){
  replacement = ((i - 1) * 20 + 1):(i * 20)
  toPlot[replacement, 1] = freeParameters$Strategy[i]
  toPlot[replacement, 2] = as.numeric(trialList[, 1])
  toPlot[replacement, 3] = as.numeric(predictions[i, 1:20])
  toPlot[replacement, 4] = new_value(trialList[, 1], as.numeric(predictions[i, 1:20])) - trialList[, 1]
}
colnames(toPlot) = c('Strategy', 'Initial_Allocation', 'Tax_Rate', 'Payout_Change')

ggplot(data = toPlot, aes(x = Initial_Allocation, y = Payout_Change, group = Strategy, color = Strategy)) + geom_smooth(se = TRUE) +
  scale_color_manual(values = c(rgb(0,130,229, maxColorValue = 255), rgb(255,25,0, maxColorValue = 255), rgb(174,0,255, maxColorValue = 255))) +
  lims(x = c(0, 10), y = c(0, 30)) + labs(x = 'Initial Allocation', 'Payout Change')

ggplot(data = toPlot, aes(x = Initial_Allocation, y = Tax_Rate, group = Strategy, color = Strategy)) + geom_smooth(se = TRUE) +
  scale_color_manual(values = c(rgb(0,130,229, maxColorValue = 255), rgb(255,25,0, maxColorValue = 255), rgb(174,0,255, maxColorValue = 255))) +
  lims(x = c(0, 10), y = c(0, 30)) + labs(x = 'Initial Allocation', 'Tax Rate')

################ 2_1_0

### Preallocating and Defining Functions

all_subjects = read.csv2('C:/Users/DELL/Downloads/prolific_export_648c19e5420c9b10a79589a4.csv', sep = ',')
potentially_excluded_subjects = read.csv2('C:/Users/DELL/Downloads/prolific_export_6437e471bec005411b1503ea.csv', sep = ',')

k = 0
for (i in 1:length(potentially_excluded_subjects$Participant.id)){
  k = k + sum(potentially_excluded_subjects$Participant.id[i] == all_subjects$Participant.id)
}
actual_excluded_subjects = vector('character', k)
j = 0
for (i in 1:length(potentially_excluded_subjects$Participant.id)){
  if (sum(potentially_excluded_subjects$Participant.id[i] == all_subjects$Participant.id) == 1) {
    j = j + 1
    actual_excluded_subjects[j] = which(potentially_excluded_subjects$Participant.id[i] == all_subjects)
  }
}
included_subjects = all_subjects[-actual_excluded_subjects]

parentfolder = 'C:/Users/DELL/rdb_all_modded/' #the parentfolder where the subject folder/file is
restoffilepath = '.csv' #everything after the subject folder/file name

trialData = data.frame()
subjectData = data.frame()

conditions = c('merit', 'entitlement', 'corruption', 'luck')

### Recover Free Parameters and Define Predicted Decisions for these Free Parameters

for (i in 1:length(included_subjects)){
  datafile = paste(parentfolder, included_subjects[i], restoffilepath, sep = '') # produces a character vector 'parentfolder/included_subjects[i]restoffilepath'
  fullData = read.csv2(datafile)
  
  thetaPerCondition = vector('numeric', length(conditions))
  phiPerCondition = vector('numeric', length(conditions))
  strategyPerCondition = vector('numeric', length(conditions))
  SSPerCondition = vector('numeric', length(conditions))
  
  for (j in 1:length(conditions)){
    
    df = fullData[which(fullData$condition == conditions[j]), c(49, 40:48, 33)] #49 is subject's initial allocation, 40:48 are players 1:9 initial allocation, 33 is redistribution rate
    df$redistributionRate = df$redistributionRate/100 #converting to a decimal from a percent
    
    result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                     lb = lower_bounds, ub = upper_bounds,
                     df = df)
    thetaPerCondition[j] = result$par[1]
    phiPerCondition[j] = result$par[2]
    closestPoint = which(as.numeric(freeParameters[,1]) == round(result$par[1], 2) & as.numeric(freeParameters[,2]) == round(result$par[2], 2))
    strategyPerCondition[j] = freeParameters$Strategy[closestPoint]
    
    #Just Added
    
    df$predictedRR = vector('numeric')
    df$predictedOutcome = vector('numeric')
    df$actualOutcome = vector('numeric')
    
    for (k in 1:length(df$redistributionRate)){
      Utility = vector('numeric', length(Choices))
      for (n in 1:length(Choices)){
        Utility[n] = utility(theta = thetaPerCondition[j],
                             phi = phiPerCondition[j],
                             Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                             Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                             Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
      }
      correct_choice = which(Utility == max(Utility))
      df$predictedRR[k] = Choices[correct_choice[sample(length(correct_choice), 1)]]
      df$predictedOutcome[k] = new_value(df$myself[k], df$predictedRR[k])
      df$actualOutcome[k] = new_value(df$myself[k], df$redistributionRate[k])
    }
    SSPerCondition[j] = sum((df$predictedOutcome - df$actualOutcome)**2)
    fullData[which(fullData$condition == conditions[j]), 76] = df$predictedOutcome
    fullData[which(fullData$condition == conditions[j]), 77] = df$predictedRR
    fullData[which(fullData$condition == conditions[j]), 78] = df$actualOutcome
    fullData[which(fullData$condition == conditions[j]), 79] = freeParameters$Strategy[closestPoint]
  }
  
  subjectData[i, 1:17] = c(included_subjects[i], thetaPerCondition, phiPerCondition, strategyPerCondition, SSPerCondition)
  
  start = length(subjectData[, 1]) + 1
  end = start + length(fullData$redistributionRate)
  trialData[start:end, 1] = included_subjects[i]
  trialData[start:end, 2] = fullData$redistributionRate/100
  trialData[start:end, 3] = fullData[,78]
  trialData[start:end, 4] = fullData[,77]
  trialData[start:end, 5] = fullData[,76]
  trialData[start:end, 6] = fullData$myself
  trialData[start:end, 7] = fullData[, 79]
  trialData[start:end, 9] = rep(seq(1, 4), each = length(df$redistributionRate))
  trialData[start:end, 10] = seq(1, length(fullData$redistributionRate))
  trialData[start:end, 11] = fullData$condition
  trialData[start:end, 12:20] = fullData[, c40:48]
}
colnames(subjectData) = c('SubjectID', 'thetaMerit', 'thetaEntitlement', 'thetaCorruption', 'thetaLuck', 'phiMerit', 'phiEntitlement', 'phiCorruption', 'phiLuck', 'strategyMerit', 'strategyEntitlement', 'strategyCorruption', 'strategyLuck', 'SSMerit', 'SSEntitlement', 'SSCorruption', 'SSLuck')
colnames(trialData) = c('SubjectID', 'observedTaxRate', 'observedOutcome', 'predictedTaxRate', 'predictedOutcome', 'initialAllocation', 'strategy', 'blockNumber', 'trialNumber', 'condition',
                        'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9')

################ 2_2_0

N = length(df[, 1])
k = 2
subjectData$AICMerit = N * log(subjectData$SSMerit/N) + 2*k
subjectData$AICEntitlement = N * log(subjectData$SSEntitlement/N) + 2*k
subjectData$AICCorruption = N * log(subjectData$SSCorruption/N) + 2*k
subjectData$AICLuck = N * log(subjectData$SSLuck/N) + 2*k

################ 2_3_0

### Create Utility Equations for Our Alternative Models

utilityP = function(Payout){
  return(Payout)
}
utilityL = function(Equality){
  return(Equality)
}
utilityT = function(Equity){
  return(Equity)
}
utilityPL = function(Payout, Equality, theta){
  return(theta * Payout + (1 - theta) * Equality)
}
utilityPT = function(Payout, Equity, theta){
  return(theta * Payout + (1 - theta) * Equity)
}
utilityLT = function(Equity, Equality, phi){
  return(phi * Equality + (1 - phi) * Equity)
}

### Preallocating and Defining Functions

obj_functionPL = function(params, df, method = "OLS") {
  Theta = params[1]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  choices = seq(0, 1, 0.1)
  
  for (k in 1:length(df[,1])){
    Utility = vector('numeric', length(choices))
    for (n in 1:length(choices)){
      Utility[n] = utilityPL(theta = Theta,
                             Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                             Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[which(df[k, 11] == choices)]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
}

obj_functionPT = function(params, df, method = "OLS") {
  Theta = params[1]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  choices = seq(0, 1, 0.1)
  
  for (k in 1:length(df[,1])){
    Utility = vector('numeric', length(choices))
    for (n in 1:length(choices)){
      Utility[n] = utilityPT(theta = Theta,
                             Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                             Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[which(df[k, 11] == choices)]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
}

obj_functionLT = function(params, df, method = "OLS") {
  Phi = params[1]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  choices = seq(0, 1, 0.1)
  
  for (k in 1:length(df[,1])){
    Utility = vector('numeric', length(choices))
    for (n in 1:length(choices)){
      Utility[n] = utilityLT(phi = Phi,
                             Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                             Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[which(df[k, 11] == choices)]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
}

altSubjectData = data.frame()

### Recover Free Parameters and Determine Predicted Decisions Per Condition

for (i in 1:length(included_subjects)){
  datafile = paste(parentfolder, included_subjects[i], restoffilepath, sep = '') # produces a character vector 'parentfolder/included_subjects[i]**.filetype'
  fullData = read.csv2(datafile)
  
  parametersPerCondition = vector('numeric', length(conditions) * 3)
  SSPerCondition = vector('numeric', length(conditions) * 6)
  
  for (j in 1:length(conditions)){
    df = fullData[which(fullData$condition == conditions[j]), c(49, 40:48), 33] #49 is subject's initial allocation, 40:48 are players 1:9 initial allocation, 33 is redistribution rate
    df$redistributionRate = df$redistributionRate/100 #converting to a decimal from a percent
    
    resultPL = optim(obj_functionPL, par = 0.5, lower = 0, upper = 1, df = df, method = 'L-BFGS-B')
    resultPT = optim(obj_functionPT, par = 0.5, lower = 0, upper = 1, df = df, method = 'L-BFGS-B')
    resultLT = optim(obj_functionLT, par = 0.5, lower = 0, upper = 1, df = df, method = 'L-BFGS-B')
    
    parametersPerCondition[(((j - 1) * 3) + 1):(j * 3)] = c(resultPL$par[1], resultPT$par[1], resultLT$par[1])
    
    #Just Added
    
    df$PredictionP = vector('numeric')
    df$PredictionL = vector('numeric')
    df$PredictionT = vector('numeric')
    df$PredictionPT = vector('numeric')
    df$PredictionPL = vector('numeric')
    df$PredictionLT = vector('numeric')
    for (k in 1:length(df$Decisions)){
      UtilityP = vector('numeric', length(Choices))
      UtilityL = vector('numeric', length(Choices))
      UtilityT = vector('numeric', length(Choices))
      UtilityPT = vector('numeric', length(Choices))
      UtilityPL = vector('numeric', length(Choices))
      UtilityLT = vector('numeric', length(Choices))
      for (n in 1:length(Choices)){
        UtilityP[n] = utilityP(Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
        UtilityL[n] = utilityL(Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]))
        UtilityT[n] = utilityT(Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]))
        UtilityPT[n] = utilityPT(theta = resultPT$par[1],
                                 Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                                 Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
        UtilityPL[n] = utilityPL(theta = resultPL$par[1],
                                 Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                                 Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
        UtilityLT[n] = utilityLT(phi = resultLT$par[1],
                                 Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                                 Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]))
      }
      correct_choiceP = which(UtilityP == max(UtilityP))
      correct_choiceL = which(UtilityL == max(UtilityL))
      correct_choiceT = which(UtilityT == max(UtilityT))
      correct_choicePT = which(UtilityPT == max(UtilityPT))
      correct_choicePL = which(UtilityPL == max(UtilityPL))
      correct_choiceLT = which(UtilityLT == max(UtilityLT))
      
      df$predictedP[k] = new_value(df$myself[k], Choices[correct_choiceP[sample(length(correct_choiceP), 1)]])
      df$predictedL[k] = new_value(df$myself[k], Choices[correct_choiceL[sample(length(correct_choiceL), 1)]])
      df$predictedT[k] = new_value(df$myself[k], Choices[correct_choiceT[sample(length(correct_choiceT), 1)]])
      df$predictedPL[k] = new_value(df$myself[k], Choices[correct_choicePT[sample(length(correct_choicePT), 1)]])
      df$predictedPT[k] = new_value(df$myself[k], Choices[correct_choicePL[sample(length(correct_choiceTL), 1)]])
      df$predictedLT[k] = new_value(df$myself[k], Choices[correct_choiceLT[sample(length(correct_choiceLT), 1)]])
    }
    df$Outcome = new_value(df$myself, df$redistributionRate)
    SSPerCondition[(((j - 1) * 6) + 1):(j*6)] = c(sum((df$Outcome - df$predictedP)**2),
                                                  sum((df$Outcome - df$predictedL)**2),
                                                  sum((df$Outcome - df$predictedT)**2),
                                                  sum((df$Outcome - df$predictedPL)**2),
                                                  sum((df$Outcome - df$predictedPT)**2),
                                                  sum((df$Outcome - df$predictedLT)**2))
  }
  altSubjectData[i, 1:37] = c(included_subjects[i], parametersPerCondition, SSPerCondition)
}
colnames(altSubjectData) = c('SubjectID',
                             'thetaPLMerit', 'thetaPTMerit', 'phiLTMerit',
                             'thetaPLEntitlement', 'thetaPTEntitlement', 'phiLTEntitlement',
                             'thetaPLCorruption', 'thetaPTCorruption', 'phiLTCorruption',
                             'thetaPLLuck', 'thetaPTLuck', 'phiLTLuck',
                             'SSPMerit', 'SSLMerit', 'SSTMerit', 'SSPLMerit', 'SSPTMerit', 'SSLTMerit',
                             'SSPEntitlement', 'SSLEntitlement', 'SSTEntitlement', 'SSPLEntitlement', 'SSPTEntitlement', 'SSLTEntitlement',
                             'SSPCorruption', 'SSLCorruption', 'SSTCorruption', 'SSPLCorruption', 'SSPTCorruption', 'SSLTCorruption',
                             'SSPLuck', 'SSLLuck', 'SSTLuck', 'SSPLLuck', 'SSPLuck', 'SSLTLuck')

### Now Compute BIC Per Model Per Condition

N = length(df[, 1])
k = 2

altSubjectData$AICPMerit = N * log(altSubjectData$SSPMerit/N) + 2*k
altSubjectData$AICLMerit = N * log(altSubjectData$SSLMerit/N) + 2*k
altSubjectData$AICTMerit = N * log(altSubjectData$SSTMerit/N) + 2*k
altSubjectData$AICPLMerit = N * log(altSubjectData$SSPLMerit/N) + 2*k
altSubjectData$AICPTMerit = N * log(altSubjectData$SSPTMerit/N) + 2*k
altSubjectData$AICLTMerit = N * log(altSubjectData$SSLTMerit/N) + 2*k

altSubjectData$AICPEntitlement = N * log(altSubjectData$SSPEntitlement/N) + 2*k
altSubjectData$AICLEntitlement = N * log(altSubjectData$SSLEntitlement/N) + 2*k
altSubjectData$AICTEntitlement = N * log(altSubjectData$SSTEntitlement/N) + 2*k
altSubjectData$AICPLEntitlement = N * log(altSubjectData$SSPLEntitlement/N) + 2*k
altSubjectData$AICPTEntitlement = N * log(altSubjectData$SSPTEntitlement/N) + 2*k
altSubjectData$AICLTEntitlement = N * log(altSubjectData$SSLTEntitlement/N) + 2*k

altSubjectData$AICPCorruption = N * log(altSubjectData$SSPCorruption/N) + 2*k
altSubjectData$AICLCorruption = N * log(altSubjectData$SSLCorruption/N) + 2*k
altSubjectData$AICTCorruption = N * log(altSubjectData$SSTCorruption/N) + 2*k
altSubjectData$AICPLCorruption = N * log(altSubjectData$SSPLCorruption/N) + 2*k
altSubjectData$AICPTCorruption = N * log(altSubjectData$SSPTCorruption/N) + 2*k
altSubjectData$AICLTCorruption = N * log(altSubjectData$SSLTCorruption/N) + 2*k

altSubjectData$AICPLuck = N * log(altSubjectData$SSPLuck/N) + 2*k
altSubjectData$AICLLuck = N * log(altSubjectData$SSLLuck/N) + 2*k
altSubjectData$AICTLuck = N * log(altSubjectData$SSTLuck/N) + 2*k
altSubjectData$AICPLLuck = N * log(altSubjectData$SSPLLuck/N) + 2*k
altSubjectData$AICPTLuck = N * log(altSubjectData$SSPTLuck/N) + 2*k
altSubjectData$AICLTLuck = N * log(altSubjectData$SSLTLuck/N) + 2*k

### Now Compare BIC

excludedM = which(is.infinite(subjectData$AICMerit)) #any perfectly predicted subjects by alternative models should also be perfectly predicted by the most complex model as well
excludedE = which(is.infinite(subjectData$AICEntitlement))
excludedC = which(is.infinite(subjectData$AICCorruption))
excludedL = which(is.infinite(subjectData$AICLuck))

models = c('Favored Model', 'Equality Only', 'Equity Only', 'Payout Only', 'Payout-Equality', 'Payout-Equity', 'Equality-Equity')

bestModelMerit = models[which(c(mean(subjectData$AICMerit[-excluded]),
                                mean(altSubjectData$AICPMerit[-excluded]),
                                mean(altSubjectData$AICLMerit[-excluded]),
                                mean(altSubjectData$AICTMerit[-excluded]),
                                mean(altSubjectData$AICPLMerit[-excluded]),
                                mean(altSubjectData$AICPTMerit[-excluded]),
                                mean(altSubjectData$AICLTMerit[-excluded])) == max(c(mean(subjectData$AICMerit[-excluded]),
                                                                                     mean(altSubjectData$AICPMerit[-excluded]),
                                                                                     mean(altSubjectData$AICLMerit[-excluded]),
                                                                                     mean(altSubjectData$AICTMerit[-excluded]),
                                                                                     mean(altSubjectData$AICPLMerit[-excluded]),
                                                                                     mean(altSubjectData$AICPTMerit[-excluded]),
                                                                                     mean(altSubjectData$AICLTMerit[-excluded]))))]

bestModelEntitlement = models[which(c(mean(subjectData$AICEntitlement[-excluded]),
                                      mean(altSubjectData$AICPEntitlement[-excluded]),
                                      mean(altSubjectData$AICLEntitlement[-excluded]),
                                      mean(altSubjectData$AICTEntitlement[-excluded]),
                                      mean(altSubjectData$AICPLEntitlement[-excluded]),
                                      mean(altSubjectData$AICPTEntitlement[-excluded]),
                                      mean(altSubjectData$AICLTEntitlement[-excluded])) == max(c(mean(subjectData$AICEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICPEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICLEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICTEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICPLEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICPTEntitlement[-excluded]),
                                                                                                 mean(altSubjectData$AICLTEntitlement[-excluded]))))]

bestModelCorruption = models[which(c(mean(subjectData$AICCorruption[-excluded]),
                                     mean(altSubjectData$AICPCorruption[-excluded]),
                                     mean(altSubjectData$AICLCorruption[-excluded]),
                                     mean(altSubjectData$AICTCorruption[-excluded]),
                                     mean(altSubjectData$AICPLCorruption[-excluded]),
                                     mean(altSubjectData$AICPTCorruption[-excluded]),
                                     mean(altSubjectData$AICLTCorruption[-excluded])) == max(c(mean(subjectData$AICCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICPCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICLCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICTCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICPLCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICPTCorruption[-excluded]),
                                                                                               mean(altSubjectData$AICLTCorruption[-excluded]))))]

bestModelLuck = models[which(c(mean(subjectData$AICLuck[-excluded]),
                               mean(altSubjectData$AICPLuck[-excluded]),
                               mean(altSubjectData$AICLLuck[-excluded]),
                               mean(altSubjectData$AICTLuck[-excluded]),
                               mean(altSubjectData$AICPLLuck[-excluded]),
                               mean(altSubjectData$AICPTLuck[-excluded]),
                               mean(altSubjectData$AICLTLuck[-excluded])) == max(c(mean(subjectData$AICLuck[-excluded]),
                                                                                   mean(altSubjectData$AICPLuck[-excluded]),
                                                                                   mean(altSubjectData$AICLLuck[-excluded]),
                                                                                   mean(altSubjectData$AICTLuck[-excluded]),
                                                                                   mean(altSubjectData$AICPLLuck[-excluded]),
                                                                                   mean(altSubjectData$AICPTLuck[-excluded]),
                                                                                   mean(altSubjectData$AICLTLuck[-excluded]))))]

################ 2_4_0

### Assessing Model Performance

### Model Performance

taxRatePredictions = lm(data = trialData, observedTaxRate ~ predictedTaxRate)
summary(taxRatePredictions) # R-squared predicting the subject's tax rate

payoutPredictions = lm(data = trialData, observedOutcome ~ predictedOutcome)
summary(payoutPredictions) # R-squared predicting the subject's payout, we're going to focus on this

### Checking Assumptions

qplot(x = trialData$predictedOutcome, y = trialData$observedOutcome, geom = 'smooth') +
  geom_abline(slope = 1, intercept = 0) + labs(x = 'Prediction', y = 'Observed') + lims(x = c(0, 30), y = c(0, 30)) # linearity

normvals = rnorm(1000, mean = 0, sd = sd(trialData$predictedOutcome - trialData$observedOutcome))
qplot(x = (trialData$predictedOutcome - trialData$observedOutcome), geom = 'density', bw = 1, color = 'Actual') +
  geom_density(aes(x = normvals, color = 'Predicted'), bw = 1) + labs(x = 'Prediction - Observed', y = 'Density') #normality

qplot(x = trialData$initialAllocation, y = (trialData$predictedOutcome - trialData$observedOutcome), group = as.factor(trialData$Multiplier), color = as.factor(trialData$Multiplier), geom = 'smooth')  +
  labs(x = 'Initial Allocation', y = 'Prediction - Observed', color = 'Multiplier') #independence of error and homoscedasticity

### Assessing Independence

library(lme4)
library(MuMin)
ris_model = lmer(data = trialData, observedOutcome ~ predictedOutcome + (1 + predictedOutcome | SubjectID))
summary(ris_model) #model should have issues
r.squaredGLMM(ris_model)

cri_model = lmer(data = trialData, observedOutcome ~ predictedOutcome + condition + (1 | SubjectID))
summary(cri_model) #model should have issues
r.squaredGLMM(cri_model)

cric_model = lmer(data = trialData, observedOutcome ~ predictedOutcome + condition + (1 + condition | SubjectID))
summary(cric_model) #model should have issues
r.squaredGLMM(cric_model)

ri_model = lmer(data = trialData, observedOutcome ~ predictedOutcome + (1 | SubjectID))
summary(ri_model)
r.squaredGLMM(ri_model) #if conditional Rsq is between 0 and 0.05 lower than the multiple Rsq, that's good enough

### Fivefold Validation

fivefold = data.frame()

for (i in 1:length(included_subjects)){
  datafile = paste(parentfolder, included_subjects[i], restoffilepath, sep = '') # produces a character vector 'parentfolder/included_subjects[i]**.filetype'
  fullData = read.csv2(datafile)
  
  thetaPerCondition = vector('numeric', length(conditions) * 5)
  phiPerCondition = vector('numeric', length(conditions) * 5)
  SSPerCondition = vector('numeric', length(conditions))
  
  for (j in 1:length(conditions)){
    
    df = fullData[which(fullData$condition == conditions[j]), c(49, 40:48), 33] #49 is subject's initial allocation, 40:48 are players 1:9 initial allocation, 33 is redistribution rate
    df$redistributionRate = df$redistributionRate/100 #converting to a decimal from a percent
    
    order = sample(20)
    
    for (z in 1:5){
      j = (z - 1) * 4 + 1
      n = z * 4
      withheld = order[j:n]
      m = ((j - 1) * 5) + z
      result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                       lb = lower_bounds, ub = upper_bounds,
                       df = df[-withheld,])
      
      thetaPerCondition[m] = result$par[1]
      phiPerCondition[m] = result$par[2]
      for (k in 1:length(withheld)){
        y = withheld[k]
        Utility = vector('numeric', length(Choices))
        for (n in 1:length(Choices)){
          Utility[n] = utility(theta = thetaPerCondition[m],
                               phi = phiPerCondition[m],
                               Equity = equity(new_value(df[y, 1:10], choices[n]), df[y, 1:10], choices[n]),
                               Equality = equality(new_value(df[y, 1:10], choices[n]), df[y, 1:10], choices[n]),
                               Payout = payout(new_value(df[y, 1], choices[n]), df[y, 1], choices[n]))
        }
        correct_choice = which(Utility == max(Utility))
        predictedOutcome = new_value(df$myself[y], Choices[correct_choice[sample(length(correct_choice), 1)]])
        observedOutcome = new_value(df$myself[y], df$redistributionRate)
        SSPerCondition = SSPerCondition + (observedOutcome - predictedOutcome)**2
      }
    }
  }
  
  fivefold[i, 1:45] = c(included_subjects[i], thetaPerCondition, phiPerCondition, SSPerCondition)
}
colnames(fivefold) = c('SubjectID',
                       'thetaMeritF1', 'thetaMeritF2', 'thetaMeritF3', 'thetaMeritF4', 'thetaMeritF5',
                       'thetaEntitlementF1', 'thetaEntitlementF2', 'thetaEntitlementF3', 'thetaEntitlementF4', 'thetaEntitlementF5',
                       'thetaCorruptionF1', 'thetaCorruptionF2', 'thetaCorruptionF3', 'thetaCorruptionF4', 'thetaCorruptionF5',
                       'thetaLuckF1', 'thetaLuckF2', 'thetaLuckF3', 'thetaLuckF4', 'thetaLuckF5',
                       'phiMeritF1', 'phiMeritF2', 'phiMeritF3', 'phiMeritF4', 'phiMeritF5',
                       'phiEntitlementF1', 'phiEntitlementF2', 'phiEntitlementF3', 'phiEntitlementF4', 'phiEntitlementF5',
                       'phiCorruptionF1', 'phiCorruptionF2', 'phiCorruptionF3', 'phiCorruptionF4', 'phiCorruptionF5',
                       'phiLuckF1', 'phiLuckF2', 'phiLuckF3', 'phiLuckF4', 'phiLuckF5',
                       'SSMerit', 'SSEntitlement', 'SSCorruption', 'SSLuck')

sqrt(mean(fivefold$SSMerit)/length(df$redistributionRate)) - sqrt(mean(subjectData$SSMerit)/length(df$redistributionRate)) #the change in root mean squared error, per trial in merit condition
sqrt(mean(fivefold$SSEntitlement)/length(df$redistributionRate)) - sqrt(mean(subjectData$SSEntitlement)/length(df$redistributionRate)) #the change in root mean squared error, per trial in entilement condition
sqrt(mean(fivefold$SSCorruption)/length(df$redistributionRate)) - sqrt(mean(subjectData$SSCorruption)/length(df$redistributionRate)) #the change in root mean squared error, per trial in corruption condition
sqrt(mean(fivefold$SSLuck)/length(df$redistributionRate)) - sqrt(mean(subjectData$SSLuck)/length(df$redistributionRate)) #the change in root mean squared error, per trial in luck condition

fivefold$AICMerit = length(df$redistributionRate) * log(fivefold$SSMerit/length(df$redistributionRate)) + 2 * 2
fivefold$AICEntitlement = length(df$redistributionRate) * log(fivefold$SSEntitlement/length(df$redistributionRate)) + 2 * 2
fivefold$AICCorruption = length(df$redistributionRate) * log(fivefold$SSCorruption/length(df$redistributionRate)) + 2 * 2
fivefold$AICLuck = length(df$redistributionRate) * log(fivefold$SSLuck/length(df$redistributionRate)) + 2 * 2

t.test(fivefold$AICMerit, subjectData$AIC, paired = T) #test fivefold MFI against normal MFI for this model
t.test(fivefold$AICEntitlement, subjectData$AICEntitlement, paired = T) #test fivefold MFI against normal MFI for this model
t.test(fivefold$AICCorruption, subjectData$AICCorruption, paired = T) #test fivefold MFI against normal MFI for this model
t.test(fivefold$AICLuck, subjectData$AICLuck, paired = T) #test fivefold MFI against normal MFI for this model

library(lsa)
cosines = vector('numeric', length = length(conditions) * 2 * 5) #number of conditions times number of parameters times the number of folds
for (i in 1:5){
  for (j in 1:length(conditions)){
    k = ((i - 1) * 5) + j
    cosines[k] = cosine(subjectData$Parameter1, fivefold[, (k + 1)])
    cosines[(k+20)] = cosine(subjectData$Parameter1, fivefold[, (k + 21)])
  }
}

mean(cosines[1:5]) #cosine similarity of theta in merit condition
mean(cosines[6:10]) #cosine similarity of theta in entitlement condition
mean(cosines[11:15]) #cosine similarity of theta in corruption condition
mean(cosines[16:20]) #cosine similarity of theta in luck condition
mean(cosines[21:25]) #cosine similarity of phi in merit condition
mean(cosines[26:30]) #cosine similarity of phi in entitlement condition
mean(cosines[31:35]) #cosine similarity of phi in corruption condition
mean(cosines[36:40]) #cosine similarity of phi in luck condition


mean(cosines[1:20]) #cosine similarity of theta across all conditions
mean(cosines[21:40]) #cosine similarity of phi across all conditions

################ 3_1_0

alpha = 0.05/length(conditions)

# Merit Condition
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICPMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICLMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICTMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICPLMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICPTMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICMerit[-excludedM], altSubjectData$AICLTMerit[-excludedM], paired = T, conf.level = (1 - alpha), alternative = 'less')


# Corruption Condition
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICPEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICLEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICTEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICPLEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICPTEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICEntitlement[-excludedE], altSubjectData$AICLTEntitlement[-excludedE], paired = T, conf.level = (1 - alpha), alternative = 'less')

# Corruption Condition
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICPCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICLCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICTCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICPLCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICPTCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICCorruption[-excludedC], altSubjectData$AICLTCorruption[-excludedC], paired = T, conf.level = (1 - alpha), alternative = 'less')

# Luck Condition
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICPLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICLLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICTLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICPLLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICPTLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')
t.test(subjectData$AICLuck[-excludedL], altSubjectData$AICLTLuck[-excludedL], paired = T, conf.level = (1 - alpha), alternative = 'less')

################ 3_2_0

### Robust and Reliable Free Parameters

subjectData$phiMeritAdjusted = (subjectData$phiMerit - 0.5) #subtracting 0.5 puts 0 where phi weights on equality-seeking and equity-seeking in the same way (i.e. indifference point)
subjectData$phiMeritAdjusted = (subjectData$phiMeritAdjusted * (1 - subjectData$thetaMerit)) #the greater theta is, the closer to 0 this adjusted value is
subjectData$phiMeritAdjusted = subjectData$phiMeritAdjusted + 0.5 #returning to initial scale value (i.e. 0 to 1 with 0.5 beingthe indifference point)

subjectData$phiEntitlementAdjusted = ((subjectData$phiEntitlement - 0.5) * (1 - subjectData$thetaEntitlement)) + 0.5
subjectData$phiCorruptionAdjusted = ((subjectData$phiCorruption - 0.5) * (1 - subjectData$thetaCorruption)) + 0.5
subjectData$phiLuckAdjusted = ((subjectData$phiLuck - 0.5) * (1 - subjectData$thetaLuck)) + 0.5

### Meaningful Condition Differences

subjectData$SSAAO = vector('numeric')

for (i in 1:length(included_subjects)){
  datafile = paste(parentfolder, included_subjects[i], restoffilepath, sep = '') # produces a character vector 'parentfolder/included_subjects[i]**.filetype'
  df = read.csv2(datafile)
  df = df[, c(49, 40:48, 33)]
  df$redistributionRate = df$redistributionRate/100 #converting to a decimal from a percent
  result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   df = df)
  
  df$PredictionAAO = vector('numeric')
  df$ObservedOutcome = new_value(df$myself, df$redistributionRate)
  for (k in 1:length(df$redistributionRate)){
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(theta = result,
                           phi = phiPerCondition[j],
                           Equity = equity(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                           Equality = equality(new_value(df[k, 1:10], choices[n]), df[k, 1:10], choices[n]),
                           Payout = payout(new_value(df[k, 1], choices[n]), df[k, 1], choices[n]))
    }
    correct_choice = which(Utility == max(Utility))
    df$PredictionAAO[k] = new_value(df$myself[k], Choices[correct_choice[sample(length(correct_choice), 1)]])
  }
  subjectData$SSAAO[i] = sum((df$PredictionAAO - df$ObservedOutcome)**2)
}
N = length(df[, 1])
k = 2

subjectData$AICAAO = N * log(subjectData$SSAAO/N) + 2*2
subjectData$SSCond = subjectData$SSMerit + subjectData$SSEntitlement + subjectData$SSCorruption + subjectData$SSLuck
subjectData$AICCond = N * log(subjectData$SSCond/N) + 2*2*length(conditions)

excluded = which(subjectData$SSCond == 0)
t.test(subjectData$AICCond[-excluded], subjectData$AICAAO, paired = T, alternative = 'less')

### Testing a Modulatory Hypothesis

omnibusData = data.frame(Theta = c(subjectData$thetaMerit, subjectData$thetaEntitlement, subjectData$thetaCorruption, subjectData$thetaLuck),
                         Phi = c(subjectData$phiMeritAdjusted, subjectData$phiEntitlementAdjusted, subjectData$phiCorruptionAdjusted, subjectData$phiLuckAdjusted),
                         Condition = rep(conditions, each = length(subjectData$SubjectID)),
                         SubjectID = rep(subjectData$SubjectID, times = length(conditions)))

library(lme4)

omnibusModulatoryEffectTheta = lmer(data = omnibusData, Theta ~ Condition + (1 | SubjectID))
omnibusModulatoryEffectPhi = lmer(data = omnibusData, Phi ~ Condition + (1 | SubjectID))

### Using Categorical Clusters to Test Modulatory Hypotheses

subjectData$strategyMerit = as.factor(subjectData$strategyMerit)
subjectData$strategyEntitlement = as.factor(subjectData$strategyEntitlement)
subjectData$strategyCorruption = as.factor(subjectData$strategyCorruption)
subjectData$strategyLuck = as.factor(subjectData$strategyLuck)
strategies = levels(subjectData)

columns = 10:13 #columns of subject data correspond, in order, to the conditions in the conditions vector

for (i in 1:length(conditions)){
  for (j in 1:length(strategies)){
    group_by_condition[i, j] = sum(subjectData[, columns[i]] == strategies[j])
  }
}

colnames(group_by_condition) = strategies
rownames(group_by_condition) = conditions
chisq.test(group_by_condition)

################ 3_3_0

trialData$PredictedNID = vector('numeric', length(trialData$SubjectID))

for (k in 1:length(conditions)){
  
  z = which(trialData$condition == conditions[j])
  df = trialData[z, c(6, 12:20, 2)]
  
  result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   df = df)
  
  for (i in 1:length(df)){
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(theta = result$par[1],
                           phi = result$par[2],
                           Equity = equity(new_value(df[i, 1:10], choices[n]), df[i, 1:10], choices[n]),
                           Equality = equality(new_value(df[i, 1:10], choices[n]), df[i, 1:10], choices[n]),
                           Payout = payout(new_value(df[i, 1], choices[n]), df[i, 1], choices[n]))
    }
    trialData$PredictedNID[z[i]] = Choices[which(Utility == max(Utility))]
  }
}

subjectData$SSNIDMerit = vector('numeric', length(subjectDatasubjectData$SubjectID))
subjectData$SSNIDEntitlement = vector('numeric', length(subjectDatasubjectData$SubjectID))
subjectData$SSNIDCorruption = vector('numeric', length(subjectDatasubjectData$SubjectID))
subjectData$SSNIDLuck = vector('numeric', length(subjectDatasubjectData$SubjectID))

for (i in 1:length(subjectData$SubjectID)){
  trialsMerit = which(subjectData$SubjectID[i] == trialData$SubjectID && 'merit' == trialData$Condition)
  trialsEntitlement = which(subjectData$SubjectID[i] == trialData$SubjectID && 'entitlement' == trialData$Condition)
  trialsCorruption = which(subjectData$SubjectID[i] == trialData$SubjectID && 'corruption' == trialData$Condition)
  trialsLuck = which(subjectData$SubjectID[i] == trialData$SubjectID && 'luck' == trialData$Condition)
  
  subjectData$SSNIDMerit[i] = sum((trialData$observedTaxRate[trialsMerit] - trialData$PredictedNID[trialsMerit])**2)
  subjectData$SSNIDEntitlement[i] = sum((trialData$observedTaxRate[trialsEntitlement] - trialData$PredictedNID[trialsEntitlement])**2)
  subjectData$SSNIDCorruption[i] = sum((trialData$observedTaxRate[trialsCorruption] - trialData$PredictedNID[trialsCorruption])**2)
  subjectData$SSNIDLuck[i] = sum((trialData$observedTaxRate[trialsLuck] - trialData$PredictedNID[trialsLuck])**2)
}

subjectData$AICNIDMerit = length(trialsMerit) * log(subjectData$SSNIDMerit/length(trialsMerit)) + 2 * (2/length(subjectData$SubjectID))
subjectData$AICNIDEntitlement = length(trialsEntitlement) * log(subjectData$SSNIDEntitlement/length(trialsEntitlement)) + 2 * (2/length(subjectData$SubjectID))
subjectData$AICNIDCorruption = length(trialsCorruption) * log(subjectData$SSNIDCorruption/length(trialsCorruption)) + 2 * (2/length(subjectData$SubjectID))
subjectData$AICNIDLuck = length(trialsLuck) * log(subjectData$SSNIDLuck/length(trialsLuck)) + 2 * (2/length(subjectData$SubjectID))

t.test(subjectData$AICMerit, subjectData$AICNIDMerit, paired = T, alternative = 'less', conf.int = (1 - alpha))
t.test(subjectData$AICEntitlement, subjectData$AICNIDEntitlement, paired = T, alternative = 'less', conf.int = (1 - alpha))
t.test(subjectData$AICCorruption, subjectData$AICNIDCorruption, paired = T, alternative = 'less', conf.int = (1 - alpha))
t.test(subjectData$AICLuck, subjectData$AICNIDLuck, paired = T, alternative = 'less', conf.int = (1 - alpha))