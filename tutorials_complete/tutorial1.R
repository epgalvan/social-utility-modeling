################ 1_4_0

### Create Construct Value Functions

payout_maximization = function(investment, multiplier, returned){
  return(((investment * multiplier) - returned)/(investment * multiplier))
}
inequity = function(investment, multiplier, returned, endowment){
  return(((((investment * multiplier) - returned)/(endowment - investment + (investment * multiplier))) - 0.5)**2)
}
guilt = function(investment, believed_multiplier, returned, multiplier){
  return((((investment * believed_multiplier * 0.5) - returned)/(investment * multiplier))**2)
}

################ 1_5_0

### Preallocating and Defining
trialList = data.frame(Investment = rep(seq(1, 10, 1), times = 6), 
                       Multiplier = rep(c(2, 4, 6), each = 20), 
                       Believed_Multiplier = rep(4, 60), 
                       Endowment = rep(10, 60))

utility = function(theta, phi, guilt, inequity, payout){
  return(theta*payout - (1-theta)* (min(guilt + phi, inequity - phi)))
}

freeParameters = data.frame(theta = rep(seq(0, 0.5, 0.005), each = 101),
                            phi = rep(seq(-0.1, 0.1, 0.002), times = 101))

predictions = data.frame()

### Define All Loops

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
                           guilt = guilt(I, B, Choices[n], M),
                           inequity = inequity(I, M, Choices[n], E),
                           payout = payout_maximization(I, M, Choices[n]))
    }
    correct_choice = which(Utility == max(Utility))
    if (length(correct_choice) > 1){
      correct_choice = correct_choice[sample(seq(length(correct_choice)), 1)]
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
  observed_utility = vector('numeric', length(trialList[,1]))
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
  this_idx = which(round(theta_true[i] * 2, 2)/2 == round(freeParameters$theta * 2, 2)/2 & phi_true[i] == (round((freeParameters$phi + 0.1) * 5, 2)) -0.1)
  result = fmincon(obj_function,x0 = initial_params, lb = lower_bounds, ub = upper_bounds,
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

################ 1_7_0

### Cluster Your Data Using HAC

distance_mat = dist(predictions, method = 'euclidean')
set.seed(240)
hierarchical = hclust(distance_mat, method = 'average')
plot(hierarchical)
fit = cutree(hierarchical, k = 4)

### Identify Where Your Clusters Are

freeParameters$Strategy = as.character(fit)
freeParameters$Strategy[which(freeParameters$Strategy == freeParameters$Strategy[1])] = 'Guilt-Averse'
freeParameters$Strategy[which(freeParameters$Strategy == freeParameters$Strategy[10101])] = 'Greedy'
freeParameters$Strategy[which(freeParameters$Strategy == freeParameters$Strategy[100])] = 'Inequity-Averse'
freeParameters$Strategy[which(freeParameters$Strategy != 'Inequity-Averse' & freeParameters$Strategy != 'Greedy' & freeParameters$Strategy != 'Guilt-Averse')] = 'Moral Opportunists'; 
freeParameters$Strategy = as.factor(freeParameters$Strategy) #Strategy clusters
model_space = ggplot(data = freeParameters, aes(x = theta, y = phi, color = Strategy)) + 
  labs(x = 'Theta', y = 'Phi', color = 'Strategy') + geom_point(size = 2.5) + 
  scale_color_manual(values = c(rgb(50,50,200, maxColorValue = 255), rgb(230,157,54, maxColorValue = 255), rgb(57,193,59, maxColorValue = 255), 
                                rgb(200,50,50, maxColorValue = 255))); model_space

### Examine Model Predictions Efficiently

toPlot = data.frame()
for (i in 1:length(freeParameters[,1])){
  replacement = ((i - 1) * 60 + 1):(i * 60)
  toPlot[replacement, 1] = freeParameters$Strategy[i]
  toPlot[replacement, 2] = trialList$Investment
  toPlot[replacement, 3] = trialList$Multiplier
  toPlot[replacement, 4] = as.numeric(predictions[i,])
}
colnames(toPlot) = c('Strategy', 'Investment', 'Multiplier', 'Return')

ggplot(data = toPlot[which(toPlot$Multiplier==2),], aes(x = Investment, y = Return, group = Strategy, color = Strategy)) + geom_smooth(se = TRUE) + 
  scale_color_manual(values = c(rgb(50,50,200, maxColorValue = 255), rgb(230,157,54, maxColorValue = 255), rgb(57,193,59, maxColorValue = 255), 
                                rgb(200,50,50, maxColorValue = 255))) + lims(x = c(0, 10), y = c(0, 30))

ggplot(data = toPlot[which(toPlot$Multiplier==4),], aes(x = Investment, y = Return, group = Strategy, color = Strategy)) + geom_smooth(se = TRUE) + 
  scale_color_manual(values = c(rgb(50,50,200, maxColorValue = 255), rgb(230,157,54, maxColorValue = 255), rgb(57,193,59, maxColorValue = 255), 
                                rgb(200,50,50, maxColorValue = 255))) + lims(x = c(0, 10), y = c(0, 30))

ggplot(data = toPlot[which(toPlot$Multiplier==6),], aes(x = Investment, y = Return, group = Strategy, color = Strategy)) + geom_smooth(se = TRUE) + 
  scale_color_manual(values = c(rgb(50,50,200, maxColorValue = 255), rgb(230,157,54, maxColorValue = 255), rgb(57,193,59, maxColorValue = 255), 
                                rgb(200,50,50, maxColorValue = 255))) + lims(x = c(0, 10), y = c(0, 30))

################ 2_1_0

### Preallocating and Defining Functions

included_subjects = c(t(read.csv2('C:/Users/DELL/Downloads/tutorial1_Data/subjectsIncluded_batch1.csv', sep=',', header = F)), 
                      t(read.csv2('C:/Users/DELL/Downloads/tutorial1_Data/subjectsIncluded_batch2.csv', sep=',', header = F)))

trialData = read.csv2("C:/Users/DELL/Downloads/tutorial1_Data/allDataLong.csv", sep =',')
trialData$Prediction = vector('numeric', length(trialData$Subject))
trialData$Strategy = vector('numeric', length(trialData$Subject))
trialData = trialData[-which(trialData$Investment == 0), ]
trialData = trialData[-which(is.na(as.numeric(trialData$Returned))),]


subjectData = data.frame()
trialList = read.csv2("C:/Users/DELL/Downloads/tutorial1_Data/trialSet.csv", sep =',', )
trialList = trialList[,2:3]
trialList$Believed_Multiplier = 4
trialList$Endowment = 10

obj_function = function(params, df, method = "OLS") {
  Theta = params[1]
  Phi = params[2]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  chosen = as.numeric(df[,4]) + 1
  for (k in 1:length(df[,1])){
    I = df[k, 2]
    M = df[k, 3]
    B = 4
    E = 10
    if (I > 10) {Choices = seq(0, (I*M), round((I*M)/10))} else {Choices = seq(0, (I * M), 1)}
    
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(theta = Theta, 
                           phi = Phi, 
                           guilt = guilt(I, B, Choices[n], M), 
                           inequity = inequity(I, M, Choices[n], E), 
                           payout = payout_maximization(I, M, Choices[n]))
    }
    predicted_utility[k] = max(Utility)
    observed_utility[k] = Utility[chosen[k]]
  }
  if (method == "OLS"){
    return(sum((predicted_utility - observed_utility)**2))
  } else if (method == "MLE"){
    return(-1 * sum(dnorm(observed_utility, mean = predicted_utility, sd = sd, log = TRUE)))
  }
} #must redefine since we are using a diferent trialList

### Recover Free Parameters and Determine Predicted Decisions for these Free Parameters

for (i in 1:length(included_subjects)){
  df = trialData[which(included_subjects[i] == trialData$Subject), ]
  result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                   lb = lower_bounds, ub = upper_bounds,
                   df = df)
  
  closestPoint = which(as.numeric(freeParameters[,1], 3) == (round(result$par[1]/2, 2))*2 & 
                         ((round((as.numeric(freeParameters[,2]) + 0.1)*5, 2))/5) - 0.1 == ((round((result$par[2] + 0.1)*5, 2))/5) - 0.1)
  for (k in 1:length(df$Returned)){
    I = df$Investment[k]
    M = df$Multiplier[k]
    R = df$Returned[k]
    B = 4
    E = 10
    if (I > 10) {Choices = seq(0, (I*M), round((I*M)/10))} else {Choices = seq(0, (I * M), 1)}
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility(theta = result$par[1], 
                           phi = result$par[2], 
                           guilt = guilt(I, B, Choices[n], M), 
                           inequity = inequity(I, M, Choices[n], E), 
                           payout = payout_maximization(I, M, Choices[n]))
    }
    correct_choice = which(Utility == max(Utility))
    if (length(correct_choice) > 1){
      correct_choice = correct_choice[sample(1:length(correct_choice), 1)]
    }
    df$Prediction[k] = Choices[correct_choice]
  }
  
  model_NLL = -2 * log(sum(dnorm(as.numeric(df$Returned), mean = df$Prediction)))
  model_SS = sum((as.numeric(df$Returned) - df$Prediction)**2)
  
  subjectData[i, 1:6] = c(included_subjects[i], result$par[1], result$par[2], freeParameters$Strategy[closestPoint], model_NLL, model_SS) 

  trialData$Prediction[which(included_subjects[i] == trialData$Subject)] = df$Prediction
  trialData$Strategy[which(included_subjects[i] == trialData$Subject)] = freeParameters$Strategy[closestPoint]
}
colnames(subjectData) = c('SubjectID', 'Theta', 'Phi', 'Strategy', 'modelNLL', 'modelSS')

################ 2_2_0

### Compute Model Fit Index for Each Subject

N = length(trialList[, 1])
k = 2
subjectData$modelAIC = N * log(subjectData$modelSS/N) + 2*k

################ 2_3_0

### Preallocating and Defining Functions
utility_greed = function(greed){
  return(greed)
}
utility_guilt = function(theta, greed, guilt){
  return(theta * greed - (1 - theta) * guilt)
}
utility_inequity = function(theta, greed, inequity){
  return(theta * greed - (1 - theta) * inequity)
}

obj_function_guilt = function(params, df, method = "OLS") {
  Theta = params[1]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  chosen = as.numeric(df[,4]) + 1
  for (k in 1:length(df[,1])){
    I = df[k, 2]
    M = df[k, 3]
    B = 4
    E = 10
    if (I > 10) {Choices = seq(0, (I*M), round((I*M)/10))} else {Choices = seq(0, (I * M), 1)}
    
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility_guilt(theta = Theta, 
                                 greed = payout_maximization(I, M, Choices[n]), 
                                 guilt = guilt(I, B, Choices[n], M))
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
obj_function_inequity = function(params, df, method = "OLS") {
  Theta = params[1]
  
  predicted_utility = vector('numeric', length(df[,1]))
  observed_utility = vector('numeric', length(df[,1]))
  chosen = as.numeric(df[,4]) + 1
  for (k in 1:length(df[,1])){
    I = df[k, 2]
    M = df[k, 3]
    B = 4
    E = 10
    if (I > 10) {Choices = seq(0, (I*M), round((I*M)/10))} else {Choices = seq(0, (I * M), 1)}
    
    Utility = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      Utility[n] = utility_inequity(theta = Theta, 
                                    greed = payout_maximization(I, M, Choices[n]), 
                                    inequity = inequity(I, M, Choices[n], E))
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

altSubjectData = data.frame()

### Recover Free Parameters and Determine Predicted Decisions Per Subject

for (i in 1:length(included_subjects)){
  df = trialData[which(included_subjects[i] == trialData$Subject), ]
  
  result_guilt = optim(fn = obj_function_guilt, par = 0.5, lower = 0, upper = 1, df = df, method = 'L-BFGS-B')
  result_inequity = optim(fn = obj_function_inequity, par = 0.5, lower = 0, upper = 1, df = df, method = 'L-BFGS-B')
  
  df$PredictionGreed = vector('numeric', length(df$Subject))
  df$PredictionGuilt = df$PredictionAlt1
  df$PredictionInequity = df$PredictionAlt1
  for (k in 1:length(df$Returned)){
    I = df$Investment[k]
    M = df$Multiplier[k]
    R = df$Returned[k]
    B = 4
    E = 10
    if (I > 10) {Choices = seq(0, (I*M), round((I*M)/10))} else {Choices = seq(0, (I * M), 1)}
    UtilityGreed = vector('numeric', length(Choices))
    UtilityGuilt = vector('numeric', length(Choices))
    UtilityInequity = vector('numeric', length(Choices))
    for (n in 1:length(Choices)){
      UtilityGreed[n] = utility_greed(greed = payout_maximization(I, M, Choices[n]))
      UtilityGuilt[n] = utility_guilt(theta = result_guilt$par[1],
                                     greed = payout_maximization(I, M, Choices[n]),
                                     guilt = guilt(I, B, Choices[n], M))
      UtilityInequity[n] = utility_inequity(theta = result_guilt$par[1],
                                        greed = payout_maximization(I, M, Choices[n]),
                                        inequity = inequity(I, M, Choices[n], E))
    }
    correct_choice_greed = which(UtilityGreed == max(UtilityGreed))
    correct_choice_guilt = which(UtilityGuilt == max(UtilityGuilt))
    correct_choice_inequity = which(UtilityInequity == max(UtilityInequity))
    if (length(correct_choice_greed) > 1){
      correct_choice_greed = correct_choice_greed[sample(correct_choice_greed, 1)]
    }
    if (length(correct_choice_guilt) > 1){
      correct_choice_guilt = correct_choice_guilt[sample(correct_choice_guilt, 1)]
    }
    if (length(correct_choice_inequity) > 1){
      correct_choice_inequity = correct_choice_inequity[sample(correct_choice_inequity, 1)]
    }
    df$PredictionGreed[k] = Choices[correct_choice_greed]
    df$PredictionGuilt[k] = Choices[correct_choice_guilt]
    df$PredictionInequity[k] = Choices[correct_choice_inequity]
  }
  
  model_NLL_Greed = -2 * log(sum(dnorm(as.numeric(df$Returned), mean = df$PredictionGreed)))
  model_SS_Greed = sum((as.numeric(df$Returned) - df$PredictionGreed)**2)
  model_NLL_Guilt = -2 * log(sum(dnorm(as.numeric(df$Returned), mean = df$PredictionGuilt)))
  model_SS_Guilt = sum((as.numeric(df$Returned) - df$PredictionGuilt)**2)
  model_NLL_Inequity = -2 * log(sum(dnorm(as.numeric(df$Returned), mean = df$PredictionInequity)))
  model_SS_Inequity = sum((as.numeric(df$Returned) - df$PredictionInequity)**2)
  
  altSubjectData[i, 1:9] = c(included_subjects[i], result_guilt$par[1], result_inequity$par[1], model_NLL_Greed, model_SS_Greed, model_NLL_Guilt, model_SS_Guilt, model_NLL_Inequity, model_SS_Inequity)
}
colnames(altSubjectData) = c('SubjectID', 'guilt_theta', 'inequity_theta', 'greed_modelNLL', 'greed_modelSS', 'guilt_modelNLL', 'guilt_modelSS', 'inequity_ModelNLL', 'inequity_ModelSS')

### Compute Model Fit Indices

altSubjectData$modelAICGreed = N * log(altSubjectData$greed_modelSS/N) + 2*0
altSubjectData$modelAICGuilt = N * log(altSubjectData$guilt_modelSS/N) + 2*1
altSubjectData$modelAICInequity = N * log(altSubjectData$inequity_ModelSS/N) + 2*1

### Compare Model Performance

#subjects excluded from analyses not included in dataset
averageAIC = c(mean(subjectData$modelAIC), mean(altSubjectData$modelAICGreed), mean(altSubjectData$modelAICGuilt), mean(altSubjectData$modelAICInequity))
fullAIC = length(trialData$Subject) * log(sum(subjectData$modelSS)/length(trialData$Subject)) + (2 * k * length(subjectData$Subject))
fullAICGreed = length(trialData$Subject) * log(sum(altSubjectData$greed_modelSS)/length(trialData$Subject)) + (2 * 0 * length(subjectData$Subject))
fullAICGuilt = length(trialData$Subject) * log(sum(altSubjectData$guilt_modelSS)/length(trialData$Subject)) + (2 * 1 * length(subjectData$Subject))
fullAICInequity = length(trialData$Subject) * log(sum(altSubjectData$inequity_ModelSS)/length(trialData$Subject)) + (2 * 1 * length(subjectData$Subject))

bestModel = c("Moral Strategies Model", "Greed Model", "Guilt Model", "Inequity Model")[which(averageAIC == min(averageAIC))] #best model based on average performance per subject
bestModelFullDataset = c("Moral Strategies Model", "Greed Model", "Guilt Model", "Inequity Model")[which(c(fullAIC, fullAICGreed, fullAICGuilt, fullAICInequity) == min(c(fullAIC, fullAICGreed, fullAICGuilt, fullAICInequity)))] #best model based on all data observations

################ 2_4_0

### Model Performance

modelPredictions = lm(data = trialData, Returned ~ Prediction)
summary(modelPredictions) # R-squared

### Visually Checking Assumptions

qplot(x = trialData$Prediction, y = as.numeric(trialData$Returned), geom = 'smooth') + 
  geom_abline(slope = 1, intercept = 0) + labs(x = 'Prediction', y = 'Observed') + lims(x = c(0, 30), y = c(0, 30)) # linearity

normvals = rnorm(1000, mean = 0, sd = sd(trialData$Prediction - as.numeric(trialData$Returned)))
qplot(x = trialData$Prediction - as.numeric(trialData$Returned), geom = 'density', bw = 1, color = 'Actual') + 
  geom_density(aes(x = normvals, color = 'Predicted'), bw = 1) + labs(x = 'Prediction - Observed', y = 'Density') #normality

qplot(x = trialData$Investment, y = (trialData$Prediction-as.numeric(trialData$Returned)), group = as.factor(trialData$Multiplier), color = as.factor(trialData$Multiplier), geom = 'smooth')  + 
  labs(x = 'Investment', y = 'Prediction - Observed', color = 'Multiplier') #independence of error

qplot(x = trialData$Investment, y = (trialData$Prediction-as.numeric(trialData$Returned)), geom = 'smooth') + 
  labs(x = 'Investment', y = 'Prediction - Observed')#homoscedasticity

## Fivefold Validation

fivefold = data.frame() #preallocate for parameters and errors from the fivefold validation to go into

for (i in 1:length(included_subjects)){
  df = trialData[which(included_subjects[i] == trialData$Subject), ]
  order = sample(length(df$Returned))
  Theta_ff = vector('numeric', length = 5)
  Phi_ff = vector('numeric', length = 5)
  SS_ff = 0
  Prediction_ff = vector('numeric', length(df$Returned))
  for (z in 1:5){
    j = round((z - 1) * (length(df$Returned)/5) + 1)
    n = round(z * (length(df$Returned)/5))
    withheld = order[j:n]
    m = ((i - 1) * 5) + z
    
    result_ff = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                        lb = lower_bounds, ub = upper_bounds,
                        df = df[-withheld,])
    
    Theta_ff[z] = result_ff$par[1]
    Phi_ff[z] = result_ff$par[2]
    for (n in 1:length(withheld)){
      if (df$Investment[withheld[n]] > 10) {
        Choices = seq(0, (df$Investment[withheld[n]] * df$Multiplier[withheld[n]]), round((df$Investment[withheld[n]]*df$Multiplier[withheld[n]])/10))
      } else {
        Choices = seq(0, (df$Investment[withheld[n]] * df$Multiplier[withheld[n]]), 1)
      }
      utility_choices = vector('numeric', length(Choices))
      for (q in 1:length(Choices)){
        utility_choices[q] = utility(theta = result_ff$par[1], 
                             phi = result_ff$par[2],
                             guilt = guilt(df$Investment[withheld[n]], believed_multiplier = 4, Choices[q], df$Multiplier[withheld[n]]),
                             payout = payout_maximization(df$Investment[withheld[n]], df$Multiplier[withheld[n]], Choices[q]),
                             inequity = inequity(df$Investment[withheld[n]], df$Multiplier[withheld[n]], Choices[q], endowment = 10))
      }
      Prediction_ff[withheld[n]] = Choices[which(utility_choices == max(utility_choices))[1]]
    }
  }
  SS_ff = sum((as.numeric(df$Returned) - Prediction_ff)**2)
  fivefold[i, 1:12] = c(SS_ff, Theta_ff, Phi_ff, included_subjects[i])
}
colnames(fivefold) = c('SS', 'Par1_fold1', 'Par1_fold2', 'Par1_fold3', 'Par1_fold4', 'Par1_fold5', 
                       'Par2_fold1', 'Par2_fold2', 'Par2_fold3', 'Par2_fold4', 'Par2_fold5', 'SubjectID')

sqrt(mean(fivefold$SS)/length(df$Investment)) - sqrt(mean(subjectData$modelSS)/length(df$Investment)) #the change in root mean squared error, per trial
fivefold$AIC = length(df$Investment) * log(fivefold$SS/length(df$Investment)) + 2 * 2
t.test(fivefold$AIC, subjectData$modelAIC, paired = T) #test fivefold MFI against normal MFI for this model

library(lsa)
cosines = vector('numeric', length = 10)
for (i in 1:5){
  cosines[i] = cosine(subjectData$Theta, fivefold[, (i + 1)]) #to get the correct columns in the fivefold dataframe (2-6)
  cosines[(i+5)] = cosine(subjectData$Phi, fivefold[, (i + 6)]) #to get the correct columns in the fivefold dataframe (7-11)
}

mean(cosines[1:5]) #cosine similarity of parameter 1
mean(cosines[6:10]) #cosine similarity of parameter 2

################ 3_1_0

t.test(subjectData$modelAIC, altSubjectData$modelAICGreed, paired = T) #negative mean difference/t value means that first term is lower than second term and lower MFI is better
t.test(subjectData$modelAIC, altSubjectData$modelAICGuilt, paired = T)
t.test(subjectData$modelAIC, altSubjectData$modelAICInequity, paired = T)

library(ggsignif)
aic = c(mean(subjectData$modelAIC), mean(altSubjectData$modelAICGreed), mean(altSubjectData$modelAICGuilt), mean(altSubjectData$modelAICInequity))
qplot(y = aic,
      x = as.factor(c('Moral Strategies Model', 'Greed Model', 'Guilt Model', 'Inequity Model')), 
      fill = as.factor(c('Moral Strategies Model', 'Greed Model', 'Guilt Model', 'Inequity Model')), 
      color = '',
      geom = 'col') + 
  labs(x = 'Model', y = 'AIC', fill = NULL, color = NULL) + 
  theme_minimal() + scale_color_manual(values = c(rgb(0, 0, 0, maxColorValue = 255))) + 
  scale_fill_manual(values = c(rgb(132.5, 132.5, 132.5, maxColorValue = 255), 
                               rgb(132.5, 132.5, 132.5, maxColorValue = 255), 
                               rgb(132.5, 132.5, 132.5, maxColorValue = 255), 
                               rgb(218, 165, 32, maxColorValue = 255))) + 
  geom_signif(comparisons = list(c('Moral Strategies Model', 'Greed Model')), y = 341, textsize = 0)+
  geom_signif(comparisons = list(c('Moral Strategies Model', 'Guilt Model')), y = 328, textsize = 0)+
  geom_signif(comparisons = list(c('Moral Strategies Model', 'Inequity Model')), y = 315, textsize = 0) +
  annotate("text", x = 2.5, label = "***", y = 358, size = 4)+ # 3 for p < 0.001, 2 for p < 0.01, 1 for p < 0.05 usually
  annotate("text", x = 3, label = "***", y = 345, size = 4)+
  annotate("text", x = 3.5, label = "***", y = 332, size = 4)

################ 3_2_0

#Not Applicable

################ 3_3_0

result = fmincon(obj_function,x0 = initial_params, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
                 lb = lower_bounds, ub = upper_bounds,
                 df = trialData)

trialData$PredictedNID = vector('numeric', length(trialData$Subject))
for (i in 1:length(trialData$Subject)){
  if (trialData$Investment[i] > 10) {
    Choices = seq(0, (trialData$Investment[i] * trialData$Multiplier[i]), round((trialData$Investment[i]*trialData$Multiplier[i])/10))
  } else {
    Choices = seq(0, (trialData$Investment[i] * trialData$Multiplier[i]), 1)
  }
  Utility = vector('numeric', length(Choices))
  for (n in 1:length(Choices)){
    Utility[n] = utility(theta = result$par[1],
                         phi = result$par[2],
                         payout = payout_maximization(trialData$Investment[i], trialData$Multiplier[i], Choices[n]),
                         guilt = guilt(trialData$Investment[i], believed_multiplier = 4, Choices[n], trialData$Multiplier[i]),
                         inequity = inequity(trialData$Investment[i], trialData$Multiplier[i], Choices[n], endowment = 10))
  }
  trialData$PredictedNID[i] = Choices[which(Utility == max(Utility))]
}

subjectData$SS_NID = vector('numeric', length(subjectData$SubjectID))
for (i in 1:length(subjectData$SubjectID)){
  trials = which(subjectData$SubjectID[i] == trialData$Subject)
  subjectData$SS_NID[i] = sum((as.numeric(trialData$Returned[trials]) - trialData$PredictedNID[trials])**2)
}

# number of parameter divided by the number of people (i.e. number of parameters for each person)
subjectData$AIC_NID = length(df$Investment) * log(subjectData$SS_NID/length(df$Investment)) + 2 * (2/length(subjectData$SubjectID))

t.test(subjectData$modelAIC, subjectData$AIC_NID, paired = T)

summary(lm(data = trialData, Returned ~ PredictedNID))

aic_id = c(mean(subjectData$modelAIC), mean(subjectData$AIC_NID))

qplot(y = aic_id,
      x = as.factor(c('Individual Differences', 'No Individual Differences')), 
      fill = as.factor(c('Individual Differences', 'No Individual Differences')),
      color = '',
      geom = 'col') + 
  labs(x = 'Model', y = 'AIC', fill = NULL) + lims(y = c(-5, 155)) +
  theme_minimal() + 
  geom_signif(comparisons = list(c('Individual Differences', 'No Individual Differences')), y_position = 135, textsize = 0, tip_length = 0.15)+
  annotate("text", x = 1.5, y = 135, 
           label = "p = 0.18", vjust = -1, size = 4) + 
  scale_fill_manual(values = c(rgb(218, 165, 32, maxColorValue = 255), rgb(132.5, 132.5, 132.5, maxColorValue = 255))) + 
  scale_color_manual(values = c(rgb(0, 0, 0, maxColorValue = 255)))

################ 3_4_0