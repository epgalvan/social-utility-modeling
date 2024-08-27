################ 1_4_0

### Create Construct Value Functions

inequality = function(a1, b1, a2, b2){
  d_Inequality = abs(a2 - b2) - abs(a1 - b1)
  return(d_Inequality)
}
harm = function(a0, b0, a1, b1, a2, b2){
  initial = c(a0, b0)
  choice1 = c(a1, b1)
  choice2 = c(a2, b2)
  if (a0 == b0){advantaged = 1} else {advantaged = which(initial == max(initial))}
  d_lossAdvantaged = (initial[advantaged] - choice1[advantaged]) - (initial[advantaged] - choice2[advantaged])
  return(d_lossAdvantaged)
}
rankReverse = function(a0, b0, a1, b1, a2, b2){
  if (a0 == b0){return(0)}
  d_initial = a0 - b0
  d_choice1 = a1 - b1
  d_choice2 = a2 - b2
  if (d_initial > 0){
    if (d_choice1 < 0){choice1Reversed = 1} else {choice1Reversed = 0}
    if (d_choice2 < 0){choice2Reversed = 1} else {choice2Reversed = 0}
  } else if (d_initial < 0){
    if (d_choice1 > 0){choice1Reversed = 1} else {choice1Reversed = 0}
    if (d_choice2 > 0){choice2Reversed = 1} else {choice2Reversed = 0}
  }
  return(choice1Reversed - choice2Reversed)
}

################ 1_5_0

### Preallocating and Defining Functions, TrialList, and Parameters

trialList = data.frame(a0 = sample(10:20, 100, replace = T))
trialList$b0 = 20 - trialList$a0
for (i in 1:nrow(trialList)){
  trialList$a1 = sample(5:trialList$a0[i], 100, replace = T)
}
trialList$b1 = 20 - trialList$a1

for (i in 1:nrow(trialList)){
  trialList$a2 = sample(5:trialList$a0[i], 100, replace = T); trialList$a2[which(trialList$a2 == trialList$a1)] = 10
}
trialList$b2 = 20 - trialList$a2
trialList = rbind(trialList, trialList[,c(2,1,4,3,6,5)])

#trialList = read.csv2("C:/Users/DELL/Downloads/Data/Data/HPP_fMRI_beh_data_for_lmm.csv", sep = ',')
#trialList = trialList[which(trialList$subject == 101 & trialList$trail_type == 3), c(6,7,13,14,21,22)]
#colnames(trialList) = c('a0', 'b0', 'a1', 'b1', 'a2', 'b2')

utility = function(pars, IVs){
  IVS = as.numeric(IVs)
  a0 = IVs[1]; b0 = IVs[2]; a1 = IVs[3]; b1 = IVs[4]; a2 = IVs[5]; b2 = IVs[6]
  alpha = pars[1]; delta = pars[2]; rho = pars[3]
  inequality = inequality(a1, b1, a2, b2)
  harm = harm(a0, b0, a1, b1, a2, b2)
  rank = rankReverse(a0, b0, a1, b1, a2, b2)
  return((alpha * inequality) - (delta * harm) - (rho * rank))
}

probability = function(pars, utilitydiff){
  beta = pars[length(pars)-2]; epsilon = pars[length(pars)-1]; gamma = pars[length(pars)]
  prob = 1 / (1 + exp(-(beta * utilitydiff)))
  prob = prob * (1 - 2 * epsilon) + epsilon + gamma * (2 * epsilon)
  return(prob)
}

freeParameters = data.frame(alpha = rep(seq(0, 1.9, 0.1), times = 60) + sample(seq(0, 0.1, 0.001), 20*6*10, replace = T), #ranging from 0 to 2
                            delta = rep(seq(0, 1.9, 0.1), times = 60) + sample(seq(0, 0.1, 0.001), 20*6*10, replace = T), #ranging from 0 to 2
                            rho = rep(seq(0, 1.9, 0.1), times = 60) + sample(seq(0, 0.1, 0.001), 20*6*10, replace = T), #ranging from 0 to 2
                            beta = sample(seq(0,10), 20*6*10, replace = T), # stochasticity
                            epsilon = rep(rep(seq(0, 0.5, 0.1), times = 10), each = 20), # noise
                            gamma = rep(seq(-0.5, 0.4, 0.1), each = 20*6) + sample(seq(0, 0.1, 0.001), 20*6*10, replace = T)) # bias
predictions = data.frame()

generatePredictions = function(parameters, df){
  pred = vector('numeric', nrow(df))
  for (i in 1:nrow(df)){
    thisTrialIVs = as.numeric(df[i, ])
    utilityDiff = utility(parameters, thisTrialIVs)
    pred[i] = max(min(probability(parameters, utilityDiff), 0.9999999999), 0.00000000001)
  }
  return(pred)
}

### Determine Predictions

for (i in 1:nrow(freeParameters)){
  pars = freeParameters[i, ]
  predictions[i, 1:nrow(trialList)] = as.numeric(round(generatePredictions(pars, trialList) * 1000)> sample(1000,1))
}

################ 1_6_0

### Objective Functions

obj_function = function(params, df, optimMethod = "MLE") { 
  Prob1 = generatePredictions(params, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}

### Optimizers

library(pracma)

initial_params = c(1,1,1, 4, 0.25, 0)
lower_bounds = c(0, 0, 0, 0, 0, -0.5)
upper_bounds = c(2,2,2, 10, 0.5, 0.5)
freeParameters$alphaRecovered = 0
freeParameters$deltaRecovered = 0
freeParameters$rhoRecovered = 0
freeParameters$betaRecovered = 0
freeParameters$epsilonRecovered = 0
freeParameters$gammaRecovered = 0

optimize = function(obj, initial_params, lower_bounds, upper_bounds, df){
  tryCatch({
    fmincon(obj, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "MLE", tol = 1e-08)
  }, error = function(e){
    fmincon(obj, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = df, optimMethod = "OLS", tol = 1e-08)
  })
}

for (i in 1:nrow(freeParameters)) {
  trialList$Predictions = as.numeric(predictions[i,])
  result = optimize(obj_function, initial_params, lower_bounds, upper_bounds, trialList)
  freeParameters$alphaRecovered[i] = result$par[1]
  freeParameters$deltaRecovered[i] = result$par[2]
  freeParameters$rhoRecovered[i] = result$par[3]
  freeParameters$betaRecovered[i] = result$par[4]
  freeParameters$epsilonRecovered[i] = result$par[5]
  freeParameters$gammaRecovered[i] = result$par[6]
  if (mod(i, nrow(freeParameters)/100) == 0){message(round(100* (i/(nrow(freeParameters)))), '% there', sep = '')}
}

### Verifying Free Parameter Recovery Process

library(ggplot2)

freeParameters$Epsilon = factor(freeParameters$epsilon)

qplot(data = freeParameters, x = alpha, y = alphaRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline()
qplot(data = freeParameters, x = delta, y = deltaRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline()
qplot(data = freeParameters, x = rho, y = rhoRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline()
qplot(data = freeParameters, x = beta, y = betaRecovered, color = Epsilon) + geom_smooth(se = F)+ geom_abline()
qplot(data = freeParameters, x = epsilon, y = epsilonRecovered) + geom_smooth(method = 'lm') + geom_abline()
qplot(data = freeParameters, x = gamma, y = gammaRecovered, color = Epsilon) + geom_smooth(se = F) + geom_abline() 

################ 2_1_0

### Preallocating and Defining Functions

trialData = read.csv2("C:/Users/DELL/Downloads/Data/Data/HPP_fMRI_beh_data_for_lmm.csv", sep = ',')
trialData = trialData[which(trialData$trail_type==3), c(1,6,7,13,14,21,22,36)]
colnames(trialData) = c('SubjectID', 'a0', 'b0', 'a1', 'b1', 'a2', 'b2', 'Chose1')
trialData$Chose1 = trialData$Chose1 - 1
trialData$Prob1 = 0
included_subjects = levels(factor(trialData$SubjectID))

subjectData = data.frame()

grab_data = function(subject){
  return(trialData[which(trialData$SubjectID == subject), -1])
}

addPredictions = function(trialData, subject, predictions){
  trialData$Prob1[which(trialData$SubjectID == subject)] = predictions
  return(trialData)
}

### Recover Free Parameters and Define Predicted Decisions for these Free Parameters

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  result = optimize(obj_function, initial_params, lower_bounds, upper_bounds, df)
  
  # Just Added
  df$Prob1 = generatePredictions(result$par, df)
  
  model_SS = sum((df$Chose1 - df$Prob1)**2)
  model_NLL = -2*sum(df$Chose1 * log(df$Prob1) + (1 - df$Chose1) * log(1 - df$Prob1))
  
  subjectData[i, 1:9] = c(included_subjects[i], result$par[1], result$par[2], result$par[3], result$par[4], result$par[5], result$par[6],
                     model_SS, model_NLL)
  
  trialData = addPredictions(trialData, included_subjects[i], df$Prob1)
}

colnames(subjectData) = c("subjectID", "Alpha", "Delta", "Rho", "Beta", "Epsilon", "Gamma", "SS", "Deviance")

################ 2_2_0

subjectData$BIC = as.numeric(subjectData$Deviance) + log(65) * 6

################ 2_3_0

### Preallocating and Defining Functions

#adjusting

of_alphaOnly = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(1, 4:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_deltaOnly = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(2, 4:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_rhoOnly = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(3:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_ad = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(1:2, 4:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_ar = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(1, 3:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_dr = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(2:6)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_noEpsilon = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(1:4)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_noGamma = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[c(1:5)] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}
of_GammaOnly = function(params, df, optimMethod = "MLE"){
  params_new = rep(0, 6); params_new[5] = 0.5; params_new[6] = params
  Prob1 = generatePredictions(params_new, df); Chose1 = df[,7]
  if (optimMethod == "OLS"){return(sum((Chose1 - Prob1)**2)) 
  } else if (optimMethod == "MLE"){return(-sum(Chose1 * log(Prob1) + (1 - Chose1) * log(1 - Prob1)))}
}

ofs = list(of_alphaOnly, of_deltaOnly, of_rhoOnly, of_ad, of_ar, of_dr, of_noEpsilon, of_noGamma, of_GammaOnly)
idxs = list(c(1,4:6), c(2,4:6), c(3:6), c(1:2, 4:6), c(1, 3:6), c(2:6), c(1:4), c(1:5), c(6))

altSubjectData = data.frame()
altTrialData = trialData[,-9]
altTrialData$alphaOnly_Prob1 = 0
altTrialData$deltaOnly_Prob1 = 0
altTrialData$rhoOnly_Prob1 = 0
altTrialData$ad_Prob1 = 0
altTrialData$ar_Prob1 = 0
altTrialData$dr_Prob1 = 0
altTrialData$noEpsilon_Prob1 = 0
altTrialData$noGamma_Prob1 = 0
altTrialData$gammaOnly_Prob1 = 0

### Recover Free Parameters and Determine Predicted Decisions 

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  outputs = vector('numeric', length = sum(length(c(c(1,4:6), c(2,4:6), c(3:6), c(1:2, 4:6), c(1, 3:6), c(2:6), c(1:4), c(1:5), c(6)))) + 2*9)
  j = 0
  for (k in 1:length(idxs)){
    idx = idxs[k][[1]]; initials = initial_params[idx]; uppers = upper_bounds[idx]; lowers = lower_bounds[idx]; of = ofs[k][[1]]
    if (length(idx) > 1){
      result = optimize(of, initials, lowers, uppers, df)
      } else {
      result = optim(par = initials, fn = of, method = "L-BFGS-B", lower = lowers, upper = uppers, df = df)
    }
    pars = rep(0, times = 6)
    pars[idx] = result$par
    df$Prob1 = 0; df$Prob1 = generatePredictions(pars, df)
    
    model_SS = sum((df$Chose1 - df$Prob1)**2)
    model_NLL = -2*sum(df$Chose1 * log(df$Prob1) + (1 - df$Chose1) * log(1 - df$Prob1))
    outputs[(j+1):(j+2+length(result$par))] = c(result$par, model_SS, model_NLL)
    j = j+2+length(result$par)
    altTrialData[which(altTrialData$SubjectID == included_subjects[i]), (8+k)] = df$Prob1
  }
  
  altSubjectData[i, 1:56] = c(included_subjects[i], outputs)
}

colnames(altSubjectData) = c('SubjectID', 
                             'Alpha_M1', 'Beta_M1', 'Epsilon_M1', 'Gamma_M1', 'SS_M1', 'Deviance_M1',
                             'Delta_M2', 'Beta_M2', 'Epsilon_M2', 'Gamma_M2', 'SS_M2', 'Deviance_M2', 
                             'Rho_M3', 'Beta_M3', 'Epsilon_M3', 'Gamma_M3', 'SS_M3', 'Deviance_M3', 
                             'Alpha_M4', 'Delta_M4', 'Beta_M4', 'Epsilon_M4', 'Gamma_M4', 'SS_M4', 'Deviance_M4',
                             'Alpha_M5', 'Rho_M5', 'Beta_M5', 'Epsilon_M5', 'Gamma_M5', 'SS_M5', 'Deviance_M5',
                             'Delta_M6', 'Rho_M6', 'Beta_M6', 'Epsilon_M6', 'Gamma_M6', 'SS_M6', 'Deviance_M6',
                             'Alpha_M7', 'Delta_M7', 'Rho_M7', 'Beta_M7', 'SS_M7', 'Deviance_M7',
                             'Alpha_M8', 'Delta_M8', 'Rho_M8', 'Beta_M8', 'Epsilon_M8', 'SS_M8', 'Deviance_M8',
                             'Gamma_M9', 'SS_M9', 'Deviance_M9')

for (i in 2:ncol(altSubjectData)){altSubjectData[,i] = as.numeric(altSubjectData[,i])}

### Now Compute BIC

altSubjectData$BIC_M1 = (altSubjectData$Deviance_M1) + log(65) * 4
altSubjectData$BIC_M2 = (altSubjectData$Deviance_M2) + log(65) * 4
altSubjectData$BIC_M3 = (altSubjectData$Deviance_M3) + log(65) * 4
altSubjectData$BIC_M4 = (altSubjectData$Deviance_M4) + log(65) * 5
altSubjectData$BIC_M5 = (altSubjectData$Deviance_M5) + log(65) * 5
altSubjectData$BIC_M6 = (altSubjectData$Deviance_M6) + log(65) * 5
altSubjectData$BIC_M7 = (altSubjectData$Deviance_M7) + log(65) * 4
altSubjectData$BIC_M8 = (altSubjectData$Deviance_M8) + log(65) * 5
altSubjectData$BIC_M9 = (altSubjectData$Deviance_M9) + log(65) * 1

### Now Compare BIC

modelBIC = c(sum(subjectData$BIC), sum(altSubjectData$BIC_M1), sum(altSubjectData$BIC_M2), sum(altSubjectData$BIC_M3), sum(altSubjectData$BIC_M4), sum(altSubjectData$BIC_M5), sum(altSubjectData$BIC_M6), sum(altSubjectData$BIC_M7), sum(altSubjectData$BIC_M8), sum(altSubjectData$BIC_M9))
which(modelBIC == min(modelBIC))

################ 2_4_0

### Assessing Model Performance

sum(trialData$Chose1 == round(trialData$Prob1))/nrow(trialData)

qplot(subjectData$BIC, geom = 'density')
worstExplained = which(subjectData$BIC > as.numeric(summary(subjectData$BIC)[4]))

qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = a0-a1, y = Chose1 - Prob1, group = SubjectID) + geom_smooth(se=F)
qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = a0-a2, y = Chose1 - Prob1, group = SubjectID) + geom_smooth(se=F)
qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = b0-b1, y = Chose1 - Prob1, group = SubjectID) + geom_smooth(se=F)
qplot(data = trialData[which(trialData$SubjectID %in% altSubjectData$SubjectID[worstExplained]),], x = b0-b2, y = Chose1 - Prob1, group = SubjectID) + geom_smooth(se=F)


### Checking Assumptions

ggplot(data = trialData) + geom_smooth(aes(x = Prob1, y = Chose1, group = a0 < b0, color = factor(a0 < b0))) 
ggplot(data = trialData) + geom_density(aes(x = Chose1, group = a0 < b0, fill = factor(a0 < b0), alpha = 0.5))
normvals = rnorm(1000, mean = 0, sd = sd(trialData$Prob1 - trialData$Chose1))
qplot(x = trialData$Prob1 - trialData$Chose1, geom = 'density', bw = sd(trialData$Prob1 - trialData$Chose1), color = 'Actual') +
  geom_density(aes(x = normvals, color = 'Predicted'), bw = sd(trialData$Prob1 - trialData$Chose1))
qplot(x = abs(trialData$a1 - trialData$a2), y = trialData$Prob1 - trialData$Chose1, geom = 'jitter') + geom_smooth(method = 'loess')
qplot(x = abs(trialData$a1 - trialData$a2), y = trialData$Prob1 - trialData$Chose1, color = factor(trialData$a0 > trialData$b0), geom = 'jitter') + geom_smooth(method = 'loess')

### Assessing Independence

library(lme4)
library(MuMIn)

ris_model = glmer(data = trialData, Chose1 ~ Prob1 + (1 + Prob1 | SubjectID), family = "binomial")
r.squaredGLMM(ris_model)

ri_model = glmer(data = trialData, Chose1 ~ Prob1 + (1 | SubjectID), family = "binomial")
r.squaredGLMM(ri_model)

ric_model = glmer(data = trialData, Chose1 ~ Prob1 + factor(a0<b0) + (1 | SubjectID), family = "binomial")
r.squaredGLMM(ric_model)

### Fivefold Validation

fivefold = data.frame() #preallocate for parameters and errors from the fivefold validation to go into
trialData$Prob1_ff = 0

for (i in 1:length(included_subjects)){
  df = grab_data(included_subjects[i])
  df$Prob1 = 0
  
  order = sample(nrow(df))
  A_ff = vector('numeric', length = 5)
  D_ff = vector('numeric', length = 5)
  R_ff = vector('numeric', length = 5)
  B_ff = vector('numeric', length = 5)
  E_ff = vector('numeric', length = 5)
  G_ff = vector('numeric', length = 5)
  Deviance_ff = 0
  df$Pred = 0
  for (z in 1:5){
    j = round((z - 1) * (nrow(df)/5) + 1)
    n = round(z * (nrow(df)/5))
    withheld = order[j:n]
    
    result_ff = optimize(obj_function, initial_params, lower_bounds, upper_bounds, df[-withheld,])
    
    A_ff[z] = result_ff$par[1]
    D_ff[z] = result_ff$par[2]
    R_ff[z] = result_ff$par[3]
    B_ff[z] = result_ff$par[4]
    E_ff[z] = result_ff$par[5]
    G_ff[z] = result_ff$par[6]
    df$Prob1[withheld] = generatePredictions(result_ff$par, df[withheld, ])
  }
  Deviance_ff = -2*sum(df$Chose1 * log(df$Prob1) + (1 - df$Chose1) * log(1 - df$Prob1))
  fivefold[i, 1:32] = c(included_subjects[i], Deviance_ff, A_ff, D_ff, R_ff, B_ff, E_ff, G_ff)
  trialData$Prob1_ff[which(trialData$SubjectID == included_subjects[i])] = df$Prob1
}
colnames(fivefold) = c('SubjectID', 'Deviance', 'A_F1', 'A_F2', 'A_F3', 'A_F4', 'A_F5',
                       'D_F1','D_F2','D_F3','D_F4','D_F5', 'R_F1', 'R_F2', 'R_F3', 'R_F4', 'R_F5',
                       'B_F1','B_F2','B_F3','B_F4','B_F5', 'E_F1', 'E_F2', 'E_F3', 'E_F4', 'E_F5', 
                       'G_F1', 'G_F2', 'G_F3', 'G_F4', 'G_F5')

sum(round(trialData$Prob1_ff) == trialData$Chose1)/nrow(trialData)
fivefold$BIC = as.numeric(fivefold$Deviance) + log(65) * 6
t.test(fivefold$BIC, subjectData$BIC, paired = T) #test fivefold MFI against normal MFI for this model

library(lsa)
cosines = vector('numeric', length = 60)
for (i in 1:5){
  cosines[i] = cosine(as.numeric(subjectData$Alpha), as.numeric(fivefold[, (i + 2)]))
  cosines[(i+5)] = cosine(as.numeric(subjectData$Delta), as.numeric(fivefold[, (i + 7)]))
  cosines[(i+10)] = cosine(as.numeric(subjectData$Rho), as.numeric(fivefold[, (i + 12)]))
  cosines[(i+15)] = cosine(as.numeric(subjectData$Beta), as.numeric(fivefold[, (i + 17)]))
  cosines[(i+20)] = cosine(as.numeric(subjectData$Epsilon), as.numeric(fivefold[, (i + 22)]))
  cosines[(i+25)] = cosine(as.numeric(subjectData$Gamma), as.numeric(fivefold[, (i + 27)]))
}

mean(cosines[1:5]) 
mean(cosines[6:10])
mean(cosines[11:15]) 
mean(cosines[16:20]) 
mean(cosines[21:25]) 
mean(cosines[26:30])  

################ 3_1_0

t.test(subjectData$BIC, altSubjectData$BIC_M1, paired = T)
t.test(subjectData$BIC, altSubjectData$BIC_M2, paired = T)
t.test(subjectData$BIC, altSubjectData$BIC_M3, paired = T)
t.test(subjectData$BIC, altSubjectData$BIC_M4, paired = T)
t.test(subjectData$BIC, altSubjectData$BIC_M5, paired = T)
t.test(subjectData$BIC, altSubjectData$BIC_M6, paired = T)

################ 3_3_0

resultNID = fmincon(obj_function, x0 = initial_params, lb = lower_bounds, ub = upper_bounds, df = trialData[2:8], optimMethod = "MLE")
trialData$Prob1_NID = generatePredictions(resultNID$par, trialData)
sum(trialData$Chose1 == round(trialData$Prob1_NID))/nrow(trialData)

altSubjectData$Deviance_NID = 0

for (i in 1:length(included_subjects)){ 
  trials = which(included_subjects[i] == trialData$SubjectID)
  df = trialData[trials, ]
  altSubjectData$Deviance_NID[i] = -2*sum(df$Chose1 * log(df$Prob1_NID) + (1 - df$Chose1) * log(1 - df$Prob1_NID))
}

altSubjectData$BIC_NID = as.numeric(altSubjectData$Deviance_NID) + log(65) * 6/(length(included_subjects))

t.test(subjectData$BIC, altSubjectData$BIC_NID, paired = T)
which(modelBIC > sum(altSubjectData$BIC_NID))
