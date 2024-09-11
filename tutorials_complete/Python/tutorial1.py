import numpy as np
import pandas as pd

def payout_maximization(investment, multiplier, returned):
    return ((investment * multiplier) - returned) / (investment * multiplier)

def inequity(investment, multiplier, returned, endowment):
    return((((investment * multiplier - returned)/(investment * multiplier + endowment - investment)) - 0.5)**2)

def guilt(investment, believed_multiplier, returned, multiplier):
    return ((((investment * believed_multiplier)/2 - returned) / (investment * multiplier))**2)

example = pd.DataFrame()

trialList = pd.DataFrame({
    'Investment': np.repeat(np.arange(1, 11), repeats=6),
    'Multiplier':  np.repeat([2, 4, 6], repeats=20),
    'Believed_Multiplier':  np.repeat(4, 60),
    'Endowment':  np.repeat(10, 60)
})

def utility(theta, phi, guilt, inequity, payout):
    return(theta*payout - (1-theta)*min(guilt + phi, inequity - phi))


freeParameters = pd.DataFrame({
    'theta': np.repeat(np.arange(0, 0.505, 0.005), repeats=101),
    'phi': np.tile(np.arange(-0.1, 0.102, 0.002), 101)
})

predictions = pd.DataFrame()

import random
non_specific = np.zeros(len(freeParameters))

for i in range(len(freeParameters)):
    [Theta, Phi] = freeParameters.iloc[i]

    for k in range(len(trialList)):
        [I, M, B, E] = trialList.iloc[k]
        Choices = list(range(0, I * M + 1, 1))

        Utility = [0] * len(Choices)
        for n in range(len(Choices)):
            Utility[n] = utility(theta=Theta,
                                phi=Phi,
                                guilt=guilt(I, B, Choices[n], M),
                                inequity=inequity(I, M, Choices[n], E),
                                payout=payout_maximization(I, M, Choices[n]))
        
        correct_choice = []
        for n in range(len(Choices)): 
            if Utility[n] == max(Utility):
                correct_choice.append(n)
        if len(correct_choice) > 1:
            correct_choice = random.sample(range(len(correct_choice)), 1)
            non_specific[i] += 1

        predictions[i, k] = Choices[correct_choice[0]]

from scipy.stats import norm

def obj_function(params, decisions, method="OLS"):
    [Theta, Phi] = params
    
    predicted_utility = np.zeros(len(trialList))
    observed_utility = np.zeros(len(trialList))

    for k in range(len(trialList)):
        [I, M, B, E] = trialList.iloc[k]
        R = int(decisions[k])
        Choices = np.arange(0, I * M + 1, 1)
        
        Utility = np.zeros(len(Choices))
        for n in range(len(Choices)):
            Utility[n] = utility(Theta, Phi, 
                                 guilt(I, B, Choices[n], M), 
                                 inequity(I, M, Choices[n], E), 
                                 payout_maximization(I, M, Choices[n]))
        
        predicted_utility[k] = np.max(Utility)
        observed_utility[k] = Utility[R]

    if method == "OLS":
        return np.sum((predicted_utility - observed_utility)**2)
    elif method == "MLE":
        sd = np.std(observed_utility - predicted_utility)
        return -1 * np.sum(norm.logpdf(observed_utility, loc=predicted_utility, scale=sd))
    
from scipy.optimize import minimize

initial_params = [0, 0]
lower_bounds = [0, -0.1]
upper_bounds = [0.5, 0.1]

theta_recovered = np.zeros(11**2)
phi_recovered = np.zeros(11**2)
theta_true = np.repeat(np.arange(0, 0.55, 0.05), 11)  
phi_true = np.tile(np.arange(-0.1, 0.12, 0.02), 11)

for i in range(len(theta_true)):
    this_idx = np.where((theta_true[i] == freeParameters['theta']) & (phi_true[i] == freeParameters['phi']))[0][0]
    decisions = predictions.loc[this_idx, :]
    decisionsList = []
    for j in range(len(trialList)):
        decisionsList.append(decisions.iloc[0, j])
    result = minimize(obj_function,
                      x0=initial_params,  
                      args=(decisionsList,),  
                      bounds=list(zip(lower_bounds, upper_bounds))) 

    theta_recovered[i] = result.x[0] 
    phi_recovered[i] = result.x[1]