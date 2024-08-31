################ 1_4_0

### Create Construct Value Functions

import numpy as np
import pandas as pd

def inequality(a1, b1, a2, b2):
    d_Inequality = abs(a2 - b2) - abs(a1 - b1)
    return d_Inequality

def harm(a0, b0, a1, b1, a2, b2):
    initial = np.array([a0, b0])
    choice1 = np.array([a1, b1])
    choice2 = np.array([a2, b2])
    advantaged = 1 if a0 == b0 else np.argmax(initial)
    d_lossAdvantaged = (initial[advantaged] - choice1[advantaged]) - (initial[advantaged] - choice2[advantaged])
    return d_lossAdvantaged

def rankReverse(a0, b0, a1, b1, a2, b2):
    if a0 == b0:
        return 0
    d_initial = a0 - b0
    d_choice1 = a1 - b1
    d_choice2 = a2 - b2
    if d_initial > 0:
        choice1Reversed = 1 if d_choice1 < 0 else 0
        choice2Reversed = 1 if d_choice2 < 0 else 0
    else:
        choice1Reversed = 1 if d_choice1 > 0 else 0
        choice2Reversed = 1 if d_choice2 > 0 else 0
    return choice1Reversed - choice2Reversed

################ 1_5_0

### Preallocating and Defining Functions, TrialList, and Parameters

trialList = pd.DataFrame({
    'a0': np.random.randint(10, 21, 100)
})
trialList['b0'] = 20 - trialList['a0']
trialList['a1'] = [np.random.randint(5, a0 + 1) for a0 in trialList['a0']]
trialList['b1'] = 20 - trialList['a1']
trialList['a2'] = [np.random.randint(5, a0 + 1) for a0 in trialList['a0']]
trialList.loc[trialList['a2'] == trialList['a1'], ('a2')] = 10
trialList['b2'] = 20 - trialList['a2']
trialList = pd.concat([trialList, trialList[['b0', 'a0', 'b1', 'a1', 'b2', 'a2']]], ignore_index=True)

def utility(pars, IVs):
    IVS = np.array(IVs)
    a0, b0, a1, b1, a2, b2 = IVs[:6]
    alpha, delta, rho = pars[:3]
    ineq = inequality(a1, b1, a2, b2)
    hrm = harm(a0, b0, a1, b1, a2, b2)
    rank = rankReverse(a0, b0, a1, b1, a2, b2)
    return (alpha * ineq) - (delta * hrm) - (rho * rank)

def probability(pars, utilitydiff):
    beta, epsilon, gamma = pars[-3:]
    prob = 1 / (1 + np.exp(-(beta * utilitydiff)))
    prob = prob * (1 - 2 * epsilon) + epsilon + gamma * (2 * epsilon)
    return prob

freeParameters = pd.DataFrame({
    'alpha': np.tile(np.arange(0, 2, 0.1), 60) + np.random.choice(np.arange(0, 0.1, 0.001), 20 * 6 * 10),
    'delta': np.tile(np.arange(0, 2, 0.1), 60) + np.random.choice(np.arange(0, 0.1, 0.001), 20 * 6 * 10),
    'rho': np.tile(np.arange(0, 2, 0.1), 60) + np.random.choice(np.arange(0, 0.1, 0.001), 20 * 6 * 10),
    'beta': np.random.choice(np.arange(0, 11), 20 * 6 * 10),
    'epsilon': np.repeat(np.repeat(np.arange(0, 0.6, 0.1), 10), 20),
    'gamma': np.repeat(np.arange(-0.5, 0.5, 0.1), 20 * 6) + np.random.choice(np.arange(0, 0.1, 0.001), 20 * 6 * 10)
})
predictions = pd.DataFrame()

def generatePredictions(parameters, df):
    pred = np.zeros(len(df))
    for i in range(len(df)):
        thisTrialIVs = df.iloc[i].to_numpy()
        utilityDiff = utility(parameters, thisTrialIVs)
        pred[i] = max(min(probability(parameters, utilityDiff), 0.9999999999), 0.00000000001)
    return pred

### Determine Predictions

for i in range(len(freeParameters)):
    pars = freeParameters.iloc[i]
    predictions.loc[i, range(0, len(trialList))] = generatePredictions(pars, trialList)

################ 1_6_0

### Objective Functions

def obj_function(params, df, optimMethod="MLE"):
    Prob1 = generatePredictions(params, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1) ** 2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

from scipy.optimize import minimize

def optimize(obj, initial_params, lower_bounds, upper_bounds, df):
    try:
        result = minimize(obj, initial_params, args=(df,), bounds=list(zip(lower_bounds, upper_bounds)), tol=1e-08)
    except:
        result = minimize(obj, initial_params, args=(df,), bounds=list(zip(lower_bounds, upper_bounds)), tol=1e-08, method="L-BFGS-B")
    return result

freeParameters[['alphaRecovered', 'deltaRecovered', 'rhoRecovered', 'betaRecovered', 'epsilonRecovered', 'gammaRecovered']] = 0.0

initial_params = [1, 1, 1, 4, 0.25, 0]
lower_bounds = [0, 0, 0, 0, 0, -0.5]
upper_bounds = [2, 2, 2, 10, 0.5, 0.5]

for i in range(len(freeParameters)):
    trialList['Predictions'] = predictions.loc[i].to_numpy()
    result = optimize(obj_function, initial_params, lower_bounds, upper_bounds, trialList)
    freeParameters.loc[i, ['alphaRecovered', 'deltaRecovered', 'rhoRecovered', 'betaRecovered', 'epsilonRecovered', 'gammaRecovered']] = result.x
    if i % (len(freeParameters) // 100) == 0:
        print(f"{round(100 * (i / len(freeParameters)))}% there")

### Verifying Free Parameter Recovery Process

import seaborn as sns
import matplotlib.pyplot as plt

freeParameters['Epsilon'] = freeParameters['epsilon'].astype('category')

sns.lmplot(data=freeParameters, x='alpha', y='alphaRecovered', hue='Epsilon', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['alpha'].min(), freeParameters['alpha'].max()],
         [freeParameters['alphaRecovered'].min(), freeParameters['alphaRecovered'].max()], 'k--')
plt.show()

sns.lmplot(data=freeParameters, x='delta', y='deltaRecovered', hue='Epsilon', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['delta'].min(), freeParameters['delta'].max()],
         [freeParameters['deltaRecovered'].min(), freeParameters['deltaRecovered'].max()], 'k--')
plt.show()

sns.lmplot(data=freeParameters, x='rho', y='rhoRecovered', hue='Epsilon', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['rho'].min(), freeParameters['rho'].max()],
         [freeParameters['rhoRecovered'].min(), freeParameters['rhoRecovered'].max()], 'k--')
plt.show()

sns.lmplot(data=freeParameters, x='beta', y='betaRecovered', hue='Epsilon', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['beta'].min(), freeParameters['beta'].max()],
         [freeParameters['betaRecovered'].min(), freeParameters['betaRecovered'].max()], 'k--')
plt.show()

sns.lmplot(data=freeParameters, x='epsilon', y='epsilonRecovered', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['epsilon'].min(), freeParameters['epsilon'].max()],
         [freeParameters['epsilonRecovered'].min(), freeParameters['epsilonRecovered'].max()], 'k--')
plt.show()

sns.lmplot(data=freeParameters, x='gamma', y='gammaRecovered', hue='Epsilon', lowess=True, ci=None, scatter_kws={'s': 50})
plt.plot([freeParameters['gamma'].min(), freeParameters['gamma'].max()],
         [freeParameters['gammaRecovered'].min(), freeParameters['gammaRecovered'].max()], 'k--')
plt.show()

################ 2_1_0

### Preallocating and Defining Functions

trialData = pd.read_csv("C:/Users/DELL/Downloads/Data/Data/HPP_fMRI_beh_data_for_lmm.csv", sep=',')
trialData = trialData.iloc[np.where(trialData['trail_type'] == 3)[0].tolist(), :]
trialData = trialData.iloc[:, [0,5,6,12,13,20,21,35]]
trialData.columns = ['SubjectID', 'a0', 'b0', 'a1', 'b1', 'a2', 'b2', 'Chose1']
trialData['Chose1'] -= 1
trialData['Prob1'] = 0.0
included_subjects = trialData['SubjectID'].unique()
subjectData = pd.DataFrame()

def grab_data(subject):
    return trialData[trialData['SubjectID'] == subject].drop(columns=['SubjectID'])

def addPredictions(trialData, subject, predictions):
    trialData.loc[trialData['SubjectID'] == subject, 'Prob1'] = predictions
    return trialData

### Recover Free Parameters and Define Predicted Decisions for these Free Parameters

for i in range(0, len(included_subjects)):
    df = grab_data(included_subjects[i])
    result = optimize(obj_function, initial_params, lower_bounds, upper_bounds, df)
    df['Prob1'] = generatePredictions(result.x, df)
    model_SS = np.sum((df['Chose1'] - df['Prob1']) ** 2)
    model_NLL = -2 * np.sum(df['Chose1'] * np.log(df['Prob1']) + (1 - df['Chose1']) * np.log(1 - df['Prob1']))
    subjectData.loc[i, range(0, 9)] = [included_subjects[i]] + result.x.tolist() + [model_SS, model_NLL]
    trialData = addPredictions(trialData, included_subjects[i], df['Prob1'])

subjectData.columns = ["subjectID", "Alpha", "Delta", "Rho", "Beta", "Epsilon", "Gamma", "SS", "Deviance"]

################ 2_2_0

subjectData['BIC'] = subjectData['Deviance'] + np.log(65) * 6

################ 2_3_0

### Preallocating and Defining Functions

#adjusting

def of_alphaOnly(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[0, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_deltaOnly(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[1, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_rhoOnly(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[2, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_ad(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[0, 1, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_ar(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[0, 2, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_dr(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[1, 2, 3, 4, 5]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_noEpsilon(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[0, 1, 2, 3]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_noGamma(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[[0, 1, 2, 3, 4]] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

def of_GammaOnly(params, df, optimMethod="MLE"):
    params_new = np.zeros(6)
    params_new[4] = 0.5
    params_new[5] = params
    Prob1 = generatePredictions(params_new, df)
    Chose1 = df.iloc[:, 6]
    if optimMethod == "OLS":
        return np.sum((Chose1 - Prob1)**2)
    elif optimMethod == "MLE":
        return -np.sum(Chose1 * np.log(Prob1) + (1 - Chose1) * np.log(1 - Prob1))

ofs = [of_alphaOnly, of_deltaOnly, of_rhoOnly, of_ad, of_ar, of_dr, of_noEpsilon, of_noGamma, of_GammaOnly]
idxs = [[0, 3, 4, 5], [1, 3, 4, 5], [2, 3, 4, 5], [0, 1, 3, 4, 5], [0, 2, 3, 4, 5], [1, 2, 3, 4, 5], [0, 1, 2, 3], [0, 1, 2, 3, 4], [5]]

altSubjectData = pd.DataFrame()
altTrialData = trialData.drop(columns=[trialData.columns[8]])

### Recover Free Parameters and Determine Predicted Decisions 

for i in range(0, len(included_subjects)):
    df = grab_data(included_subjects[i])
    outputs = []
    j = 0
    
    for k in range(0, len(idxs)):
        idx = idxs[k]
        initials = [initial_params[i] for i in idx]
        uppers = [upper_bounds[i] for i in idx]
        lowers = [lower_bounds[i] for i in idx]
        of = ofs[k]
        
        result = optimize(of, initials, lowers, uppers, df)

        pars = np.zeros(6)
        pars[idx] = result.x
        if len(idx) == 1:
            pars[4] = 0.5
        
        df['Prob1'] = generatePredictions(pars, df)
        
        model_SS = np.sum((df['Chose1'] - df['Prob1'])**2)
        model_NLL = -2 * np.sum(df['Chose1'] * np.log(df['Prob1']) + (1 - df['Chose1']) * np.log(1 - df['Prob1']))
        outputs = outputs + result.x.tolist() + [model_SS, model_NLL]
        j += 2 + len(result.x)
        
        altTrialData.loc[altTrialData['SubjectID'] == included_subjects[i], 8 + k] = df['Prob1']
    
    altSubjectData.loc[i, range(0, 56)] = [included_subjects[i]] + outputs

altSubjectData.columns = [
    'SubjectID', 
    'Alpha_M1', 'Beta_M1', 'Epsilon_M1', 'Gamma_M1', 'SS_M1', 'Deviance_M1',
    'Delta_M2', 'Beta_M2', 'Epsilon_M2', 'Gamma_M2', 'SS_M2', 'Deviance_M2', 
    'Rho_M3', 'Beta_M3', 'Epsilon_M3', 'Gamma_M3', 'SS_M3', 'Deviance_M3', 
    'Alpha_M4', 'Delta_M4', 'Beta_M4', 'Epsilon_M4', 'Gamma_M4', 'SS_M4', 'Deviance_M4',
    'Alpha_M5', 'Rho_M5', 'Beta_M5', 'Epsilon_M5', 'Gamma_M5', 'SS_M5', 'Deviance_M5',
    'Delta_M6', 'Rho_M6', 'Beta_M6', 'Epsilon_M6', 'Gamma_M6', 'SS_M6', 'Deviance_M6',
    'Alpha_M7', 'Delta_M7', 'Rho_M7', 'Beta_M7', 'SS_M7', 'Deviance_M7',
    'Alpha_M8', 'Delta_M8', 'Rho_M8', 'Beta_M8', 'Epsilon_M8', 'SS_M8', 'Deviance_M8',
    'Gamma_M9', 'SS_M9', 'Deviance_M9'
]

altTrialData.columns = ['SubjectID', 'a0', 'b0', 'a1', 'b1', 'a2', 'b2', 'Chose1','alphaOnly_Prob1', 'deltaOnly_Prob1', 'rhoOnly_Prob1', 'ad_Prob1', 'ar_Prob1', 'dr_Prob1', 'noEpsilon_Prob1', 'noGamma_Prob1', 'gammaOnly_Prob1']

for col in altSubjectData.columns[1:]:
    altSubjectData[col] = pd.to_numeric(altSubjectData[col])

### Now Compute BIC

altSubjectData['BIC_M1'] = altSubjectData['Deviance_M1'] + np.log(65) * 4
altSubjectData['BIC_M2'] = altSubjectData['Deviance_M2'] + np.log(65) * 4
altSubjectData['BIC_M3'] = altSubjectData['Deviance_M3'] + np.log(65) * 4
altSubjectData['BIC_M4'] = altSubjectData['Deviance_M4'] + np.log(65) * 5
altSubjectData['BIC_M5'] = altSubjectData['Deviance_M5'] + np.log(65) * 5
altSubjectData['BIC_M6'] = altSubjectData['Deviance_M6'] + np.log(65) * 5
altSubjectData['BIC_M7'] = altSubjectData['Deviance_M7'] + np.log(65) * 4
altSubjectData['BIC_M8'] = altSubjectData['Deviance_M8'] + np.log(65) * 5
altSubjectData['BIC_M9'] = altSubjectData['Deviance_M9'] + np.log(65) * 1

### Now Compare BIC

modelBIC = [subjectData['BIC'].sum(), altSubjectData['BIC_M1'].sum(), altSubjectData['BIC_M2'].sum(), altSubjectData['BIC_M3'].sum(), altSubjectData['BIC_M4'].sum(), altSubjectData['BIC_M5'].sum(), altSubjectData['BIC_M6'].sum(), altSubjectData['BIC_M7'].sum(), altSubjectData['BIC_M8'].sum(), altSubjectData['BIC_M9'].sum()]
print(np.argmin(modelBIC)) 

################ 2_4_0

### Assessing Model Performance

print(np.sum(altTrialData['Chose1'] == np.round(altTrialData['ad_Prob1'])) / len(altTrialData))
altTrialData['a0-a1'] = altTrialData['a0']-altTrialData['a1']
altTrialData['a0-a2'] = altTrialData['a0']-altTrialData['a2']
altTrialData['b0-b1'] = altTrialData['b0']-altTrialData['b1']
altTrialData['b0-b2'] = altTrialData['b0']-altTrialData['b2']
altTrialData['a1-a2'] = altTrialData['a1']-altTrialData['a2']
altTrialData['b1-b2'] = altTrialData['b1']-altTrialData['b2']
altTrialData['a0_less_than_b0'] = trialData['a0'] < trialData['b0']
altTrialData['Chose1 - Prob1'] = altTrialData['Chose1']-altTrialData['ad_Prob1']

plt.figure(figsize=(8, 6))
sns.kdeplot(altSubjectData['BIC_M4'])
plt.title('Density Plot of BIC Values')
plt.xlabel('BIC')
plt.ylabel('Density')
plt.show()

bic_summary = altSubjectData['BIC_M4'].describe()
worst_explained = altSubjectData.index[altSubjectData['BIC_M4'] > bic_summary['75%']].tolist()
worst_subject_ids = altSubjectData['SubjectID'].iloc[worst_explained]
filtered_trialData = altTrialData[altTrialData['SubjectID'].isin(worst_subject_ids)]

sns.lmplot(data=filtered_trialData, x='a0-a1', y='Chose1 - Prob1', hue='SubjectID', legend=False, lowess=True)
plt.title('a0 - a1 vs Chose1 - Prob1')
plt.xlabel('a0 - a1')
plt.ylabel('Chose1 - Prob1')
plt.show()

sns.lmplot(data=filtered_trialData, x='a0-a2', y='Chose1 - Prob1', hue='SubjectID', legend=False, lowess=True)
plt.title('a0 - a2 vs Chose1 - Prob1')
plt.xlabel('a0 - a2')
plt.ylabel('Chose1 - Prob1')
plt.show()

sns.lmplot(data=filtered_trialData, x='b0-b1', y='Chose1 - Prob1', hue='SubjectID', legend=False, lowess=True)
plt.title('b0 - b1 vs Chose1 - Prob1')
plt.xlabel('b0 - b1')
plt.ylabel('Chose1 - Prob1')
plt.show()

sns.lmplot(data=filtered_trialData, x='b0-b2', y='Chose1 - Prob1', hue='SubjectID', legend=False, lowess=True)
plt.title('b0 - b2 vs Chose1 - Prob1')
plt.xlabel('b0 - b2')
plt.ylabel('Chose1 - Prob1')
plt.show()

sns.lmplot(data=altTrialData, x='ad_Prob1', y='Chose1', hue='a0_less_than_b0', lowess = True)
plt.title('Smooth Plot of Prob1 vs Chose1')
plt.xlabel('Prob1')
plt.ylabel('Chose1')
plt.show()

plt.figure(figsize=(10, 6))
sns.kdeplot(data=altTrialData, x='Chose1', hue='a0_less_than_b0', fill=True, alpha=0.5, palette='coolwarm')
plt.title('Density Plot of Chose1 by a0 < b0')
plt.xlabel('Chose1')
plt.ylabel('Density')
plt.legend(title='a0 < b0', labels=['False', 'True'])
plt.show()

residuals = altTrialData['ad_Prob1'] - trialData['Chose1']
altTrialData['residuals'] = residuals
normvals = np.random.normal(loc=0, scale=np.std(residuals), size=1000)

sns.kdeplot(residuals, bw_adjust=1, label='Actual', color='blue')
sns.kdeplot(normvals, bw_adjust=1, label='Predicted', color='red')
plt.title('Density Plot of Residuals vs Normal Distribution')
plt.xlabel('Residuals')
plt.ylabel('Density')
plt.legend()
plt.show()

sns.lmplot(data=altTrialData, x='a1-a2', y='residuals', ci=None, lowess=True)
plt.title('Abs(a1 - a2) vs Residuals with Loess Smoothing')
plt.xlabel('Abs(a1 - a2)')
plt.ylabel('Residuals')
plt.show()

sns.lmplot(data=altTrialData, x='a1-a2', y='residuals', hue='a0_less_than_b0', palette='coolwarm', lowess=True)
plt.title('Abs(a1 - a2) vs Residuals Colored by a0 > b0 with Loess Smoothing')
plt.xlabel('Abs(a1 - a2)')
plt.ylabel('Residuals')
plt.legend(title='a0 > b0', labels=['False', 'True'])
plt.show()

### Assessing Independence

# can do directly in R

### Fivefold Validation

fivefold = pd.DataFrame()
trialData['Prob1_ff'] = 0.0
initial_params = [1, 1, 4, 0.25, 0]
lower_bounds = [0, 0, 0, 0, -0.5]
upper_bounds = [2, 2, 10, 0.5, 0.5]

for i in range(0, len(included_subjects)):
    df = grab_data(included_subjects[i])
    df = df.reset_index()
    df['Prob1'] = 0.0
    
    order = np.random.permutation(len(df))
    
    A_ff, D_ff, B_ff, E_ff, G_ff = np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5)
    Deviance_ff = 0.0
    
    for z in range(5):
        j = int((z) * (len(df) / 5))
        n = int((z+1) * (len(df) / 5))
        withheld = order[j:n]

        input = df.iloc[~df.index.isin(withheld)]
        input = input.loc[:, ['a0', 'b0', 'a1', 'b1', 'a2', 'b2', 'Chose1', 'Prob1']]
        
        result_ff = optimize(of_ad, initial_params, lower_bounds, upper_bounds, input)
        
        A_ff[z], D_ff[z], B_ff[z], E_ff[z], G_ff[z] = result_ff.x
        input = df.iloc[df.index.isin(withheld)]
        input = input.loc[:, ['a0', 'b0', 'a1', 'b1', 'a2', 'b2', 'Chose1', 'Prob1']]
        pars = result_ff.x[0:2].tolist() + [0] + result_ff.x[2:].tolist()
        df.loc[df.index.isin(withheld),'Prob1'] = generatePredictions(pars, input)
    
    Deviance_ff = -2 * np.sum(df['Chose1'] * np.log(df['Prob1']) + (1 - df['Chose1']) * np.log(1 - df['Prob1']))
    fivefold.loc[i, range(0, 27)] = [included_subjects[i], Deviance_ff] + list(A_ff) + list(D_ff) + list(B_ff) + list(E_ff) + list(G_ff)
    trialData.loc[trialData['SubjectID'] == included_subjects[i], 'Prob1_ff'] = df['Prob1']

fivefold.columns = ['SubjectID', 'Deviance', 'A_F1', 'A_F2', 'A_F3', 'A_F4', 'A_F5',
                    'D_F1', 'D_F2', 'D_F3', 'D_F4', 'D_F5', 'B_F1', 'B_F2', 'B_F3', 'B_F4', 'B_F5', 
                    'E_F1', 'E_F2', 'E_F3', 'E_F4', 'E_F5', 'G_F1', 'G_F2', 'G_F3', 'G_F4', 'G_F5']

fivefold['BIC'] = fivefold['Deviance'] + np.log(65) * 6
print(sum(round(trialData['Prob1_ff']) == trialData['Chose1'])/len(trialData))

from scipy.stats import ttest_rel
print(ttest_rel(fivefold['BIC'], altSubjectData['BIC_M4']))

from sklearn.metrics.pairwise import cosine_similarity

cosines = []
for i in range(5):
    cosines.append(cosine_similarity([altSubjectData['Alpha_M4'].values], [fivefold.iloc[:, i + 2].values])[0][0])
    cosines.append(cosine_similarity([altSubjectData['Delta_M4'].values], [fivefold.iloc[:, i + 7].values])[0][0])
    cosines.append(cosine_similarity([altSubjectData['Beta_M4'].values], [fivefold.iloc[:, i + 10].values])[0][0])
    cosines.append(cosine_similarity([altSubjectData['Epsilon_M4'].values], [fivefold.iloc[:, i + 17].values])[0][0])
    cosines.append(cosine_similarity([altSubjectData['Gamma_M4'].values], [fivefold.iloc[:, i + 22].values])[0][0])

cosines = np.array(cosines)
print(np.mean(cosines[:5]))
print(np.mean(cosines[5:10]))
print(np.mean(cosines[10:15]))
print(np.mean(cosines[15:20]))
print(np.mean(cosines[20:25]))

################ 3_1_0

print(ttest_rel(altSubjectData['BIC_M4'], subjectData['BIC']))
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_M1']))
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_M2']))
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_M3']))
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_M5']))
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_M6']))

################ 3_3_0

resultNID = optimize(of_ad, initial_params, lower_bounds, upper_bounds, trialData.iloc[:, 1:8])
pars = resultNID.x[0:2].tolist() + [0] + resultNID.x[2:].tolist()
trialData['Prob1_NID'] = generatePredictions(pars, trialData)
print(np.sum(trialData['Chose1'] == np.round(trialData['Prob1_NID'])) / len(trialData))
altSubjectData['Deviance_NID'] = 0.0
for i, subject in enumerate(included_subjects):
    trials = trialData['SubjectID'] == subject
    df = trialData[trials]
    altSubjectData.loc[i, 'Deviance_NID'] = -2 * np.sum(df['Chose1'] * np.log(df['Prob1_NID']) + (1 - df['Chose1']) * np.log(1 - df['Prob1_NID']))

altSubjectData['BIC_NID'] = altSubjectData['Deviance_NID'] + np.log(65) * 6 / len(included_subjects)
print(ttest_rel(altSubjectData['BIC_M4'], altSubjectData['BIC_NID']))
print(np.where(modelBIC > np.sum(altSubjectData['BIC_NID']))[0])