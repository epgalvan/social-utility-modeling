%% %%%%%%%% 1_4_0

% Create Construct Value Functions

function d_Inequality = inequality(a1, b1, a2, b2)
d_Inequality = abs(a2 - b2) - abs(a1 - b1);
end

function d_lossAdvantaged = harm(a0, b0, a1, b1, a2, b2)
initial = [a0, b0];
choice1 = [a1, b1];
choice2 = [a2, b2];
if a0 == b0
    advantaged = 1;
else
    [~, advantaged] = max(initial);
end
d_lossAdvantaged = (initial(advantaged) - choice1(advantaged)) - ...
    (initial(advantaged) - choice2(advantaged));
end

function rankReverseDiff = rankReverse(a0, b0, a1, b1, a2, b2)
if a0 == b0
    rankReverseDiff = 0;
    return;
end
d_initial = a0 - b0;
d_choice1 = a1 - b1;
d_choice2 = a2 - b2;
choice1Reversed = 0;
choice2Reversed = 0;
if d_initial > 0
    if d_choice1 < 0, choice1Reversed = 1; end
    if d_choice2 < 0, choice2Reversed = 1; end
elseif d_initial < 0
    if d_choice1 > 0, choice1Reversed = 1; end
    if d_choice2 > 0, choice2Reversed = 1; end
end
rankReverseDiff = choice1Reversed - choice2Reversed;
end


%% %%%%%%%% 1_5_0

% Preallocating and Defining Functions, TrialList, and Parameters

a0 = randi([10,20], 100, 1);
b0 = 20 - a0;
trialList = table(a0, b0);
clear a0 b0

for i = 1:length(trialList.a0)
    trialList.a1(i) = randi([5, trialList.a0(i)]);
end
trialList.b1 = 20 - trialList.a1;

for i = 1:length(trialList.a0)
    trialList.a2(i) = randi([5, trialList.a0(i)]);
    if trialList.a2(i) == trialList.a1(i)
        trialList.a2(i) = 10;
    end
end
trialList.b2 = 20 - trialList.a2;
trialList = [trialList; trialList(:, [2, 1, 4, 3, 6, 5])];

function util = utility(pars, IVs)
IVs = double(IVs);
a0 = IVs(1); b0 = IVs(2); a1 = IVs(3); b1 = IVs(4); a2 = IVs(5); b2 = IVs(6);
alpha = pars(1); delta = pars(2); rho = pars(3);
util = (alpha * inequality(a1, b1, a2, b2)) - ...
    (delta * harm(a0, b0, a1, b1, a2, b2)) - ...
    (rho * rankReverse(a0, b0, a1, b1, a2, b2));
end

function prob = probability(pars, utilitydiff)
beta = pars(end-2); epsilon = pars(end-1); gamma = pars(end);
prob = 1 / (1 + exp(-(beta * utilitydiff)));
prob = prob * (1 - 2 * epsilon) + epsilon + gamma * (2 * epsilon);
prob = max(min(prob, 0.9999999999), 0.00000000001);
end

freeParameters = struct();
mainpars = (0:3)./2;
bet = (0:3)./2;
eps = (0:3)./2;
gam = (-0.5:0.1:0.4);
for i = 1:length(mainpars)
    for j = 1:length(mainpars)
        for k = 1:length(mainpars)
            for l = 1:length(bet)
                for m = 1:length(eps)
                    for n = 1:length(gam)
                        freeParameters(i, j, k, l, m, n).alpha = mainpars(i) + rand(1,1)*0.5;
                        freeParameters(i, j, k, l, m, n).delta = mainpars(j) + rand(1,1)*0.5;
                        freeParameters(i, j, k, l, m, n).rho = mainpars(k) + rand(1,1)*0.;
                        freeParameters(i, j, k, l, m, n).beta = bet(l);
                        freeParameters(i, j, k, l, m, n).epsilon = eps(m);
                        freeParameters(i, j, k, l, m, n).gamma = gam(n) + rand(1, 1)*0.1;
                    end
                end    
            end
        end
    end    
end 

% Determine Predictions

function pred = generatePredictions(parameters, df)
pred = zeros(size(df, 1), 1);
for i = 1:size(df, 1)
    thisTrialIVs = double(df(i, :));
    utilityDiff = utility(parameters, thisTrialIVs);
    pred(i) = probability(parameters, utilityDiff);
end
end

%% %%%%%%%% 1_6_0

%% %%%%%%%% 1_7_0

%% %%%%%%%% 2_1_0

%% %%%%%%%% 2_2_0

%% %%%%%%%% 2_3_0

%% %%%%%%%% 2_4_0

%% %%%%%%%% 3_1_0

%% %%%%%%%% 3_2_0

%% %%%%%%%% 3_3_0