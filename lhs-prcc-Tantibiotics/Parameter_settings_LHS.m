% PARAMETER BASELINE VALUES
% set parameters (lytic versus chronic)
K = 1; 
etaT = 20; etaC = 20;
betaT = 100; betaC = 10; 
gamma = 1; 
delta = 4; 
fT = 0.01; fC = 0.01; 
kappa = 1;
d = 1;
% max burst size under max stress (chronic only)
beta_max = 100;
% reproduction reduction factor when chronically infected
lambda = 1;
% antibiotic intensity
A = 1.1;
% antibiotic decay rate
decay = 0.328;

% "goal" average bacteria population
Bthresh = 0.1;

% ratio of PF to PT in IC
ratio = 1.0;
% period of antibiotic dosing
T0 = 10; % initial guess for Tstar

% Parameter Labels 
PRCC_var={'lambda','etaT','etaC','kappa','A','k','delta','fT','fC','betaT','betaC','bMax','d','ratioIC'}; 

% TIME SPAN OF THE SIMULATION
tend=300; % length of the simulations

% Variables Labels
y_var_label={'T'};
