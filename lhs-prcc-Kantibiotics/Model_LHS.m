% The results should be compared to the PRCC results section in
% Supplementary Material D and Table D.1 for different N (specified by
% "runs" in the script below
clear all;
close all;

rng(1);
% Sample size N
runs=150;

% LHS MATRIX 
Parameter_settings_LHS;

% LHS_Call(xmin,xmean,xmax,xsd,nsample,distrib,threshold) 
lambda_LHS = LHS_Call(0.9, lambda, 1.1, 0 ,runs,'unif'); 
etaT_LHS = LHS_Call(19, etaT, 21, 0 ,runs,'unif'); 
etaC_LHS = LHS_Call(19, etaC, 21, 0 ,runs,'unif'); 
T_LHS = LHS_Call(9, T, 11, 0 ,runs,'unif');
A_LHS = LHS_Call(1.09, A, 1.11, 0, runs,'unif'); 
k_LHS = LHS_Call(0.32, decay, 0.34, 0 ,runs,'unif'); 
delta_LHS = LHS_Call(3.5, delta, 4.5, 0 ,runs,'unif'); 
fT_LHS = LHS_Call(0.009, fT, 0.011, 0 ,runs,'unif'); 
fC_LHS = LHS_Call(0.009, fC, 0.011, 0, runs,'unif'); 
betaT_LHS = LHS_Call(95, betaT, 105, 0 ,runs,'unif');
betaC_LHS = LHS_Call(5, betaC, 15, 0 ,runs,'unif');
betaMax_LHS = LHS_Call(95, beta_max, 105, 0 ,runs,'unif');
d_LHS = LHS_Call(0.9, d, 1.1, 0 ,runs,'unif');
ratio_LHS = LHS_Call(0.9, ratio, 1.1, 0 ,runs,'unif');

% LHS MATRIX and PARAMETER LABELS
LHSmatrix = [lambda_LHS etaT_LHS etaC_LHS T_LHS A_LHS k_LHS ...
           delta_LHS fT_LHS fC_LHS betaT_LHS betaC_LHS betaMax_LHS d_LHS ratio_LHS];
% options for fminsearch
options = optimset('TolFun',1e-3,'TolX',1e-2);
% output that needs to be measured
Kstar_lhs = zeros(1,runs);
tic
for x=1:runs %Run solution x times choosing different values
    x
    % initial condition
    IC = [1e-3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1e-7,LHSmatrix(x,14)*1e-7,0.0];
    Bdiff = @(kappaval) (averageB(tend,K,gamma,[LHSmatrix(x,2),LHSmatrix(x,3)],...
                  LHSmatrix(x,7),[LHSmatrix(x,10),LHSmatrix(x,11)],LHSmatrix(x,12),...
                  LHSmatrix(x,1),[LHSmatrix(x,8),LHSmatrix(x,9)],kappaval*ones(length(IC)-3,1),...
                  LHSmatrix(x,13),LHSmatrix(x,4),LHSmatrix(x,5),LHSmatrix(x,6),IC)-Bthresh)^2;
    
    % save averages
    Kstar_lhs(1,x) = fminsearch(Bdiff,kappa0,options);
    
end
toc
% Save the workspace
save Model_LHS.mat;
%CALCULATE PRCC
alpha=0.05;
[prccB, signB, sign_labelB]=PRCC(LHSmatrix,Kstar_lhs,1,PRCC_var,alpha);
