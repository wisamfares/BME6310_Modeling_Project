function Bave = averageB(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,T,A,decay,IC)
stress_times = T:T:tend;
% integrate system (insert date for different versions of model)
[~,Y] = integrate_lyticVSchronic_resistance_011119(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,A,decay,IC);
%T
Bave = mean(sum(Y(:,[1:13 16]),2),1);
end

