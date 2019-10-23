% lytic versus chronic with no superinfection (bacteria affected by
% antibiotic) - catches integration warnings now
function [T,Y] = integrate_lyticVSchronic_resistance_011119(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,A,decay,IC)
% integrate system
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.05,'NonNegative',ones(16,1));
%[T,Y] = ode15s(@(t,y)oderhs(t,y,K,gamma,eta,delta,beta,f,k,d,stress_times),[0,tend],IC,options);
lastwarn('')    
    [T,Y] = ode15s(@(t,y)oderhs(t,y,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,A,decay),0:0.05:tend,IC,options);
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    T = NaN*ones(floor(tend),1); Y = NaN*ones(floor(tend),length(IC));
end
end

% Right hand side
function dY = oderhs(t,Y,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,A,decay)
% Hollings halfway point
h = 1/2;
% unpack infection rates
etaT = eta(1); etaF = eta(2);
% unpack latent fractions
fT = f(1); fF = f(2);
% unpack burst sizes
betaT = beta(1); betaF = beta(2);

% unpack vector into variable names (bacteria and phage)
S=Y(1); IT=Y(2); IF=Y(3); 
LT = Y(4); CF=Y(5); LF=Y(6); 
IFT_T=Y(7); IFT_C=Y(8); IFT_F=Y(9); 
LFT_T=Y(10); CFT_T=Y(11);
LFT_F = Y(12); CFT_F = Y(13);
PT=Y(14); PF=Y(15);

% new state: antibiotic resistant eclipse (ready to burst)
IR=Y(16);

% total bacteria and phage
Btot = S+IT+IF+LT+CF+LF+IFT_T+IFT_C+IFT_F+LFT_T+CFT_T+LFT_F+CFT_F + IR;

% indicator of antibiotic resistance
k_ind = (k>0);

% dynamical system (lytic versus chronic)
   % growth             - infection (temp)      - infection (fila)      - antibiotic death
dS = gamma*S*(1-Btot/K) - hollings(etaT,h,PT,S) - hollings(etaF,h,PF,S) - k(1)*s(t,stress_times,A,decay)*S;
    % infection                    - burst    + lysis due to stress                                             - antibiotic death
dIT = hollings(etaT*(1-fT),h,PT,S) - delta*IT + s(t,stress_times,A,decay)*(k_ind(4)*LT+k_ind(10)*LFT_T+k_ind(11)*CFT_T) - k(2)*s(t,stress_times,A,decay)*IT;
    % infection                    - eclipse  - antibiotic death
dIF = hollings(etaF*(1-fF),h,PF,S) - delta*IF - k(3)*s(t,stress_times,A,decay)*IF;
    % growth              + infection                - cross infection        - lysis                - antibiotic death
dLT = gamma*LT*(1-Btot/K) + hollings(etaT*fT,h,PT,S) - hollings(etaF,h,PF,LT) - s(t,stress_times,A,decay)*LT - k(4)*s(t,stress_times,A,decay)*LT;
    % growth                                                     + chronic  - cross infection        - antibiotic death
dCF = gamma_stress(lambda*gamma,s(t,stress_times,A,decay))*CF*(1-Btot/K) + delta*IF - hollings(etaT,h,PT,CF) - k(5)*s(t,stress_times,A,decay)*CF;
    % growth              + infection                - cross infection        - antibiotic death
dLF = gamma*LF*(1-Btot/K) + hollings(etaF*fF,h,PF,S) - hollings(etaT,h,PT,LF) - k(6)*s(t,stress_times,A,decay)*LF;
       % infection                     - eclipse     - antibiotic death
dIFT_T = hollings(etaF*(1-fF),h,PF,LT) - delta*IFT_T - k(7)*s(t,stress_times,A,decay)*IFT_T;
       % infection                     - burst       + induce lysis                      - antibiotic death
dIFT_C = hollings(etaT*(1-fT),h,PT,CF) - delta*IFT_C + s(t,stress_times,A,decay)*k_ind(13)*CFT_F - k(8)*s(t,stress_times,A,decay)*IFT_C;
       % infection                     - burst       + induce lysis                      - antibiotic death
dIFT_F = hollings(etaT*(1-fT),h,PT,LF) - delta*IFT_F + s(t,stress_times,A,decay)*k_ind(12)*LFT_F - k(9)*s(t,stress_times,A,decay)*IFT_F;
       % growth                 + cross infection          - induce lysis            - antibiotic death
dLFT_T = gamma*LFT_T*(1-Btot/K) + hollings(etaF*fF,h,PF,LT) - s(t,stress_times,A,decay)*LFT_T - k(10)*s(t,stress_times,A,decay)*LFT_T;
       % growth                                                        + chronic     - induce lysis            - antibiotic death
dCFT_T = gamma_stress(lambda*gamma,s(t,stress_times,A,decay))*CFT_T*(1-Btot/K) + delta*IFT_T - s(t,stress_times,A,decay)*CFT_T - k(11)*s(t,stress_times,A,decay)*CFT_T;
       % growth                 + cross infection           - induce lysis            - antibiotic death
dLFT_F = gamma*LFT_F*(1-Btot/K) + hollings(etaT*fT,h,PT,LF) - s(t,stress_times,A,decay)*LFT_F - k(12)*s(t,stress_times,A,decay)*LFT_F;
       % growth                                                        + cross infection           - induce lysis            - antibiotic death
dCFT_F = gamma_stress(lambda*gamma,s(t,stress_times,A,decay))*CFT_F*(1-Btot/K) + hollings(etaT*fT,h,PT,CF) - s(t,stress_times,A,decay)*CFT_F - k(13)*s(t,stress_times,A,decay)*CFT_F;
    % burst                           - absorbtion                    - degradation
dPT = betaT*delta*(IT+IFT_C+IFT_F+IR) - hollings(etaT,h,PT,(S+CF+LF)) - d*PT;
    % rise phase                                                     - absorbtion                 - degradation
dPF = beta_stress(betaF,beta_max,s(t,stress_times,A,decay))*(CF+CFT_T+CFT_F) - hollings(etaF,h,PF,(S+LT)) - d*PF;

    % - burst    + lysis due to stress (1-k_ind=0 if vulnerable to antibiotics)
dIR = - delta*IR + s(t,stress_times,A,decay)*((1-k_ind(4))*LT+(1-k_ind(10))*LFT_T+(1-k_ind(11))*CFT_T+(1-k_ind(12))*LFT_F+(1-k_ind(13))*CFT_F);

% pack variables into vector
dY = [dS;dIT;dIF;dLT;dCF;dLF;dIFT_T;dIFT_C;dIFT_F;dLFT_T;dCFT_T;dLFT_F;dCFT_F;dPT;dPF; dIR];

end

function eta_new = hollings(eta,h,P,B)
    % Option 1: mass action
    %eta_new = eta*P*B;
    % Option 2: bacteria is limiting factor (Michaelis-Menten)
    eta_new = eta*P*B/(h+B);
end

function beta_new = beta_stress(beta,beta_max,s)
    % Option 1: burst size independent of stress
    %beta_new = beta;
    % Option 2: linear response to stress
    %beta_new = beta + (beta_max-beta)*s;
    % Option 3: hollings type response to stress
    h = 1;
    beta_new = beta+s/(h+s)*(beta_max-beta);
end

function gamma_new = gamma_stress(gamma,s)
    % Option 1: growth rate independent of stress
    %gamma_new = gamma;
    % Option 2: linear response to stress
    %gamma_new = gamma*(1-s); 
    % Option 3: hollings type response to stress
    h = 1;
    gamma_new = gamma*(1-s/(h+s));
end

function stress = s(t,stress_times,A,decay)
    % t is current time, stress_times is vector of stress times
    %A = 1.1; k = 0.328;
    % Option 1: top hat function (stress on or off)
    %stress = peak*sum(heaviside(t-stress_times).*heaviside(stress_times+duration-t));
    % Option 2: spike with exponential decay
    stress = A*sum(heaviside(t-stress_times).*exp(-decay*(t-stress_times)));
    % threshold amount of stress in system
    %stress = min(stress, 1);
end

