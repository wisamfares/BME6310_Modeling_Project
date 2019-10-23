% lytic versus chronic with no superinfection (bacteria affected by
% antibiotic) - catches integration warnings now
function [T,Y] = integrate_lyticVSchronic_resistance_noPF_091918(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,IC)
% integrate system
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.05,'NonNegative',ones(16,1));
%[T,Y] = ode15s(@(t,y)oderhs(t,y,K,gamma,eta,delta,beta,f,k,d,stress_times),[0,tend],IC,options);
lastwarn('')    
    [T,Y] = ode15s(@(t,y)oderhs(t,y,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times),0:0.05:tend,IC,options);
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    T = NaN*ones(floor(tend),1); Y = NaN*ones(floor(tend),length(IC));
end
end

% Right hand side
function dY = oderhs(t,Y,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times)
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
dS = gamma*S*(1-Btot/K) - hollings(etaT,h,PT,S) - hollings(etaF,h,PF,S) - k(1)*s(t,stress_times)*S;
    % infection                    - burst    + lysis due to stress                                             - antibiotic death
dIT = hollings(etaT*(1-fT),h,PT,S) - delta*IT + s(t,stress_times)*(k_ind(4)*LT+k_ind(10)*LFT_T+k_ind(11)*CFT_T) - k(2)*s(t,stress_times)*IT;
    % infection                    - eclipse  - antibiotic death
dIF = 0.0;
    % growth              + infection                - cross infection        - lysis                - antibiotic death
dLT = gamma*LT*(1-Btot/K) + hollings(etaT*fT,h,PT,S) - hollings(etaF,h,PF,LT) - s(t,stress_times)*LT - k(4)*s(t,stress_times)*LT;
    % growth                                                     + chronic  - cross infection        - antibiotic death
dCF = 0.0;
    % growth              + infection                - cross infection        - antibiotic death
dLF = 0.0;
       % infection                     - eclipse     - antibiotic death
dIFT_T = 0.0;
       % infection                     - burst       + induce lysis                      - antibiotic death
dIFT_C = 0.0;
       % infection                     - burst       + induce lysis                      - antibiotic death
dIFT_F = 0.0;
       % growth                 + cross infection          - induce lysis            - antibiotic death
dLFT_T = 0.0;
       % growth                                                        + chronic     - induce lysis            - antibiotic death
dCFT_T = 0.0;
       % growth                 + cross infection           - induce lysis            - antibiotic death
dLFT_F = 0.0;
       % growth                                                        + cross infection           - induce lysis            - antibiotic death
dCFT_F = 0.0;
    % burst                           - absorbtion                    - degradation
dPT = betaT*delta*(IT+IFT_C+IFT_F+IR) - hollings(etaT,h,PT,(S+CF+LF)) - d*PT;
    % rise phase                                                     - absorbtion                 - degradation
dPF = 0.0;

    % - burst    + lysis due to stress (1-k_ind=0 if vulnerable to antibiotics)
dIR = 0.0;

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

function stress = s(t,stress_times)
    % t is current time, stress_times is vector of stress times
    duration = 20; A = 1.1; k = 0.328;% 0.1;
    % Option 1: top hat function (stress on or off)
    %stress = peak*sum(heaviside(t-stress_times).*heaviside(stress_times+duration-t));
    % Option 2: spike with exponential decay
    stress = A*sum(heaviside(t-stress_times).*exp(-k*(t-stress_times)));
    % threshold amount of stress in system
    %stress = min(stress, 1);
end

