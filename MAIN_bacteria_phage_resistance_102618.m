%% set parameters (lytic versus chronic)
K = 1; % carrying capacity
eta = [20,20]; % infection rate (temperate, chronic)
beta = [100,10]; % burst size/production rate (temperate, chronic)
gamma = 1; % growth rate
delta = 4; % production delay (rise and eclipse phase)
f = [0.01,0.01]; % latency rate (temperate, chronic)
kappa = 1; % death rate due to antibiotics
d = 1; % phage degradation rate
beta_max = 100; % max burst size under max stress (chronic only)
lambda = 1; % reproduction reduction factor when chronically infected

T = 24*60*0.0051; % period of antibiotic dosing
tend = 150; % end time of simulation

% initial condition (both phage, temperate only, chronic only)
IC = [1e-3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1e-7,1e-7,0.0];
%IC = [1e-3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1e-7,0.0,0.0];
%IC = [1e-3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1e-7,0.0];

% times at which antibiotic will be taken
stress_times = T:T:tend;
% uncomment below if no antibiotics given
stress_times = [];

% antibiotic susceptibility (all bacterial strains)
k = kappa*ones(length(IC)-3,1);
% antibiotic resistance scenarios
scenario = 0;
switch scenario
   case 0
      % no antibiotic resistance (leave alone) 
   case 1
      k(10) = 0; % LFT_T survives
   case 2
      k(11) = 0; % CFT_T survives
   case 3
      k(13) = 0; % CFT_F survives
   case 4
      k(12) = 0; % LFT_F survives
   case 5
      k([4 7 10 11]) = zeros(4,1); % LT, IFT_T, LFT_T, CFT_T survive
   case 6
      k([6 9 12]) = zeros(3,1); % LF, IFT_F, LFT_F survive
   case 7
      k([5 8 13]) = zeros(3,1); % CF, IFT_C, CFT_F survive
   case 8
      k = zeros(length(IC)-3,1); % all survive
   otherwise
      % default no antibiotic resistance (leave alone)
end

% integrate system
 [T,Y] = integrate_lyticVSchronic_resistance_091318(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,IC);

% integrate system with no chronic phage
%[T,Y] = integrate_lyticVSchronic_resistance_noPF_091918(tend,K,gamma,eta,delta,beta,beta_max,lambda,f,k,d,stress_times,IC);

%%
% important bacteria to plot
ind = [1 4 5 6 10 11 12 13];

figure('pos',[10 10 800 800])
subplot(2,1,1)
plot(T,Y(:,ind(1:7)),'LineWidth',2)
hold on
plot(T,Y(:,ind(end)),'--k','LineWidth',2)
ylim([0 0.8])
xlabel('time','FontSize',20)
ylabel('population size','FontSize',20)
legend('S','L_T','P_C','L_C','LCT_T','PCT_T','LCT_C','PCT_C','Location','bestoutside');


subplot(2,1,2)
plot(T,Y(:,14:15),'LineWidth',2)
ylim([0,36])
xlabel('time','FontSize',20)
ylabel('population size','FontSize',20)
legend('   V_T  ','   V_C  ','Location','bestoutside');

 
% %plot bacteria versus phage
% figure
% subplot(2,1,1)
% scatter(sum(Y(:,1:11),2),sum(Y(:,12:13),2),'LineWidth',2)
% set(gca,'xscale','log','yscale','log')
% ylabel('phage','FontSize',20)
% 
% %plot bacteria versus phage/bacteria (total)
% subplot(2,1,2)
% scatter(sum(Y(:,1:11),2),sum(Y(:,12:13),2)./sum(Y(:,1:11),2),'LineWidth',2)
% set(gca,'xscale','log','yscale','log')
% xlabel('bacteria','FontSize',20)
% ylabel('phage/bacteria','FontSize',20)

