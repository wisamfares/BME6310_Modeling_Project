load('TversusB_finalresults.mat')
freq = 0.001:0.01:0.501;
scenarios = 0:8;

%%
I = ~isnan(TversusB0);
figure
plot(freq(I(1,:)),TversusB0(1,I(1,:)),freq(I(2,:)),TversusB0(2,I(2,:)),...
     freq(I(3,:)),TversusB0(3,I(3,:)),freq(I(4,:)),TversusB0(4,I(4,:)),...
     freq(I(5,:)),TversusB0(5,I(5,:)),freq(I(6,:)),TversusB0(6,I(6,:)),...
     freq(I(7,:)),TversusB0(7,I(7,:)),'LineWidth',2)
hold on
plot(freq(I(8,:)),TversusB0(8,I(8,:)),'--',freq(I(9,:)),TversusB0(9,I(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend(num2str(scenarios'),'Location','bestoutside')
set(gca, 'YScale', 'log')

%%%%%
scenarios = [0 8];
I = ~isnan(TversusB1);
figure
plot(freq(I(1,:)),TversusB1(1,I(1,:)),freq(I(9,:)),TversusB1(9,I(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend(num2str(scenarios'),'Location','bestoutside')
set(gca, 'YScale', 'log')

%%%%%
scenarios = [0 5 8];
I = ~isnan(TversusB2);
figure
plot(freq(I(1,:)),TversusB2(1,I(1,:)),freq(I(6,:)),TversusB2(6,I(6,:)),...
     freq(I(9,:)),TversusB2(9,I(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend(num2str(scenarios'),'Location','bestoutside')
set(gca, 'YScale', 'log')

%%%%%
scenarios = [0 6 7 8];
I = ~isnan(TversusB3);
figure
plot(freq(I(1,:)),TversusB3(1,I(1,:)),freq(I(7,:)),TversusB3(7,I(7,:)),'LineWidth',2)
hold on
plot(freq(I(8,:)),TversusB3(8,I(8,:)),'--',freq(I(9,:)),TversusB3(9,I(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend(num2str(scenarios'),'Location','bestoutside')
set(gca, 'YScale', 'log')







%% Plot PF impact versus no phage impact (scenarios 0,8)
I0 = ~isnan(TversusB0); I1 = ~isnan(TversusB1); I2 = ~isnan(TversusB2); I3 = ~isnan(TversusB3);
figure
plot(freq(I3(1,:)),TversusB3(1,I3(1,:)),freq(I1(1,:)),TversusB1(1,I1(1,:)),...
     freq(I3(9,:)),TversusB3(9,I3(9,:)),freq(I1(9,:)),TversusB1(9,I1(9,:)),'LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend('P_F only - vulnerable','no phage - vulnerable','P_F only - resistant','no phage - resistant','Location','bestoutside')
set(gca, 'YScale', 'log')

% Plot PT impact versus no phage impact (scenarios 0,8)
I0 = ~isnan(TversusB0); I1 = ~isnan(TversusB1); I2 = ~isnan(TversusB2); I3 = ~isnan(TversusB3);
figure
plot(freq(I2(1,:)),TversusB2(1,I2(1,:)),freq(I1(1,:)),TversusB1(1,I1(1,:)),...
     freq(I2(9,:)),TversusB2(9,I2(9,:)),freq(I1(9,:)),TversusB1(9,I1(9,:)),'LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend('P_T only - vulnerable','no phage - vulnerable','P_T only - resistant','no phage - resistant','Location','bestoutside')
set(gca, 'YScale', 'log')

% Plot phage impact versus no phage impact (scenarios 0,8)
I0 = ~isnan(TversusB0); I1 = ~isnan(TversusB1); I2 = ~isnan(TversusB2); I3 = ~isnan(TversusB3);
figure
plot(freq(I0(1,:)),TversusB0(1,I0(1,:)),freq(I1(1,:)),TversusB1(1,I1(1,:)),...
     freq(I0(9,:)),TversusB0(9,I0(9,:)),freq(I1(9,:)),TversusB1(9,I1(9,:)),'LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend('both phage - vulnerable','no phage - vulnerable','both phage - resistant','no phage - resistant','Location','bestoutside')
set(gca, 'YScale', 'log')








%% Plot phage impact only (all antibiotic resistant)
scenarios = 8;
I0 = ~isnan(TversusB0); I1 = ~isnan(TversusB1); I2 = ~isnan(TversusB2); I3 = ~isnan(TversusB3);
figure
plot(freq(I0(scenarios+1,:)),TversusB0(scenarios+1,I0(scenarios+1,:)),...
     freq(I2(scenarios+1,:)),TversusB2(scenarios+1,I2(scenarios+1,:)),...
     freq(I3(scenarios+1,:)),TversusB3(scenarios+1,I3(scenarios+1,:)),...
     freq(I1(scenarios+1,:)),TversusB1(scenarios+1,I1(scenarios+1,:)),'LineWidth',2)
%hold on
%plot(freq(I2(9,:)),TversusB2(9,I2(9,:)),'--',freq(I3(9,:)),TversusB3(9,I3(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend('both phage','only P_T','only P_F','no phage','Location','bestoutside')
set(gca, 'YScale', 'log')

%% Plot phage impact only (all antibiotic vulnerable)
scenarios = 0;
I0 = ~isnan(TversusB0); I1 = ~isnan(TversusB1); I2 = ~isnan(TversusB2); I3 = ~isnan(TversusB3);
figure
plot(freq(I0(scenarios+1,:)),TversusB0(scenarios+1,I0(scenarios+1,:)),...
     freq(I2(scenarios+1,:)),TversusB2(scenarios+1,I2(scenarios+1,:)),...
     freq(I3(scenarios+1,:)),TversusB3(scenarios+1,I3(scenarios+1,:)),...
     freq(I1(scenarios+1,:)),TversusB1(scenarios+1,I1(scenarios+1,:)),'LineWidth',2)
%hold on
%plot(freq(I2(9,:)),TversusB2(9,I2(9,:)),'--',freq(I3(9,:)),TversusB3(9,I3(9,:)),'--','LineWidth',2)
xlabel('frequency (1/T) of antibiotic dosing','FontSize',16)
ylabel('average bacteria population size','FontSize',16)
xlim([0 0.3])
ylim([1e-2 1])
legend('both phage','only P_T','only P_F','no phage','Location','bestoutside')
set(gca, 'YScale', 'log')

