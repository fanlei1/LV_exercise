close all
clear all
clc

figure(1)
[ha, pos] = tight_subplot(3,1,[.05 .0],[.13 .03],[.13 .02]);
axes(ha(1));
a = 10;
a1 = 20;
b = 20;  
BCL0 = 800;
BCL1 = 700;
BCL2 = 550;
BCL3 = 429;
BCL4 = 375;
dt = 0.5;
deltat=0.001;
endp = BCL1*a+BCL2*a1+BCL3*b;
tend = (0:dt:endp-1);
endpoint=BCL1*a+BCL2*a1+BCL3*b;%+BCL4*b;
ts = (0:dt:endpoint-1);
load stenosisFFR8_reg.mat
plot(ts,PVQ(1:1:endpoint,2),'-','LineWidth', 2)
hold on
plot(ts,PVQ(1:1:endpoint,3),'-','LineWidth', 2)
h=legend('$P_{LV}$','$P_{art}$','FontSize',18,'FontWeight','bold','Interpreter','latex','NumColumns',2);
ylabel('P (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(2));
plot(ts,PVQ(1:1:endpoint,4),'-','LineWidth', 2)
ylabel('IMP (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(3));
load stenosisFFR8_pas.mat
plot(ts,PVQ(1:1:endpoint,5)*60,'-','LineWidth', 2)
load stenosisFFR8_reg.mat
hold on
plot(ts,PVQ(1:1:endpoint,5)*60,'-','LineWidth', 2)
h=legend('$Q_{pas}$','$Q_{reg}$','FontSize',18,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('Time (ms)','FontSize',18,'Interpreter','latex')
ylabel('Q (ml/min)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
box off

figure(2)
[ha, pos] = tight_subplot(2,3,[.0 .0],[.13 .03],[.13 .02]);
load stenosisFFR8_reg.mat
axes(ha(1));
t1 = [BCL1/dt*(a-1)-5:1:BCL1/dt*a-6];
plot(ts(1:1:BCL1/dt),PVQ(t1,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL1/dt),PVQ(t1,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL1/dt),PVQ(t1,4),'-','LineWidth', 2)
h=legend('$P_{LV}$','$P_{art}$','IMP','FontSize',18,'FontWeight','bold','Interpreter','latex');
ylabel('P (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 BCL0])
ylim([0 260])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(2));
t2 = [a*BCL1/dt+(a1-2)*BCL2/dt+1+150:1:a*BCL1/dt+(a1-1)*BCL2/dt+150];
plot(ts(1:1:BCL2/dt),PVQ(t2,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,4),'-','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([0 260])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(3));
t3 = [a*BCL1/dt+a1*BCL2/dt+(b-2)*BCL3/dt+1+0:1:a*BCL1/dt+a1*BCL2/dt+(b-1)*BCL3/dt];
plot(ts(1:1:BCL3),PVQ(t3,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3),PVQ(t3,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3),PVQ(t3,4),'-','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([0 260])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(4));
load stenosisFFR8_pas.mat
plot(ts(1:1:BCL1/dt),PVQ(t1,5)*60,'-','LineWidth', 2)
load stenosisFFR8_reg.mat
hold on
plot(ts(1:1:BCL1/dt),PVQ(t1,5)*60,'--','LineWidth', 2)
h=legend('$Q_{pas}$','$Q_{reg}$','FontSize',18,'FontWeight','bold','Interpreter','latex');
ylabel('Q (ml/min)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 BCL0])
ylim([-300 500])
set(gca,'FontSize',18)
box off

axes(ha(5));
load stenosisFFR8_pas.mat
plot(ts(1:1:BCL2/dt),PVQ(t2,5)*60,'-','LineWidth', 2)
load stenosisFFR8_reg.mat
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,5)*60,'--','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([-300 500])
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

axes(ha(6));
load stenosisFFR8_pas.mat
plot(ts(1:1:BCL3/dt),PVQ(t3,5)*60,'-','LineWidth', 2)
load stenosisFFR8_reg.mat
hold on
plot(ts(1:1:BCL3/dt),PVQ(t3,5)*60,'--','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([-300 500])
xlabel('Time (ms)','FontSize',18,'Interpreter','latex')
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

figure(3)
plot(PVQ(t1,1),PVQ(t1,2),'LineWidth', 2)  %/1e3/60
hold on
plot(PVQ(t2,1),PVQ(t2,2),'LineWidth', 2)
hold on
plot(PVQ(t3,1),PVQ(t3,2),'LineWidth', 2)
h=legend('85 bpm','110 bpm','140 bpm','FontSize',28,'FontWeight','bold','Interpreter','latex','NumColumns',1);
xlabel('$V_{LV}$ (ml)','FontSize',28,'Interpreter','latex')
ylabel('$P_{LV}$ (mmHg)','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(4)
x1 = [75,100,120,140];
x2 = [85,110,140];
subplot(1,2,1)
Ees = [4.0025, 4.4756, 4.7396, 5.3274];
%% Stenosis 90%
%Ees1 = [4.4796,4.6991,5.2094,5.2027];
%% Stenosis 80%
Ees1 = [4.2225,4.4678,4.9736];
%% Stenosis 70%
% Ees1 = [4.3308,4.4101,4.6701,4.8592];
%% Stenosis 60%
% Ees1 = [4.1142,4.1471,4.3604,4.5048];
plot(x1, Ees, 'k*:','LineWidth', 3, 'markersize', 18);
hold on
plot(x2, Ees1, 'd:','LineWidth', 3, 'markersize', 18);
h=legend('base','FFR 0.8','FontSize',28,'FontWeight','bold','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$E_{es}$ (mmHg/ml)','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

subplot(1,2,2)
PVA = [8105*75 11128*100 12815*120 12933*140];
Fmeta = [0.0636,0.4499,0.6959,0.9789];
plot(PVA/1000, Fmeta, 'k*--','LineWidth', 3, 'markersize', 18);
PVA1 = [7516*85,10839*110,13508*140];
Fmeta1 = [0.0915,0.5031,1];
hold on
plot(PVA1/1000, Fmeta1, 'd--','LineWidth', 3, 'markersize', 18);
xlabel('PVA$\cdot$HR (mmHg$\cdot$L$\cdot$bpm)','FontSize',28,'Interpreter','latex')
ylabel('$F_{meta}$','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
box off

figure(5)
Qpas = [3.8444,2.9886,2.5028,2.2808];
Qreg = [0.9216,1.42,1.60,2.09];
Qpas1 = [2.3627,1.6454,1.4979];
Qreg1 = [1.1858,1.2188,1.4841];
plot(x1,Qpas,'k<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1,Qreg,'ks:','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1,'r<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qreg1,'rs:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('$\bar{Q}_{pas}$,base','$\bar{Q}_{reg}$,base','$\bar{Q}_{pas}$,FFR 0.8','$\bar{Q}_{reg}$,FFR 0.8','FontSize',28,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$\bar{Q}$ (ml/min/g)','FontSize',28,'Interpreter','latex')
xlim([70 145])
ylim([0.8 4.7])
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(6)
plot(x1,Qpas./Qreg,'kd--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1./Qreg1,'rd--','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('Aden,base','Aden,FFR 0.8','FontSize',28,'FontWeight','bold','Interpreter','latex');
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('CFR','FontSize',28,'Interpreter','latex')
xlim([70 145])
set(h,'Box','off');
set(gca,'FontSize',28)
box off