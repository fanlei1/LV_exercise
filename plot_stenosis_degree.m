close all
clear all
%Stenosis 100%
Qpas = [3.8444,2.9886,2.5028,2.2808];
Qreg = [0.9216,1.42,1.60,2.09];
%% Stenosis 80%
Qpas1 = [2.3627,1.6454,1.4979];
Qreg1 = [1.1858,1.2188,1.4841];
%% Stenosis 60%
Qpas4 = [1.3066,0.4019,0.3316];
Qreg4 = [1.0939,0.3509,0.3193];
x1 = [75,100,120,140];
x2 = [85,110,140];
x3 = x2;

figure(1)
plot(x1,Qpas,'*--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1,'<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x3,Qpas4,'s--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1,Qreg,'*:','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qreg1,'<:','LineWidth', 3, 'markersize', 14);
hold on
plot(x3,Qreg4,'s:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$\bar{Q}$ (ml/min/g)','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
xlim([70 145])
box off

%% stenosis1
figure(2)
plot(x1,Qpas./Qreg,'*--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1./Qreg1,'<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x3,Qpas4./Qreg4,'s--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1,1.5,'k:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('CFR','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
xlim([70 145])
box off

%% Stenosis 100%
PVA = [8105*75 11128*100 12815*120 12933*140];
Fmeta = [0.0636,0.4499,0.6959,0.9789];
%% Stenosis 80%
PVA1 = [7516*85,10839*110,13508*140];
Fmeta1 = [0.0915,0.5031,1];
%% Stenosis 60%
PVA4 = [8023*95,9988*120,13088*130];
Fmeta4 = [0.1944,0.5154,0.9897];
figure(3)
plot(x1,Fmeta,'*--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Fmeta1,'<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x3,Fmeta4,'s--','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('FFR 1.0','FFR 0.8','FFR 0.6','FontSize',22,'FontWeight','bold','Interpreter','latex','NumColumns',1);
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$F_{meta}$','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
xlim([70 145])
ylim([0 1.1])
box off

%% Stenosis 100%
Ees = [4.0025, 4.4756, 4.7396, 5.3274];
%% Stenosis 80%
Ees1 = [4.2225,4.4678,4.9736];
%% Stenosis 60%
Ees4 = [4.1678,4.2629,4.2733];%,4.5048];
figure(4)
plot(x1,Ees,'*--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Ees1,'<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x3,Ees4,'s--','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$E_{es}$','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
xlim([70 145])
box off