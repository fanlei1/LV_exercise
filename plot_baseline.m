close all
clear all
clc

figure(1)
[ha, pos] = tight_subplot(3,1,[.05 .0],[.13 .03],[.13 .02]);
axes(ha(1));
a = 10;
b = 20;  
BCL1 = 800;
BCL2 = 600;
BCL3 = 536;
BCL4 = 430;
dt = 0.5;
endpoint=BCL1*a+BCL2*b+BCL3*b+BCL4*b;
ts = (0:dt:endpoint-0.5);
load baseline_reg.mat
plot(ts,PVQ(:,2),'-','LineWidth', 2)
hold on
plot(ts,PVQ(:,3),'-','LineWidth', 2)
h=legend('$P_{LV}$','$P_{art}$','FontSize',18,'FontWeight','bold','Interpreter','latex','NumColumns',2);
ylabel('P (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(2));
plot(ts,PVQ(:,4),'-','LineWidth', 2)
ylabel('IMP (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(3));
load baseline_pas.mat
plot(ts,PVQ(:,5)*60,'-','LineWidth', 2)
load baseline_reg.mat
hold on
plot(ts,PVQ(:,5)*60,'-','LineWidth', 2)
h=legend('$Q_{pas}$','$Q_{reg}$','FontSize',18,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('Time (ms)','FontSize',18,'Interpreter','latex')
ylabel('Q (ml/min)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 endpoint])
set(gca,'FontSize',18)
box off

figure(2)
[ha, pos] = tight_subplot(2,4,[.0 .0],[.13 .03],[.13 .02]);
load baseline_reg.mat
axes(ha(1));
t1 = [BCL1/dt*(a-1)+1:1:BCL1/dt*a];
plot(ts(1:1:BCL1/dt),PVQ(t1,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL1/dt),PVQ(t1,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL1/dt),PVQ(t1,4),'-','LineWidth', 2)
h=legend('$P_{LV}$','$P_{art}$','IMP','FontSize',18,'FontWeight','bold','Interpreter','latex');
ylabel('P (mmHg)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 BCL1])
ylim([0 250])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(2));
t2 = [a*BCL1/dt+(b-2)*BCL2/dt+1+350/dt:1:a*BCL1/dt+(b-1)*BCL2/dt+350/dt];
plot(ts(1:1:BCL2/dt),PVQ(t2,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,4),'-','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL1])
ylim([0 250])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(3));
t3 = [a*BCL1/dt+b*BCL2/dt+(b-2)*BCL3/dt+340/dt+1:1:a*BCL1/dt+b*BCL2/dt+(b-1)*BCL3/dt+340/dt];
plot(ts(1:1:BCL3/dt),PVQ(t3,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3/dt),PVQ(t3,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3/dt),PVQ(t3,4),'-','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL1])
ylim([0 250])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(4));
t4 = [endpoint/dt-BCL4*2/dt+201/dt-1:1:endpoint/dt-BCL4/dt+200/dt];
plot(ts(1:1:BCL4/dt),PVQ(t4,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL4/dt),PVQ(t4,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL4/dt),PVQ(t4,4),'-','LineWidth', 2)
xlim([0 BCL1])
ylim([0 250])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(5));
load baseline_pas.mat
plot(ts(1:10:BCL1/dt),PVQ(t1(1):10:t1(end),5)*60,'-','LineWidth', 2)
load baseline_reg.mat
hold on
plot(ts(1:10:BCL1/dt),PVQ(t1(1):10:t1(end),5)*60,'--','LineWidth', 2)
h=legend('$Q_{pas}$','$Q_{reg}$','FontSize',18,'FontWeight','bold','Interpreter','latex');
ylabel('Q (ml/min)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 BCL1])
ylim([-300 500])
set(gca,'FontSize',18)
box off

axes(ha(6));
load baseline_pas.mat
plot(ts(1:10:BCL2/dt),PVQ(t2(1):10:t2(end),5)*60,'-','LineWidth', 2)
load baseline_reg.mat
hold on
plot(ts(1:10:BCL2/dt),PVQ(t2(1):10:t2(end),5)*60,'--','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL1])
ylim([-300 500])
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

axes(ha(7));
load baseline_pas.mat
plot(ts(1:16:BCL3/dt),PVQ(t3(1):16:t3(end),5)*60,'-','LineWidth', 2)
load baseline_reg.mat
hold on
plot(ts(1:16:BCL3/dt),PVQ(t3(1):16:t3(end),5)*60,'--','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL1])
ylim([-300 500])
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

axes(ha(8));
load baseline_pas.mat
plot(ts(1:10:BCL4/dt),PVQ(t4(1):10:t4(end),5)*60,'-','LineWidth', 2)
load baseline_reg.mat
hold on
plot(ts(1:10:BCL4/dt),PVQ(t4(1):10:t4(end),5)*60,'--','LineWidth', 2)
xlabel('Time (ms)','FontSize',18,'Interpreter','latex')
set(h,'Box','off');
xlim([0 BCL1])
ylim([-300 500])
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

figure(3)
% PVA_(1) = polyarea(PVQ(t1,1),PVQ(t1,2))+29*116/2;
% PVA_(2) = polyarea(PVQ(t2,1),PVQ(t2,2))+31*133/2;
% PVA_(3) = polyarea(PVQ(t3,1),PVQ(t3,2))+32*142/2;
% PVA_(4) = polyarea(PVQ(t4,1),PVQ(t4,2))+33*175/2;
plot(PVQ(t1,1),PVQ(t1,2),'LineWidth', 2)  %/1e3/60
hold on
plot(PVQ(t2,1),PVQ(t2,2),'LineWidth', 2)
hold on
plot(PVQ(t3,1),PVQ(t3,2),'LineWidth', 2)
hold on
plot(PVQ(t4,1),PVQ(t4,2),'LineWidth', 2)
h=legend('75 bpm','100 bpm','120 bpm','140 bpm','FontSize',28,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('$V_{LV}$ (ml)','FontSize',28,'Interpreter','latex')
ylabel('$P_{LV}$ (mmHg)','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
xlim([30 140])
ylim([0 280])
set(gca,'FontSize',28)
box off

figure(4)
x = categorical({'75','100','120','140'});
x = reordercats(x,{'75','100','120','140'});
subplot(1,2,1)
Ees = [4.0025, 4.4756, 4.7396, 5.3274];
plot(x, Ees, 'k+:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$E_{es}$ (mmHg/ml)','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
box off

subplot(1,2,2)
PVA = [8105*75 11128*100 12815*120 12933*140];
Fmeta = [0.0636,0.4499,0.6959,0.9789];
plot(PVA/1000, Fmeta, 'k*--','LineWidth', 3, 'markersize', 14);
xlabel('PVA$\cdot$HR (mmHg$\cdot$L$\cdot$bpm)','FontSize',24,'Interpreter','latex')
ylabel('$F_{meta}$','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(5)
Qpas = [3.8444,2.9886,2.5028,2.2808];
Qreg = [0.9216,1.42,1.60,2.09];
Qregdata = [0.9,1.332,1.494,1.94];
plot(x,Qpas,'k<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x,Qreg,'ks:','LineWidth', 3, 'markersize', 14);
hold on
plot(x,Qregdata,'r*','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('$\bar{Q}_{pas}$','$\bar{Q}_{reg}$','$\bar{Q}_{reg,data}$','FontSize',28,'FontWeight','bold','Interpreter','latex');
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$\bar{Q}$ (ml/min/g)','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(6)
plot(x,Qpas./Qreg,'kd--','LineWidth', 3, 'markersize', 14);
hold on
plot(x,Qreg(end)./Qreg,'k>:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('Adenosine','Stress test','FontSize',28,'FontWeight','bold','Interpreter','latex');
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('CFR','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off