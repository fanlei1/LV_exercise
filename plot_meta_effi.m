close all
clear all
clc

figure(1)
[ha, pos] = tight_subplot(3,1,[.05 .0],[.13 .03],[.13 .02]);
axes(ha(1));
a = 10;
a1 = 20;
b = 20;  
BCL1 = 800;
BCL2 = 600;
BCL3 = 536;
BCL4 = 430;
BCL0 = 800;
dt = 0.5;
endp = BCL1*a+BCL2*a1+BCL3*b;
tend = (0:dt:endp-1);
endpoint=BCL1*a+BCL2*a1+BCL3*b;
ts = (0:dt:endpoint-1);
load metabolic_reg.mat
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
load metabolic_pas.mat
plot(ts,PVQ(1:1:endpoint,5)*60,'-','LineWidth', 2)
load metabolic_reg.mat
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
load metabolic_reg.mat
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
xlim([0 BCL0])
ylim([0 260])
set(gca,'FontSize',18)
set(gca,'xtick',[])
box off

axes(ha(2));
t2 = [a*BCL1/dt+(a1-2)*BCL2/dt+1-200:1:a*BCL1/dt+(a1-1)*BCL2/dt-200];
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
t3 = [a*BCL1/dt+a1*BCL2/dt+(b-2)*BCL3/dt+1-200:1:a*BCL1/dt+a1*BCL2/dt+(b-1)*BCL3/dt-200];
plot(ts(1:1:BCL3/dt),PVQ(t3,2),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3/dt),PVQ(t3,3),'-','LineWidth', 2)
hold on
plot(ts(1:1:BCL3/dt),PVQ(t3,4),'-','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([0 260])
set(gca,'FontSize',18)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off

axes(ha(4));
load metabolic_pas.mat
plot(ts(1:1:BCL1/dt),PVQ(t1,5)*60,'-','LineWidth', 2)
load metabolic_reg.mat
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
load metabolic_pas.mat
plot(ts(1:1:BCL2/dt),PVQ(t2,5)*60,'-','LineWidth', 2)
load metabolic_reg.mat
hold on
plot(ts(1:1:BCL2/dt),PVQ(t2,5)*60,'--','LineWidth', 2)
set(h,'Box','off');
xlim([0 BCL0])
ylim([-300 500])
set(gca,'FontSize',18)
set(gca,'ytick',[])
box off

axes(ha(6));
load metabolic_pas.mat
plot(ts(1:1:BCL3/dt),PVQ(t3,5)*60,'-','LineWidth', 2)
load metabolic_reg.mat
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
plot(PVQ(t1,1),PVQ(t1,2),'LineWidth', 2)  
hold on
plot(PVQ(t2,1),PVQ(t2,2),'LineWidth', 2)
hold on
plot(PVQ(t3,1),PVQ(t3,2),'LineWidth', 2)
h=legend('100 bpm','120 bpm','140 bpm','FontSize',28,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('$V_{LV}$ (ml)','FontSize',28,'Interpreter','latex')
ylabel('$P_{LV}$ (mmHg)','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(4)
x1 = categorical({'75','100','120','140'});
x1 = reordercats(x1,{'75','100','120','140'});
x2 = categorical({'75','100','120'});
x2 = reordercats(x2,{'75','100','120'});
subplot(1,2,1)
Ees = [4.0025, 4.4756, 4.7396, 5.3274];
Ees1 = [4.0257,4.6829,5.0612];
plot(x1, Ees, 'k*:','LineWidth', 3, 'markersize', 18);
hold on
plot(x2, Ees1, 'd:','LineWidth', 3, 'markersize', 18);
h=legend('$k_2$','$k_2$+','FontSize',28,'FontWeight','bold','Interpreter','latex');
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
PVA1 = [8126*75,11458*100,13298*120];
Fmeta1 = [0.0985,0.6176,0.99];
hold on
plot(PVA1/1000, Fmeta1, 'd--','LineWidth', 3, 'markersize', 18);
xlabel('PVA (mmHg$\cdot$L)','FontSize',28,'Interpreter','latex')
ylabel('$F_{meta}$','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
box off

figure(5)
Qpas = [3.8444,2.9886,2.5028,2.2808];
Qreg = [0.9216,1.42,1.60,2.09];
Qpas1 = [3.8476,2.9914,2.4904];
Qreg1 = [0.9940,1.8539,2.3732];
plot(x1,Qpas,'k<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1,Qreg,'ks:','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1,'r<--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qreg1,'rs:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('$\bar{Q}_{pas}$,$k_2$','$\bar{Q}_{reg}$,$k_2$','$\bar{Q}_{pas}$,$k_2$+','$\bar{Q}_{reg}$,$k_2$+','FontSize',28,'FontWeight','bold','Interpreter','latex','NumColumns',2);
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$\bar{Q}$ (ml/min/g)','FontSize',28,'Interpreter','latex')
ylim([0.9 4.7])
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(6)
plot(x1,Qpas./Qreg,'kd--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1./Qreg1,'rd--','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
h=legend('Aden,$k_2$','Aden,$k_2$+','FontSize',28,'FontWeight','bold','Interpreter','latex');
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('CFR','FontSize',28,'Interpreter','latex')
set(h,'Box','off');
set(gca,'FontSize',28)
box off

figure(7)
x1 = [75, 100, 120, 140]; 
x2 = [85, 110, 140]; 
x3 = [75, 100, 120];
subplot(1,2,1)
Ees = [4.0025, 4.4756, 4.7396, 5.3274];
Ees1 = [4.2911,4.5577, 5.3113];
Ees2 = [4.0257,4.6829,5.0612];
plot(x1, Ees, 'k*:','LineWidth', 3, 'markersize', 18);
hold on
plot(x2, Ees1, 'o--','LineWidth', 3, 'markersize', 18);
hold on
plot(x1(1:1:3), Ees2, 'd-','LineWidth', 3, 'markersize', 18);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$E_{es}$ (mmHg/ml)','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
box off

subplot(1,2,2)
PVA = [8105*75 11128*100 12815*120 12933*140];
Fmeta = [0.0636,0.4499,0.6959,0.9789];
plot(PVA/1000, Fmeta, 'k*:','LineWidth', 3, 'markersize', 18);
PVA1 = [9375*85,12791*110,14595*140];
Fmeta1 = [0.2372,0.6720,1];
PVA2 = [8126*75,11458*100,13298*120];
Fmeta2 = [0.0985,0.6176,0.99];
hold on
plot(PVA1/1000, Fmeta1, 'o--','LineWidth', 3, 'markersize', 18);
hold on
plot(PVA2/1000, Fmeta2, 'd-','LineWidth', 3, 'markersize', 18);
xlabel('PVA$\cdot$HR (mmHg$\cdot$L$\cdot$bpm)','FontSize',28,'Interpreter','latex')
ylabel('$F_{meta}$','FontSize',28,'Interpreter','latex')
set(gca,'FontSize',28)
box off

figure(8)
Qpas = [3.8444,2.9886,2.5028,2.2808];
Qreg = [0.9216,1.42,1.60,2.09];
Qpas1 = [2.9560,2.7547,2.2083];
Qreg1 = [1.6548,1.7470,2.190];
Qpas2 = [3.8476,2.9914,2.4904];
Qreg2 = [0.9940,1.8539,2.3732];
plot(x1,Qpas,'k*--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1,Qreg,'k*:','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1,'ro--','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qreg1,'ro:','LineWidth', 3, 'markersize', 14);
hold on
plot(x1(1:1:3),Qpas2,'bd--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1(1:1:3),Qreg2,'bd:','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('$\bar{Q}$ (ml/min/g)','FontSize',28,'Interpreter','latex')
xlim([70 145])
set(gca,'FontSize',28)
box off

figure(9)
plot(x1,Qpas./Qreg,'k*:','LineWidth', 3, 'markersize', 14);
hold on
plot(x2,Qpas1./Qreg1,'o--','LineWidth', 3, 'markersize', 14);
hold on
plot(x1(1:1:3),Qpas2./Qreg2,'d-','LineWidth', 3, 'markersize', 14);
set(gca,'TickLabelInterpreter','latex')
xlabel('HR (bpm)','FontSize',28,'Interpreter','latex')
ylabel('CFR','FontSize',28,'Interpreter','latex')
xlim([70 145])
set(gca,'FontSize',28)
box off