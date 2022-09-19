tic
clear all
close all
clc

a = 10;
b = 20;
c = 20;
dt = 0.5;
BCL = 800/dt;       
endpoint=700*a+550*b+462*c;
tstep1 = 700/dt;
tstep2 = 550/dt;
tstep3 = 462/dt;
cycle = 0;
ts = (0:dt:endpoint);
t = ts(2:1:length(ts));
ind = [1:1:400];
network_matrix = load('network_matrix_test_400vessel.txt');
fm = initial_variable(network_matrix(ind,:));
indx = load('indx.mat');
load('cond_leng.mat');
cond_leng = cond_leng/1000; 
%% lumped parameter
VLV = 115;
Vart = 1869;
Vad = 243;
Vvc = 3104;
Vla = 38;
Rao = 300;
Rvc = 900;
Rper = 18000;
Rad = 105000;     
Rmv = 1900;
Cao = 0.006;          
Cad = 0.02;           
Cvc = 0.8;
Vart0 = 1800;
Vad0 = 20;
Vvc0 = 1400;          
Ees_LV = 4.02/0.0075;
A_LV = 1/0.0075;           
B_LV = 0.027;
V0_LV = 10;
Tmax_LV = 280;
tau_LV = 25;
trans_LV = 325;
cycle_LV = 1;
Ees_la = 0.6/0.0075;
A_la = 0.44/0.0075;
B_la = 0.05;
V0_la = 10;
Tmax_la = 150;
tau_la = 25;
cycle_la = 1;
Fmeta=1;
alpha = 10;
gamma = 0.6;
deltaP_ = zeros(1, length(ind));
tau_    = zeros(1, length(ind));
C_      = zeros(1, length(ind));
G_      = zeros(1, length(ind));
Rat_    = zeros(1, length(ind));
for i=1:1:length(ts)-1
    if i == 1
        Rad = 105000;
        Cao = 0.006;
        Cad = 0.02;
        Vvc0 = 1700;
        BCL = 700;
        HR = 60/BCL*1000;  
        Ees_LV = 4.2235/0.0075;     
        VLV = 109;
        Vart = 1872;
        Vvc = 3056;
        Vla = 81;
        Vad = 250;
    elseif i == tstep1*a+1
        BCL = 550;
        HR = 60/BCL*1000;
        Rad = 105000*0.7;
        Cao = 0.006/6;
        Cad = 0.02/4;
        Vvc0 = 900;
        Ees_LV = 4.6386/0.0075;
        VLV = 42;
        Vart = 1817;
        Vvc = 3315;
        Vla = 88;
        Vad = 105;       
    elseif i == tstep1*a+tstep2*b+1
        BCL = 462;
        HR = 60/BCL*1000;
        Rad = 105000*0.5;
        Cao = 0.006/6;
        Cad = 0.02/4;
        Vvc0 = 780-1000;
        Ees_LV = 4.9597/0.0075;
        VLV = 41;
        Vart = 1815;
        Vvc = 3302;
        Vla = 95;
        Vad = 115;
        A_LV = 1/0.0075*1.2;
        Tmax_LV = 280-100;  
    end
    if(i>1 && mod(t(i-1),BCL) == 0)
        cycle      = cycle + 1;
        PVA(cycle) = polyarea(PVQ(i-BCL/dt:1:i-1,1),PVQ(i-BCL/dt:1:i-1,2));
        Qm(cycle)  = mean(PVQ(i-BCL/dt:1:i-1,5));
        Fmeta = PVA(cycle)*HR /1000*0.00076471-0.4011; 
        if Fmeta < 0
            Fmeta = 0;
        elseif Fmeta > 1
            Fmeta = 1;
        end
        Ees_LV1 = (Qm(cycle)*0.2/0.4105+1.4587/0.4105)/0.0075;
        B = log(2)/85;
        A = 1/exp(75*B);
        Ees_LV2 = (A*exp(HR/1*B)-1)/0.0075*1;
        Ees_LV = Ees_LV1 + Ees_LV2;
        Fm(cycle) = Fmeta;
        Ees(cycle) = Ees_LV;
    end            
     
    Part = 1/Cao*(Vart - Vart0);
    Pad = 1/Cad*(Vad - Vad0);
    Pvc = 1/Cvc*(Vvc - Vvc0);
    %% For calculating PLV
    if (mod(ts(i),BCL) == 0)
        cycle_LV = ts(i)/BCL;
    end
    t_LV = ts(i)-cycle_LV*BCL;
    PLV = e(t_LV,Tmax_LV,tau_LV, trans_LV).*Ees_LV.*(VLV - V0_LV) + (1 - e(t_LV,Tmax_LV,tau_LV,trans_LV)).*A_LV.*(exp(B_LV*(VLV - V0_LV)) - 1);    %Pa
    %% For calculating PLA
    trans_la = 3/2*Tmax_la;
    if (mod((ts(i)),BCL) == 0)
        cycle_la = ts(i)/BCL + 1;
    end
    t_la = ts(i)-(cycle_la-1)*BCL;
    Pla = e(t_la,Tmax_la,tau_la,trans_la).*Ees_la.*(Vla - V0_la) + (1 - e(t_la,Tmax_la,tau_la,trans_la)).*A_la.*(exp(B_la*(Vla - V0_la)) - 1);
    Elv(i) = e(t_LV,Tmax_LV,tau_LV, trans_LV); 
    IMP(i) = PLV*0.6+3*Ees_LV*Elv(i)+40*(1-VLV/114)/0.0075;    
    meandeltaP = mean(deltaP_,1);
    meantau    = mean(tau_,1);
    [Qi, deltaP, tau, C, G, Rat] = flow_analysis(Fmeta, i, ts*0.001, 0.80*Part, Pvc, IMP, meandeltaP, meantau, G_, C_, Rat_, network_matrix, dt*0.001, fm, ind, indx, cond_leng);
    deltaP_(i+1,:) = deltaP;
    tau_(i+1,:)    = tau;
    C_(i+1,:)      = C;
    G_(i+1,:)      = G;
    Rat_(i,:)      = Rat;
    if(PLV <= Part)  
        Qao = 0;
    else
        Qao = 1/Rao*(PLV - Part);
    end
    if(PLV >= Pla)   
        Qla = 0;
    else
        Qla = 1/Rmv*(Pla - PLV);
    end
    Qper = 1/Rper*(Part - Pad);
    Qad = 1/Rad*(Pad - Pvc); 
    Qvc = 1/Rvc*(Pvc - Pla);   
    VLV = VLV + dt*(Qla - Qao);
    Vart = Vart + dt*(Qao - Qper - Qi(1)/60);
    Vad = Vad + dt*(Qper - Qad); 
    Vvc = Vvc + dt*(Qad - Qvc) + Qi(1)/60);
    Vla = Vla + dt*(Qvc - Qla);
    PVQ(i,1) = VLV;
    PVQ(i,2) = PLV*0.0075;
    PVQ(i,3) = Part*0.0075;
    PVQ(i,4) = IMP(i)*0.0075;
    PVQ(i,5) = Qi(:,1)/34.2105;
end

plotind = find(t >= (endpoint-BCL+1));
results.EF = (max(PVQ(plotind,1))-min(PVQ(plotind,1)))/max(PVQ(plotind,1))*100;
results.LVEDV = max(PVQ(plotind,1));
LVESV = min(PVQ(plotind,1));
Pind = find(PVQ == LVESV);
results.LVESV = min(PVQ(plotind,1));
results.LVEDP = PVQ(find(max(PVQ(plotind,1))),2);
results.LVESP = PVQ(Pind(1),2);
results.DBP = min(PVQ(plotind,3));
results.SBP = max(PVQ(plotind,3));
results.PVA = [PVA(a), PVA(a+b-1), PVA(end)];
results.Fmeta = [Fm(a),Fm(a+b-1),Fm(end)];
results.Ees = [Ees(a)*0.0075,Ees(a+b-1)*0.0075,Ees(end)*0.0075];
results.Qmean = [Qm(a), Qm(a+b-1), Qm(end)];

save('stenosisFFR8_reg.mat','PVQ') 

toc
