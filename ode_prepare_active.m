function [Pcdot,PT,dPTdt]=ode_prepare_active(t, Pc,network_matrix,n,...
    full_Pin,full_Pout,full_PT,full_dPTdt,fm,G,C)

[PT]=full_PT';%*network_matrix(:,17);
[dPTdt]=full_dPTdt';%*network_matrix(:,17);
[Pin]=full_Pin;
Pout = full_Pout*network_matrix(:,17);
% [Pcdot,Rat,dRadP,deltaP]=...
%     matrices_AandB_active(PT,dPTdt,Pin,Pout,Pc,deltaP_avg,Rss,R0,network_matrix,n,fm,kdyn);
%Pcdot = matrices_AandB_active(dPTdt,Pin,Pout,Pc,n,fm,G,C);
% A matrix and B vector
flow_matrix=zeros(n,1);
flow_vector=dPTdt;
Gm=zeros(n,1);
Gs = zeros(n,1);
nnt = length(fm.Gnt_loc);
nt = length(fm.Gt_loc);
Gs2 = zeros(nnt,1);
Gsi2_t = zeros(nt,1);
Gd1 =zeros(n,1);
Gd2 = zeros(n,1);
Gd3 = zeros(nnt,1);
Pcm = zeros(n,1);
Pcs = zeros(n,1);
Pcs2 = zeros(nnt,1);
Pcs2_t = zeros(nt,1);
Pcd1 = zeros(n,1);
Pcd2 = zeros(n,1);
Pcd3 = zeros(nnt,1);
Gm = G(fm.mo);
Gs = G(fm.si);
Gs2(fm.si2) = G(fm.si2index);
Gd1 = G(fm.da1);
Gd2 = G(fm.da2);
Gd3(fm.da3) = G(fm.da3index);
Gnt = G(fm.Gnt_loc);
Cnt = C(fm.Gnt_loc);
Gm_t  = G(fm.mo_t);
Gs_t = G(fm.si_t);
Gsi2_t(fm.si2_t) = G(fm.si2tindex);
G_t = G(fm.Gt_loc);
C_t = C(fm.Gt_loc);
Pcm = Pc(fm.mo);
Pcs = Pc(fm.si);
Pcs2(fm.si2) = Pc(fm.si2index);
Pcd1 = Pc(fm.da1);
Pcd2 = Pc(fm.da2);
Pcd3(fm.da3) = Pc(fm.da3index);
Pc_nt = Pc(fm.Gnt_loc);
Pcm_t = Pc(fm.mo_t);
Pcs_t = Pc(fm.si_t);
Pcs2_t(fm.si2_t)=Pc(fm.si2tindex);
Pc_t = Pc(fm.Gt_loc);
one_by_den1 = ones(nnt,1)./(Cnt.*(Gnt+Gm+Gs+Gs2));
one_by_den2 = ones(nnt,1)./(Cnt.*(Gnt + Gd1+ Gd2 + Gd3));
flow_matrix(fm.Gnt_loc) = (2*Gm.*Gnt.*Pcm.*one_by_den1) + ...
    (2*Gs.*Gnt.*Pcs.*one_by_den1)+...
    (2*Gs2.*Gnt.*Pcs2.*one_by_den1)+...
    (-4*Gnt.*Pc_nt./Cnt + 2*Gnt.*Gnt.*Pc_nt.*one_by_den1+ ...
    2*Gnt.*Gnt.*Pc_nt.*one_by_den2)+ ...
    (2*Gd1.*Gnt.*Pcd1.*one_by_den2)+ ...
    (2*Gd2.*Gnt.*Pcd2.*one_by_den2)+ ...
    (2*Gd3.*Gnt.*Pcd3.*one_by_den2);
one_by_den3 = ones(nt,1)./(C_t.*(G_t+Gm_t+Gs_t+Gsi2_t));
flow_matrix(fm.Gt_loc) = (2*Gm_t.*G_t.*Pcm_t.*one_by_den3) + ...
    (2*Gs_t.*G_t.*Pcs_t.*one_by_den3)+...
    (2*Gsi2_t.*G_t.*Pcs2_t.*one_by_den3)+...
    (-4*G_t.*Pc_t./C_t + 2*G_t.*G_t.*Pc_t.*one_by_den3);
flow_matrix(1) =      (-4*G(1).*Pc(1)./C(1) + ...
    2*G(1).*G(1).*Pc(1)./(C(1).*(G(1) + Gd1(1)+ Gd2(1))))+ ...
    (2*Gd1(1).*G(1).*Pcd1(1)./(C(1).*(G(1) + Gd1(1)+Gd2(1))))+ ...
    (2*Gd2(1).*G(1).*Pcd2(1)./(C(1).*(G(1)+Gd1(1)+Gd2(1))));
flow_vector(fm.Gt_loc)=dPTdt(fm.Gt_loc) + (2*G_t./C_t).*Pout(fm.Gt_loc);
flow_vector(1)=dPTdt(1)+2*G(1)*Pin/C(1);
Pcdot = flow_matrix+ flow_vector;


Pcdot=Pcdot(1:1:n,1);
PT=PT(1:1:n,1);
dPTdt=dPTdt(1:1:n,1);