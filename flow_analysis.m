function [Qi, deltaP, tau, C, G, Rat] = flow_analysis(Fmeta_t, i, t, Part, Pven, IMP, meandeltaP, meantau, G_, C_, Rat_, network_matrix, dt, fm, ind, indx, cond_leng)
index1 = find(network_matrix(:,4)==1 & network_matrix(:,12)==0);
n=length(ind);%length(network_matrix(:,1));
Fmeta = vasodilator_signal(Fmeta_t, index1, indx, cond_leng, n);
Cstr = 2.0;                                                                % Tethering constant
rhom = 1.2;%3.6;                                                                % Myogenic parametes
mu = 2.7e-9;
% Compute vessel parameters
[long_dist]=VesselParams1(network_matrix(ind,:),n);
R0 = long_dist.bp(ind) + (long_dist.ap(ind)-long_dist.bp(ind))/pi.*(pi/2 + ...
    atan((-long_dist.php(ind))./long_dist.cp(ind)));                       % passive radius
A(ind) = A_(ind, Fmeta, long_dist, zeros(400,1));                                       % tau_ave = 0
%mmHg2MPa = 0.000133322;                                                    % Convert unit
%t = linspace(0, (length(Part)-1)*dt, length(Part));                        % Create timesteps
% Scale Part, Pven and IMP to MPA;
Part = Part*1e-6;
Pven = Pven*1e-6;
PT(i,:) = IMP(i)*1e-6*ones(1,n);
if i == 1
    Pm = zeros(1,length(ind));
    deltaP(ind) = Pm(ind) - PT(ind);
    Rm_avg(ind) = Rm_avg_(ind, long_dist, rhom, deltaP(ind));
    Rat(ind) = (atan((deltaP(ind)'-long_dist.php(ind))./long_dist.cp(ind))+pi/2).*(long_dist.ap(ind)-long_dist.bp(ind))/pi+long_dist.bp(ind)-A(ind)'.*Rm_avg(ind)';
    dRatdP = dPtetherdR_(Rat(ind)', ind, long_dist, A(ind), Rm_avg(ind), R0(ind), Cstr);
    % Initialize resistance and capacitance
    R = (8.0*mu*network_matrix(ind,2)*1e-3)./(pi*Rat(ind)'.^4);
    G = ones(n,1)./R*1000;
    C = ((2.0*pi*network_matrix(ind,2)*1e-3).*Rat(ind)').*dRatdP(ind);
    dPTdt = PT/dt;
    t1 = t(i);
    t2 = t(i+1);
    Pmm = P_mid(Part(1),Pven(1),PT(ind),dPTdt(ind),network_matrix(ind,:),deltaP,n,fm,t1,t2,G(ind),C(ind,1));
    Pm(1,:) = Pmm(size(Pmm,1),:);
    
    flow = check_delta_flow_sim(network_matrix,Part(1),Pven(1),dPTdt,G',C',Pmm(size(Pmm,1),:));
    Qi(1,:) = flow.flow_in; Qo(1,:) = flow.flow_out;
    Pnode = flow.Pnode;
    tau(1,1) = (Part(1) - Pnode(1)).*Rat(1)./(2*network_matrix(1,2)*1e-3);% Calculate shear stress in each vessel % Shear stress in the source vessel
    tau(1,2:1:n)= abs((Pnode(network_matrix(2:1:n,5))-Pnode(network_matrix(2:1:n,6))))...
        .*Rat(2:1:n)./(2*network_matrix(2:1:n,2)'*1e-3);
else
    dPTdt = (PT(i) - PT(i-1))/dt*ones(1,n);
    Pmm = P_mid(Part,Pven,PT(i,:),dPTdt,network_matrix,meandeltaP(1,:),n,fm,t(i),t(i+1),G_(i,:)',C_(i,:)');
    Pm(1,:) = Pmm(size(Pmm,1),:);
    flow = check_delta_flow_sim(network_matrix,Part,Pven,dPTdt,G_(i,:),C_(i,:),Pm(1,:));
    Qi(1,:) = flow.flow_in; Qo(1,:) = flow.flow_out;
    %Pn(i+1,:) = flow.Pnode;
    Pnode = flow.Pnode;
    tau(1,1) = (Part - Pnode(2)).*Rat_(i-1,1)./(2*network_matrix(1,2)*1e-3);% Calculate shear stress in each vessel % Shear stress in the source vessel
    tau(1,2:1:n)= abs((Pnode(network_matrix(2:1:n,5))-Pnode(network_matrix(2:1:n,6))))...
        .*Rat_(i-1,2:1:n)./(2*network_matrix(2:1:n,2)'*1e-3);
    %tau_avg(i,:) = mean(tau,1);                                              % mean over current cycle?
    %Fmeta = vasodilator_signal(Fmeta_t, index1, indx, cond_leng, n);
    A(1,:) = A_(ind, Fmeta, long_dist, meantau');
    deltaP(1,:) = Pm(1,:) - PT(i,:);
    %meandeltaP(1,:) = mean(deltaP,1);
    %meanP = mean(Pi(1:i));
    Rm_avg(1,:) = Rm_avg_(ind, long_dist, rhom, meandeltaP(1,:));
    % Get vessel radius using Rmin as initial guess
        % Find minimum R
    Rmin = - A(1,:).*Rm_avg(1,:) + long_dist.bp(ind)';                       % ?
    Rmin(find(Rmin<0)) = 1e-6;
%     Rmax = pi*(long_dist.ap(ind)'-long_dist.bp(ind)')/pi - A(1,:).*Rm_avg(1,:) + long_dist.bp(ind)';    
%     Rshift = Rmax - Rmin; 
    
    Rp(1,:) = (atan((deltaP(1,:)'-long_dist.php(ind))./long_dist.cp(ind))+pi/2).*(long_dist.ap(ind)-long_dist.bp(ind))/pi+long_dist.bp(ind);
    Rat(1,:) = (atan((deltaP(1,:)'-long_dist.php(ind))./long_dist.cp(ind))+pi/2).*(long_dist.ap(ind)-long_dist.bp(ind))/pi+long_dist.bp(ind)-A(1,:)'.*Rm_avg(1,:)';
    %Rat(find(Rat<0)) = (atan((deltaP(1,find(Rat<0))'-long_dist.php(find(Rat<0)))./long_dist.cp(find(Rat<0)))+pi/2).*(long_dist.ap(find(Rat<0))-long_dist.bp(find(Rat<0)))/pi+long_dist.bp(find(Rat<0))-0.8*A(1,find(Rat<0))'.*Rm_avg(1,find(Rat<0))';
    Rat(find(Rat<0)) = Rmin(find(Rat<0));
    %Rat(1,:) = fsolve(@(y) Ptether_(y, ind, long_dist, A(1,:), Rm_avg(1,:), R0, Cstr, deltaP(1,:)), Rmin, optimoptions('fsolve','Display','off','TolFun',1e-9, 'TolX', 1e-9));
%     while(length(find(Rat(1,:) < Rmin)) >0)
%         Rat(1,find(Rat(1,:) < Rmin)) = Rat(1,find(Rat(1,:) < Rmin)) + Rshift(1,find(Rat(1,:) < Rmin));
%     end 
%     while(length(find(Rat(i+1,:) > Rmax)) >0)
%         Rat(i+1,find(Rat(i+1,:) > Rmax)) = Rat(i+1,find(Rat(i+1,:) > Rmax)) + Rshift(1,find(Rat(i+1,:) > Rmax));
%     end  
    % Compute ddeltaPdR
    dRatdP = dPtetherdR_(Rat(1,:)', ind, long_dist, A(1,:), Rm_avg(1,:), R0, Cstr);   % dPdR?
    R = (8.0*mu*network_matrix(ind,2)*1e-3)./(pi*Rat(1,:)'.^4);
    G = ones(400,1)./R*1000;
    C = ((2.0*pi*network_matrix(ind,2)*1e-3).*Rat(1,:)').*dRatdP;
    %GG(i+1,:) = G';
    % For each cycle
%     if(mod(t(i+1),BCL) == 0)
%         endind = i;
%         Qitotal = trapz(t(startind:endind), Qi(startind:endind,3))/60;
%         Qototal = trapz(t(startind:endind), Qi(startind:endind,2))/60;
%         fprintf('%s:  %-d, %s: %-d, %s: %-d\n', "Cycle", cycle, "Start index", startind, "End index", endind);
%         fprintf('%s:  %10.5f, %s: %10.5f, %s: %10.5f\n', "Qin", Qitotal, "Qout", Qototal, "MeanP", meanP/mmHg2MPa);
%         startind = i+1;
%     end
end