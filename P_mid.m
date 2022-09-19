function vessels_pressure = ...
    P_mid(Pa,Pv,PT,dPTdt,network_matrix,p_lumped_init,n,fm,t1,t2,G,C)
time_span = [t1, t2];
tolerance=1e-6;
% vessel_end_pressure = p_lumped_init;
options=odeset('RelTol',tolerance*1e-3,'AbsTol',tolerance*ones(n,1),'BDF','on','JPattern',jpattern(network_matrix,n));
clear time vessels_pressure
initial_pressure=p_lumped_init;
[time,vessels_pressure]= ...
    ode15s(@ode_prepare_active,time_span,initial_pressure,options,network_matrix,...
    n,Pa,Pv,PT,dPTdt,fm,G,C);
   
function S = jpattern(nm,n)
Jac = zeros(n,n);
Jac(1,1) = 1;
Jac(1,nm(1,12))=1;
Jac(1,nm(1,13))=1;
 for i=2:n
    if (nm(i,8)~=0 && nm(i,9)==0 && nm(i,12)==0 && nm(i,13)==0) %for the external arteries
      Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
    elseif (nm(i,9)~=0 && nm(i,12)==0 && nm(i,13)==0) %for the external arteries
      Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
      Jac(i,nm(i,9))=1;
    elseif(nm(i,8)==0 && nm(i,12)==0 && nm(i,13)==0)
        Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
    elseif(nm(i,8)==0 && nm(i,12)==0)
        Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
   elseif(nm(i,8)==0 && nm(i,13)==0)
         Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,12))=1;
    elseif(nm(i,8)==0 && nm(i,12)~=0)
        Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,12))=1;
      Jac(i,nm(i,13))=1;
    elseif(nm(i,8)~=0 && nm(i,13)==0)
      Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
      Jac(i,nm(i,12))=1;
    elseif(nm(i,14)~=0)
      Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
      Jac(i,nm(i,12))=1;
      Jac(i,nm(i,13))=1;
      Jac(i,nm(i,14))=1;
    elseif(nm(i,9)~=0)
        Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
      Jac(i,nm(i,9))=1;
      Jac(i,nm(i,12))=1;
      Jac(i,nm(i,13))=1;
    else
      Jac(i,i)=1;
      Jac(i,nm(i,7))=1;
      Jac(i,nm(i,8))=1;
      Jac(i,nm(i,12))=1;
      Jac(i,nm(i,13))=1;
    end
 end

 S = sparse(Jac);