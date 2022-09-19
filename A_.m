function A = A_(k, Fmeta, long_dist, tau_avg)
    %ftau = (long_dist.fmax(k).*tau_avg)./(long_dist.ktau(k)*35 + tau_avg);
    ftau = (long_dist.fmax(k).*tau_avg(k))./(long_dist.ktau(k)*35 + tau_avg(k));
%     ftau(1) = -0.5;
%     ftau(2) = 3;
    ftau(find(ftau<0)) = 0;
    ftau(find(ftau>1)) = 1;
    A = (1 - ftau).*(1 - Fmeta(k));
end