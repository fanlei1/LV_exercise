function Ptether = Ptether_(R, k, long_dist, A, Rm_avg, R0, Cstr, P)
    Ptether = (R < R0).*(long_dist.php(k) + long_dist.cp(k).*tan(((R + (A.*Rm_avg)' - long_dist.bp(k))*pi)./(long_dist.ap(k)-long_dist.bp(k)) - pi/2) - Cstr*(R0 - R).^2./R) + ...
              (R >= R0).*(long_dist.php(k) + long_dist.cp(k).*tan(((R + (A.*Rm_avg)' - long_dist.bp(k))*pi)./(long_dist.ap(k)-long_dist.bp(k)) - pi/2)) - P';  
end