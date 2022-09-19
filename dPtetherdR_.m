function dPtetherdR = dPtetherdR_(R, k, long_dist, A, Rm_avg, R0, Cstr)
     dPtetherdR = (R < R0).*(long_dist.cp(k).*((sec(((R + (A.*Rm_avg)' - long_dist.bp(k))*pi)./(long_dist.ap(k)-long_dist.bp(k)) - pi/2).^2).*pi)./(long_dist.ap(k)-long_dist.bp(k)) + ...
                             2*Cstr*(R0 - R)./R + Cstr*(R0 - R).^2./R.^2) + ...
                  (R >= R0).*(long_dist.cp(k).*((sec(((R + (A.*Rm_avg)' - long_dist.bp(k))*pi)./(long_dist.ap(k)-long_dist.bp(k)) - pi/2).^2).*pi)./(long_dist.ap(k)-long_dist.bp(k))); 
         
end