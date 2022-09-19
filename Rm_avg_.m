function Rm_avg = Rm_avg_(k, long_dist, rhom, deltaP_avg)
    Rm_avg = (long_dist.rhoa(k)/pi*rhom).*( 0.5*pi - atan(((deltaP_avg' - ...
                long_dist.pha(k) )./long_dist.ca(k) ).^(2.0*long_dist.ma(k))));
end