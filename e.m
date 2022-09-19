
function out = e(t, Tmax, tau, trans)
      if t <= trans
		out = 0.5.*(sin((pi/Tmax).*t - pi/2) + 1)+0;
      else
		%out = 0.5.*exp((-t + (1.5*Tmax))/tau);
        out = 0.5.*(sin((pi/Tmax).*trans - pi/2) + 1)*exp((-t + (trans))/tau);
      end 	
end