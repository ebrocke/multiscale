function [inf tau] = get_k_h(v)
    tadj = 1;
	a = 1.0/((exp((v+46.05)/5)+exp(-(v+238.4)/37.45))) / tadj;
    if v<-63.0
        tau = a;
    else
        tau = 19.0/tadj;
    end
    inf = 1.0/(1+exp((v+78)/6));
    tau = tau*1e-3;
   % alpha = 1e3*inf/tau;
   % beta = 1e3*(1-inf)/tau;
end