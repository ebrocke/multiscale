function [inf tau] = get_k_p(v)
global CONST
    tau_max = 1e-3*CONST.TAUMAX_IM_SOMA; %s
    inf =  1.0/(1+exp(-(v+35.0)/10.0));
    tau = tau_max/(3.3*exp((v+35.0)/20.0)+exp(-(v+35.0)/20.0));
   
    %alpha = inf/tau;
    %beta = (1-inf)/tau;
end