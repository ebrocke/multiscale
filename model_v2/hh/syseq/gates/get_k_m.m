function [inf tau] = get_k_m(v)
    tadj = 1;
	tau = 1e-3*1.0/((exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7))+0.37) / tadj;
	inf = 1.0 / (1+exp(-(v+60)/8.5));

    
    %alpha = 1e3*inf/tau;
    %beta = 1e3*(1-inf)/tau;
end