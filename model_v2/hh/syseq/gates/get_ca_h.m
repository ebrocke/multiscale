function [inf tau] = get_ca_h(v)
    %alpha = 1e3*1.6/(1+exp(-(5+v)/13.89));
    %beta  = 1e3*0.02*vtrap(v+18.9,5);
    
    alpha = 0.457*exp((-13.0-v)/50.0);
	beta = 6.5/(exp((-v-15.0)/28.0) + 1);
    
    tau = 1/ (alpha+beta);
    inf = alpha*tau;
end