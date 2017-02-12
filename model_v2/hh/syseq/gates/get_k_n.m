function [inf tau] = get_k_n(v)
    global CONST
    v2 = v - CONST.VTRAUB_SOMA;
    alpha = 32*vtrap((15.0-v2),5.0);
    beta  = 500*exp((10.0-v2 )/40.0);
    
    tau = 1/ (alpha+beta);
    inf = alpha*tau;

end