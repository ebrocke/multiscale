function [inf tau] = get_na_m(v)
    global CONST
    v2 = v - CONST.VTRAUB_SOMA;
    alpha = 320*vtrap((13.0-v2),4.0);
    beta  = 280*vtrap((v2-40.0),5.0);
   
    
    tau = 1/ (alpha+beta);
    inf = alpha*tau;
end