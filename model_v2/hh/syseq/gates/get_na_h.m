function [inf tau] = get_na_h(v)
    global CONST
    v2 = v - CONST.VTRAUB_SOMA;
    alpha = 128*exp((17.0-v2)/18.0);
    beta  = 4000/(exp((40.0-v2)/5.0)+1);
    
    tau = 1/ (alpha+beta);
    inf = alpha*tau;

end