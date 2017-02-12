function [inf tau] = get_ca_m(v)
% if v > -70.0
%     alpha = 0.005/exp(0.05*(v+70));
% else
%     alpha = 0.005;
% end
% beta  = 1e3*(0.005-alpha);
% alpha = alpha*1e3;

alpha = 55*vtrap(-(27.0+v),3.8);
beta = 940*exp((-75.0-v)/17.0);

tau = 1/ (alpha+beta);
inf = alpha*tau;
end