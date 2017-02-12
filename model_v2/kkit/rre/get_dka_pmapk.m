% 4. P_MAPK+K_A<=>PK=>P_K_A+P_MAPK
function dkab = get_dka_pmapk(sps,rrs)
K = rrs(4,:);
R = [sps(6)*sps(9), -sps(10), -sps(10)]; %pmapk*ka -ka_pka -ka_pka
dkab = R*K';
end