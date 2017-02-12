% 4. P_MAPK+K_A<=>PK=>P_K_A+P_MAPK
% 5. P_K_A=>K_A
function dka = get_dka(sps, rrs)
K = [rrs(4,:) rrs(5,:)];
R = [-sps(9)*sps(6), sps(10), 0.0, ...%-ka*pmapk pmapk_ka
    sps(11), 0.0, 0.0]; %pka
dka = R*K';
end