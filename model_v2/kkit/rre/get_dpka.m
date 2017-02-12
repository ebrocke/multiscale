% 4. P_MAPK+K_A<=>PK=>P_K_A+P_MAPK
% 5. P_K_A=>K_A
function dpka = get_dpka(sps, rrs)
K = [rrs(4,:) rrs(5,:)];
R = [0.0, 0.0, sps(10),... %ka_pka
     -sps(11),0.0,0.0]; %-ka_pka
dpka = R*K';
end