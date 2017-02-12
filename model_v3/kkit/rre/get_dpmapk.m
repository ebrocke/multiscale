% 2. Active_Raf+MAPK<=>AM=>P_MAPK+Active_Raf
% 3. P_MAPK+Phosphotase<=>PP=>MAPK+Phosphotase
% 4. P_MAPK+K_A<=>PK=>P_K_A+P_MAPK
% 8. P_MAPK+APC<=>PA=>AA+P_MAPK
% 9. MAPK+Act_PKC<=>MA=>P_MAPK+Act_PKC

function dpmapk = get_dpmapk(sps, rrs)
% global APC_BASE
APC_BASE = 1e-6;
K = [rrs(2,:), rrs(3,:), rrs(4,:), rrs(8,:), rrs(9,:) ];
R = [0.0, 0.0, sps(5), ...  %mapk_araf 
    -sps(6)*sps(7), sps(8), 0.0, ... %-pmapk*phsph pmapk_phsph 
    -sps(6)*sps(9), sps(10), sps(10), ...%-pmapk*ka ka_pmapk ka_pmapk 
    -sps(6)*APC_BASE, sps(15), sps(15),...%-pmapk*apc apc_pmapk apc_pmapk 
    0.0, 0.0, sps(16)]; %mapk_apkc
dpmapk = R*K';
end