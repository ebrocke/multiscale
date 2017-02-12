% 2. Active_Raf+MAPK<=>AM=>P_MAPK+Active_Raf
% 3. P_MAPK+Phosphotase<=>PP=>MAPK+Phosphotase
% 9. MAPK+Act_PKC<=>MA=>P_MAPK+Act_PKC
function dmapk = get_dmapk(sps, rrs)
K = [rrs(2,:) rrs(3,:) rrs(9,:)];
R = [-sps(3)*sps(4), sps(5), 0.0, ... %-araf*mapk, mapk_araf
    0.0, 0.0, sps(8), ...        %pmapk_phsph
    -sps(4)*sps(13), sps(16), 0.0];      %-mapk*apkc mapk_apkc
dmapk = R*K';
end