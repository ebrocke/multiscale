% 3. P_MAPK+Phosphotase<=>PP=>MAPK+Phosphotase
function dpmapkb = get_dpmapk_phsph(sps, rrs)
K = rrs(3,:);
R = [sps(6)*sps(7), -sps(8), -sps(8)]; %pmapk*phsph -pmapk_phsph -pmapk_phsph
dpmapkb = R*K';
end