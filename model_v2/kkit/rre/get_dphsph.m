% 3. P_MAPK+Phosphotase<=>PP=>MAPK+Phosphotase
function dphsph = get_dphsph(sps, rrs)
K = rrs(3,:);
R = [-sps(6)*sps(7), sps(8), sps(8)]; %-pmapk*phsph, pmapk_phsph, pmapk_phsph
dphsph = R*K';
end