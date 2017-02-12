% 2. Active_Raf+MAPK<=>AM=>P_MAPK+Active_Raf
function dmapkb = get_dmapk_araf(sps, rrs)
K = rrs(2,:);
R = [ sps(3)*sps(4), -sps(5), -sps(5)]; %araf*mapk -mapk_araf -mapk_araf
dmapkb = R*K';
end