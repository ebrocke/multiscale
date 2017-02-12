% 1. 2Ca+Raf<=>Active_Raf
% 2. Active_Raf+MAPK<=>AM=>P_MAPK+Active_Raf
function daraf = get_daraf(sps, rrs)
K = [rrs(1,:) rrs(2,:)];
R = [sps(1)*sps(1)*sps(2), -sps(3), 0.0, ... %ca*ca*raf, -araf
     -sps(3)*sps(4), sps(5), sps(5)]; %araf*mapk, mapk_araf, mapk_araf
daraf = R*K';
end