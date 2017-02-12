% 1. 2Ca+Raf<=>Active_Raf
function draf = get_draf(sps, rrs)
K = rrs(1,:);
R = [-sps(1)*sps(1)*sps(2), sps(3), 0.0]; %-ca*ca*raf araf
draf = R*K';
end