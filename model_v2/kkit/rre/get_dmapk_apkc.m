% 9. MAPK+Act_PKC<=>MA=>P_MAPK+Act_PKC
function dmapkb = get_dmapk_apkc(sps,rrs)
K=rrs(9,:);
R=[sps(4)*sps(13), -sps(16), -sps(16)];%mapk*apkc -mapk_apkc -mapk_apkc
dmapkb = R*K';
end