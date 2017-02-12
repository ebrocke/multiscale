% 6. PKC+2AA<=>Act_PKC
% 9. MAPK+Act_PKC<=>MA=>P_MAPK+Act_PKC
function dapkc = get_dapkc(sps,rrs)
K=[rrs(6,:) rrs(9,:)];
R=[ sps(12)*sps(14)*sps(14), -sps(13), 0.0, ...%pkc*aa*aa -apkc
    -sps(4)*sps(13), sps(16), sps(16)]; %-mapk*apkc mapk_apkc mapk_apkc
dapkc = R*K';
end