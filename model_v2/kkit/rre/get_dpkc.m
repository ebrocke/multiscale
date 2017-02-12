% 6. PKC+2AA<=>Act_PKC
function dpkc = get_dpkc(sps,rrs)
K = rrs(6,:);
R = [-sps(12)*sps(14)*sps(14), sps(13), 0.0]; %-pkc*aa*aa apkc
dpkc = R*K';
end