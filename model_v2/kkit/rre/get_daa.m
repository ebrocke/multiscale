% 6. PKC+2AA<=>Act_PKC
% 7. AA<=>APC
% 8. P_MAPK+APC<=>PA=>AA+P_MAPK
function daa = get_daa(sps,rrs)
% global APC_BASE
APC_BASE = 1e-6;
K=[rrs(6,:) rrs(7,:) rrs(8,:)];
R=[-2*sps(14)*sps(14)*sps(12), 2*sps(13), 0.0,... %-2*aa*aa*pkc 2*apkc
   -sps(14), APC_BASE, 0.0,... %-aa apc
   0.0,0.0, sps(15)]; % apc_pmapk
daa = R*K'; 
end