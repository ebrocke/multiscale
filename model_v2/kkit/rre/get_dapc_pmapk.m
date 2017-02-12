% 8. P_MAPK+APC<=>PA=>AA+P_MAPK
function dapcb = get_dapc_pmapk(sps,rrs)
% global APC_BASE
APC_BASE = 1e-6;
K = rrs(8,:);
R = [sps(6)*APC_BASE, -sps(15), -sps(15)]; %pmapk*apc -apc_mpapk -apc_pmapk
dapcb = R*K';
end