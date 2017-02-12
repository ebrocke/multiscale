% 10. Ca+pmca<=>pmcaCa=>pmca
function dpmcaCa = get_dpmcaCa(sps,rrs)
K = rrs(10,:);
R = [sps(1)*sps(17), -sps(18), -sps(18)];%ca*pmca, -pmcaCa 
dpmcaCa = R*K';
end