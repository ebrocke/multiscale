% 1.  2Ca+Raf<=>Active_Raf
% 10. pmca+Ca <=> pmcaCa => pmca
function dca=get_dca_pump(sps, rrs)
K = [rrs(1,:) rrs(10,:)];

R = [-2*sps(1)*sps(1)*sps(2), 2*sps(3), 0.0, ... %-2*ca*ca*raf, 2*araf
    -sps(1)*sps(17), sps(18), 0.0];     %-ca*pmca, pmcaCa
dca = R*K';
end

