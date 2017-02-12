% 1.  2Ca+Raf<=>Active_Raf
function dca=get_dca(sps, rrs)
K = [rrs(1,:) rrs(10,:)];
R = [-2*sps(1)*sps(1)*sps(2), 2*sps(3), 0.0];%-2*ca*ca*raf, 2*araf
      
dca = R*K';
end

