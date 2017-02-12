function ical = get_ica_spine(v, m, h, cai)
global CONST 
global CELL_PARAMS_CAREV0 CELL_PARAMS_CAREV1
carev = CELL_PARAMS_CAREV1 - CELL_PARAMS_CAREV0*log(cai);
%high threshold Ca L current
ical = (v-carev)*CONST.GKBAR_ICAL_SPINE*m^2*h;
end