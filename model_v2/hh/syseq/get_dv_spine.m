function dv_spine = get_dv_spine(in, ra, frac, cai)
global CONST NGATES NGATES_SOMA
v = in(NGATES+2);
v_dend = in(end);

%high threshold Ca L current
ical = get_ica_spine(v, in(NGATES_SOMA+1), in(NGATES_SOMA+2), cai);

%transient potassium current
if (frac >= 0.0 )
    gka = CONST.GKBAR_IKA_SPINE*frac;
else
    gka = CONST.GKBAR_IKA_SPINE*in(NGATES_SOMA+3)^4*in(NGATES_SOMA+4);
end
ika = (v-CONST.EK_SPINE)*gka;

%leak current
il = (v-CONST.EM_SPINE)*CONST.G_PAS_SPINE;

dv_spine = -(ical+ika+il+(v-v_dend)/ra(5))/CONST.CM_SPINE;
end