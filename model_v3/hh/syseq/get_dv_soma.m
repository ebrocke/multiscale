function dv_soma = get_dv_soma(in, ra, inj)
global CONST NGATES
v = in(NGATES+1);
%we put dendrite voltages at the end of the state vector
v_dend = in(end-CONST.NSEG_DEND+1);

%voltage-dependent Na current
ina = (v-CONST.ENA_SOMA)*CONST.GKBAR_INA_SOMA*in(1)^3*in(2);
%delayed rectifier K current
ik = (v-CONST.EK_SOMA)*CONST.GKBAR_IDR_SOMA*in(3)^4;
%a slow non-inactivating K current
ikm = (v-CONST.EK_SOMA)*CONST.GKBAR_IM_SOMA*in(4);
%leak current
il = (v-CONST.EM_SOMA)*CONST.G_PAS_SOMA;

dv_soma = -((ina+ik+ikm+il)+(v-v_dend)/ra(1)-inj)/CONST.CM_SOMA;
%dv_soma = -((ik+ina+il+ikm)-inj)/CONST.CM_SOMA;

end