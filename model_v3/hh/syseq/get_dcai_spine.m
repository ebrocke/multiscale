function dcai_spine=get_dcai_spine(in,ica)
global NGATES CONST
%cacum.mod NEURON example
depth = 1000/4; %nm
cai = in(NGATES+3);
dcai_spine = -ica/depth/CONST.FARADAY/2 * (1e10) + (CONST.CAI0 - cai)/CONST.TAU_DCAI;
end