function out = ode_cell(t, input, frac, cai)
global CELL_PARAMS_RA CELL_PARAMS_I 

inj = CELL_PARAMS_I;
%frac = 1;
%cai = 0.2e-6*1e3;
%stimuli
%strong assumption of the current protocol
if t>= inj(1,2) && t < inj(2,2)
    I = inj(2,1); %mA/cm2
else 
    I = inj(1,1);
end

%calculate d_dt
dv_soma = get_dv_soma(input, CELL_PARAMS_RA, I);
dv_dend = get_dv_dend(input, CELL_PARAMS_RA);
dv_spine = get_dv_spine(input, CELL_PARAMS_RA, frac, cai);
dp_gates = get_dp_gates(input);

%out = [dp_gates; dv_soma]; 

out = [dp_gates; dv_soma; dv_spine; dv_dend];
end