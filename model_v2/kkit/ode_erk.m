function out=ode_erk(t, input, ca_clamp, ca_flux)
% globals used when clamping is applied
global CA_BASE T_CLAMP_START T_CLAMP_END
global EQ_REFS ERK_PARAMS_RATES

%if clamp ca, we fix CA
if ca_clamp && t>=T_CLAMP_START && t<T_CLAMP_END
    species = [ca_clamp; input];
elseif ca_clamp
    species = [CA_BASE; input];
else
    species = input;
end

%loop over the odes (integrate each specie)
out = zeros(length(EQ_REFS),1);
for i = 1 : length(EQ_REFS)
    out(i) = EQ_REFS{i}(species, ERK_PARAMS_RATES);
end
out(1) = out(1)+ca_flux;%species(1)-input(1);
end