function out=ode_erk_mm(t, input, rates, eq_refs, ca_clamp, ...
    ca_flux_old, ca_flux_new, t_old, t_new)
% globals used when clamping is applied
global CA_BASE T_CLAMP_START T_CLAMP_END

%if clamp ca, we fix CA
if ca_clamp && t>=T_CLAMP_START && t<T_CLAMP_END
    species = [ca_clamp; input];
elseif ca_clamp
    species = [CA_BASE; input];
else
    species = input;
end      
%species(1)=species(1)+ca_flux;
%loop over the odes (integrate each specie)
out = zeros(length(eq_refs),1);
for i = 1 : length(eq_refs)
    out(i) = eq_refs{i}(species, rates);
end
%out(1) = out(1)+ca_flux_new;  %constant prediction
out(1) = out(1)+...
    ((t-t_old)/(t_new-t_old))*ca_flux_new+...
    ((t_new-t)/(t_new-t_old))*ca_flux_old;
end
