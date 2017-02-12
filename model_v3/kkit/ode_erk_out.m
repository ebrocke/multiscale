function  status = ode_erk_out(t, in,flag, rates, eq_refs, ca_clamp, ca_flux)
global COUNT inj
if strcmp(flag,'done')
    return
end
if  ~strcmp(flag,'init')
    if inj~=0
        COUNT=COUNT+1;
    end
    inj=0;
    
end
status=0;
end