function [ode, exch] = solve_erk_handle()
global KA_BASE 
ode = @solve_erk;
exch = @get_exch;

    function out = solve_erk(t, erk_vals, cell_vals, get_exch_vals)
        ca_flux = get_exch_vals(cell_vals,erk_vals(1));
        %ca_flux = pproc_hdl(t, ca_flux);
        out = ode_erk(t, erk_vals, 0, ca_flux);
    end
    function [frac, cai] = get_exch(vals) % in SI
        frac = vals(9)/KA_BASE;
        cai = vals(1);
    end
end