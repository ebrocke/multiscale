function [ode, exch] = solve_cell_handle()
global NGATES NGATES_SOMA CELL_PARAMS_I2M %ICA_SPINE
ode = @solve_cell;
exch = @get_exch;

    function  out = solve_cell(t, cell_vals, erk_vals, get_exch_vals)%, pproc_hdl)
            [frac cai] =  get_exch_vals(erk_vals); %in SI
            %frac = 1;
            %cai = 0.2e-6;
            out = ode_cell(t, cell_vals, frac, cai*1e3);
    end
    function ca_flux = get_exch(vals, cai)% in SI
        %cai=0.2e-6;
        if size(vals,2)>1
            if size(vals,2)==2
                m = 0.5*( 3*vals(NGATES_SOMA+1,end) - vals(NGATES_SOMA+1,end-1));
                h = 0.5*( 3*vals(NGATES_SOMA+2,end) - vals(NGATES_SOMA+2,end-1));
            else
                a1 = 15/8; a2 = -5/4; a3 = 3/8;
                m = a1*vals(NGATES_SOMA+1,end)+a2*vals(NGATES_SOMA+1,end-1)+a3*vals(NGATES_SOMA+1,end-2);
                h = a1*vals(NGATES_SOMA+2,end)+a2*vals(NGATES_SOMA+2,end-1)+a3*vals(NGATES_SOMA+2,end-2);
            end

            %m,h linear extrapolation half a step forward
            %m = 0.5*( 3*vals(NGATES_SOMA+1,end) - vals(NGATES_SOMA+1,end-1));
            %h = 0.5*( 3*vals(NGATES_SOMA+2,end) - vals(NGATES_SOMA+2,end-1));
        else
            m = vals(NGATES_SOMA+1);
            h = vals(NGATES_SOMA+2);
        end
        %get_ica_spine returns [mA/cm2]
        ical =  get_ica_spine(vals(NGATES+2),... %V in the spine
            m, h, ...%m, h gates
            cai*1e3);
        
        ca_flux = CELL_PARAMS_I2M*ical;
        
        %ICA_SPINE=ca_flux;
        
        
    end
end