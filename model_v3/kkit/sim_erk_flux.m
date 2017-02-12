function sim_erk_flux()
addpath('rre') %integration functions for each specie
% used when clamping is applied
with_pump = true;
ca_clamp=0;
%get function handles to inegrate each specie
eq_refs = get_dsps_refs(with_pump);
%init reaction rate constants
rates = get_reac_rates();
%init each specie
species = get_species(with_pump);

%read flux from file
flux = dlmread('../flux');

%init the solver
odeslv = @ode15s;
%specify MaxStep in order not to jump over the clamping interval
options  = odeset('Reltol',1e-3,'Abstol',1e-9);%,'OutputFcn',@ode_erk_out );
y(1,:)=species;
t=0;
%simulate
for i=2:length(flux(:,1))
    t2 = flux(i,1);
    t1 = flux(i-1,1);
    inj = flux(i,2);%*(t2-t1);
   % y(end,1) = y(end,1)+inj;
    [t_ y_]=odeslv(@ode_erk,[t1 t2], y(end,:), options, rates, eq_refs, ca_clamp, inj);
    t=[t;t_];
    y=[y;y_];
end

%plot Ca P_MAPK KA
plot_results(t, y, ca_clamp);
end

