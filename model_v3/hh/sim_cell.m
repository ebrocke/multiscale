%usage:
%with_buff : boolean (whether there is one more ode for the ca buffer)
function sim_cell()
global CONST NGATES NGATES_SOMA

T_SIM = 2; %s
CONST = load('modelconst.mat');
%calculate axial resistance [Ohm.cm2] and current injection [mA/cm2]
setup_cell_parameters();

%get initial values for the state variables
%false = ~coupled
%true = coupled
inits = init_cell_variables(false); 

%init the solver
odeslv = @ode15s;
options  = odeset('Reltol',1e-6,'Abstol',1e-9, 'Stats', 'on'); 

%simulate
[t y] = odeslv(@ode_cell, [0 T_SIM], inits, options, -1, .00024);
%semilogy(t(1:end-1),diff(t))
plot(t, y(:,NGATES+1));
end