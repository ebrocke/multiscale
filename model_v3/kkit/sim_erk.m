%usage:
% with_pump : boolean
% ca_clamp: double (0 if no clamping) otherwise 1e-6 M is a clamping value
% path:the path to save outputs (if no save needed leave empty string)
function sim_erk(with_pump, ca_clamp,saving)
%addpath('rre') %integration functions for each specie
% used when clamping is applied
global T_CLAMP_START T_CLAMP_END CA_BASE CONST
T_SIM = 20; %s
ca_clamp = 1e-6;
with_pump = true;
saving = true;
CONST = load('modelconst.mat');

%decalre reaction rates
setup_erk_parameters();

%init each specie
species = init_erk_variables(with_pump, ca_clamp);

%if we clamp let's fix CA
if ca_clamp>0
    CA_BASE = species(1);
    T_CLAMP_START = 1; %s
    T_CLAMP_END = 5; %s
    species = species(2:end);
end

%init the solver

% solver = @RK4;
%solver = @CrankNicSimp;
solver = @ode15s;

tic
%[TOUTi YOUTi] = solver(...
%    @ode_erk,[0 T_SIM],species,options, rates,eq_refs,ca_clamp, 0);

%specify MaxStep in order not to jump over the clamping interval

 options  = odeset('Reltol',1e-10,'Abstol',1e-12, 'MaxStep', 1, ...
      'Stats','on'  );
% %simulate
[t y] = ode15s(...
     @ode_erk,[0 T_SIM],species,options,ca_clamp, 0);
 
% %plot Ca P_MAPK KA
% y = deval(sol,TOUTi);

toc

if (saving)
    mkdir(date);
    filename = sprintf('%s/%sCa=%f.mat',date,func2str(solver),ca_clamp);
    save(filename)
end
figure(1); clf;hold on; grid on;
%plot Ca P_MAPK KA
%plot_results(t, y,ca_clamp);
plot(t,y(:,8))
% plot_results(TOUTi, y', ca_clamp);
%plot_results(TOUTi, YOUTi, ca_clamp);
end


