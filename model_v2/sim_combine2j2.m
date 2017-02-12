function sim_combine2j2()
%clear global; close all; clear all;
% dbstop warning
% dbstop error
addpath(genpath('hh'),genpath('kkit'), genpath('common'));

%global erk_iter hh_iter ode_erk_iter ode_cell_iter
global CONST NGATES %ICA T_PREV
% ICA=0;
% T_PREV=0;
info = true;
%%%% Set up time scales
T_SIM = 5; %s

%%%% Set up the parametrs
CONST = load('modelconst.mat');
setup_cell_parameters();
setup_erk_parameters();

ca_clamp = 0;
%%%% Initial sate variables
cell_init_vals = init_cell_variables(true);
erk_init_vals = init_erk_variables(true, ca_clamp);
erk_size = length(erk_init_vals);

%%%% Interface between two systems
[call_erk, erk_exch] = solve_erk_handle();
[call_cell, cell_exch] = solve_cell_handle();

tic;
options  = odeset('Reltol',1e-5,'Abstol',1e-12,...% 'MaxOrder', 2,...
    'Stats','on');
[t, y] = ode15s(@solve_sys_ode,[0 T_SIM], [erk_init_vals, cell_init_vals], options, ...
    erk_size, {call_erk erk_exch}, {call_cell cell_exch});

%semilogy(t(1:end-1),diff(t))
%[ts y15s]
runtime = toc
% y15s = deval(t,sol);
%y = y15s';
%dt = diff(ts);
%plot(t,y(:,erk_size+1));
plot_results(t,y,ca_clamp);
% Saving the data if wanted
if (info)
    dbclear all;
    mkdir(sprintf('../Data/%s',date));
    
    filename = 'ode15s';
    y=y';
    
    
    fn1 = sprintf('../Data/%s/%s.mat',date, filename);
    fn_y1 = sprintf('../Data/%s/%s_erk_ca_pmapk_ka.mat',date, filename);
    fn_y2 = sprintf('../Data/%s/%s_hh_gates.mat',date, filename);
    fn_y3 = sprintf('../Data/%s/%s_hh_v.mat',date, filename);
    
    %save(fn1);
    %save the workspace except y vector
    save(fn1,'-regexp','^(?!(y)$).')
    %save biochemical values
    erk_vars = y([1,6,9],:);
    save(fn_y1,'erk_vars');
    %save gate probabilities
    hh_gates = y((erk_size+1):(erk_size+NGATES),:);
    save(fn_y2,'hh_gates');%, 'ICA');
    %save voltages
    hh_v = y((erk_size+1+NGATES):end,:);
    save(fn_y3,'hh_v');
end
end

function out=solve_sys_ode(t, states, erk_size, erk_sys, cell_sys)
%global ICA ICA_SPINE T_PREV
% combining the two solvers to get a complete system, for stepping in time
erk_vals = states(1:erk_size);
cell_vals = states(erk_size+1:end);
erk = erk_sys{1};
cell = cell_sys{1};
out1 = erk(t, erk_vals, cell_vals, cell_sys{2});
out2 = cell(t, cell_vals, erk_vals, erk_sys{2});
% if t < T_PREV || abs(t-T_PREV)<1e-12
%     ICA(end)=ICA_SPINE;
% else
%     ICA(end+1)=ICA_SPINE;
% end
% T_PREV=t;
out = [out1; out2];
end
