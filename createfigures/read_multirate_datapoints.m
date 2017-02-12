function [errs, calls] = read_multirate_datapoints(filename, signal, refsol_path)

format long
N = 5;
S = 3;
% number of species (S) x max number of datapoints
% species: vspine, pmapk, ca
errs = zeros(S,N);
calls=zeros(S,N);
PRINT = false;
PMAPK_ID = 2;
KA_ID = 3;
% voltage in the spine
VSPINE_ID = 2;
if (strcmp(signal,'fast'))
    CA_ID = 1; % calcium concentration in the erk
elseif (strcmp(signal,'slow'))
    CA_ID = 3; % calcium concentration in the cell
else
    'unknown signal'
    return
end


%% load reference solution
ref_sol = load(strcat(refsol_path,'ode15s.mat'),'t');
ref_v = load(strcat(refsol_path,'ode15s_hh_v.mat'),'hh_v');
ref_conc = load(strcat(refsol_path,'ode15s_erk_ca_pmapk_ka.mat'),'erk_vars');
ref_sol.vspine = ref_v.hh_v(VSPINE_ID,:);
ref_sol.pmapk = ref_conc.erk_vars(KA_ID,:);
if (strcmp(signal,'fast'))
    ref_sol.ca = ref_conc.erk_vars(CA_ID,:);
elseif (strcmp(signal,'slow'))
    ref_sol.ca = ref_v.hh_v(CA_ID,:);
end
%%

files = dir(strcat(filename,'*.mat'));
nfile = 0;
for file = files'
    M = 0; % number of solution points taken from the end time
    % number of species (lines) x number of solutions (used for average)
    solM = zeros(S, M+1);
    refsolM = zeros(S, M+1);
    
    nfile = nfile + 1;
    
    %% load solution
    filename = strcat('solutions/',file.name(1:end-4));
    sol = load(file.name,'t1','t2','stats','stats_erk','stats_cell');
    sol_v = load(strcat(filename,'_hh_v.mat'),'hh_v');
    sol_conc = load(strcat(filename,'_erk_ca_pmapk_ka.mat'),'erk_vars');
    %%
    
    if(any(strcmp('stats_erk',fieldnames(sol))) && ...
        any(strcmp('stats_cell',fieldnames(sol))))
        calls(1,nfile) = sol.stats_erk.n_ode_iter + sol.stats_erk.n_ode_numjac; %erk
        calls(2,nfile) = sol.stats_cell.n_ode_iter + sol.stats_cell.n_ode_numjac; %cell
        if (strcmp(signal,'fast'))
            calls(3, nfile) = calls(1,nfile);
        else
            calls(3, nfile) = calls(2,nfile);
        end 
    elseif(any(strcmp('stats',fieldnames(sol)))) % back compatibility
        % ode iter + numjac ode calls
        calls(1,nfile) = sol.stats(1,3) + sol.stats(1,2); %erk
        calls(2,nfile) = sol.stats(2,3) + sol.stats(2,2); %cell
        if (strcmp(signal,'fast'))
            calls(3, nfile) = calls(1,nfile);
        else
            calls(3, nfile) = calls(2,nfile);
        end
    else
        'Stat data is not available'
        return
    end
    
    
    solM(1,:) = sol_v.hh_v(VSPINE_ID,end-M:end);
    solM(2,:) = sol_conc.erk_vars(KA_ID,end-M:end);
    
    refsolM(1,:) = interp1(ref_sol.t, ref_sol.vspine, sol.t1(end-M:end), 'spline');
    refsolM(2,:) = interp1(ref_sol.t, ref_sol.pmapk, sol.t2(end-M:end), 'spline');
    if (strcmp(signal,'fast'))
        refsolM(3,:) = interp1(ref_sol.t, ref_sol.ca, sol.t2(end-M:end), 'spline');
        solM(3,:) = sol_conc.erk_vars(CA_ID,end-M:end);
    elseif (strcmp(signal,'slow'))
        refsolM(3,:) = interp1(ref_sol.t, ref_sol.ca, sol.t1(end-M:end), 'spline');
        solM(3,:) = sol_v.hh_v(CA_ID,end-M:end);
    end
    
    errs(:, nfile) = 100*mean(abs((refsolM - solM)./refsolM),2);
    
    if PRINT
        fprintf('%s: Relative error of (V) finale value = %g%%,%g%%\n',...
            filename, errs(1, nfile),errs(2, nfile))
    end
    
end
