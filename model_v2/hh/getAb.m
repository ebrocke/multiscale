function [A b] = getAb(SOL)
%% x' = A(y)x + b(y)
global NGATES NGATES_SOMA REF_GATES
v_soma = SOL(NGATES+1); 
v_spine = SOL(NGATES+2);
inf = zeros(NGATES,1);
tau = zeros(NGATES,1);

for i = 1 : NGATES_SOMA
    [inf(i) tau(i)] = REF_GATES{i}(v_soma);
end

for i = 1+NGATES_SOMA : NGATES
    [inf(i) tau(i)] = REF_GATES{i}(v_spine);
end

%calculate mat A
A = -diag(1./tau');

%calculate vector b
b = inf./tau;

end