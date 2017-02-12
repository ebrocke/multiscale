function dp_dt = get_dp_gates(in)
global REF_GATES NGATES NGATES_SOMA

v_soma = in(NGATES+1);
v_spine = in(NGATES+2);
dp_dt = zeros(NGATES,1);

for i = 1 : NGATES_SOMA
    %[alpha beta] = REF_GATES{i}(v_soma);
    %dp_dt(i) = alpha*(1-in(i))-beta*in(i);
    [inf tau] = REF_GATES{i}(v_soma);
    dp_dt(i) = (inf - in(i)) / tau;
end

for i = 1+NGATES_SOMA : NGATES
    [inf tau] = REF_GATES{i}(v_spine);
    dp_dt(i) = (inf - in(i)) / tau;

%    [alpha beta] = REF_GATES{i}(v_spine);
%    dp_dt(i) = alpha*(1-in(i))-beta*in(i);
end


end