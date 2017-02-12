function out = init_cell_variables(coupled)
global REF_GATES NGATES NGATES_SOMA
global CONST

REF_GATES = {@get_na_m, @get_na_h, @get_k_n, @get_k_p,...
        @get_ca_m, @get_ca_h};
if ~coupled
    REF_GATES{end+1} = @get_k_m;
    REF_GATES{end+1} = @get_k_h;
end
NGATES = length(REF_GATES);
NGATES_SOMA = 4;

inits = zeros(1, NGATES+2+CONST.NSEG_DEND);
%inits = zeros(NGATES+1,1);

volts = zeros(NGATES,1);
volts(1:NGATES_SOMA)=CONST.EM_SOMA;
volts(NGATES_SOMA+1:end)=CONST.EM_SPINE;

for i = 1 : NGATES
    %[alpha beta] = REF_GATES{i}(volts(i));
    %inits(i) = alpha/(alpha+beta);
    [inf ~] = REF_GATES{i}(volts(i));
    inits(i) = inf;
end

dend(1:CONST.NSEG_DEND) = CONST.EM_DEND;  

%inits(NGATES+1:end) = [CONST.EM_SOMA];
inits(NGATES+1:end) = [CONST.EM_SOMA, CONST.EM_SPINE, dend];
out = inits;
%end
end