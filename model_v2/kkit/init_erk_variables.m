%in [Molar]
function species = init_erk_variables(ca_clamp)
with_pump = true;
global APC_BASE KA_BASE CONST
global EQ_REFS 
%get function handles to inegrate each specie
EQ_REFS = get_dsps_refs(with_pump);
APC_BASE = 1e-6;
KA_BASE = 1e-6;
%calculate coefficient to convert number of molecules [#] to Molar [M]
N2M = 1/(CONST.AVOGADRO*CONST.VOLL_SPINE);

species = [0.2e-6, ... %ca 1.
    1e-6, ... %raf      2.  
    0, ...  %araf         3.
    1e-6, ...%mapk        4.
    0, ...%mapk_araf      5.
    0, ...%pmapk          6.
    0.5e-6, ... %phsph    7.
    0, ...%pmapk_phsph    8.
    1e-6, ...%ka        9.
    0, ...%ka_pmapk       10.
    0, ...%pka            11.
    1e-6,...%pkc          12.
    0, ... %apkc          13.
    0, ...%aa             14. 
    0, ...%apc_pmapk      15.
    0];%mapk_apkc         16.
if with_pump
    species = [species,...
    1950*N2M, ... %pmca    17.
    375*N2M ]; %pmcaCa     18.
end
if ca_clamp
    EQ_REFS = EQ_REFS(2:end);
end
end
