% 1. 2Ca+Raf<=>Active_Raf
% 2. Active_Raf+MAPK<=>AM=>P_MAPK+Active_Raf
% 3. P_MAPK+Phosphotase<=>PP=>MAPK+Phosphotase
% 4. P_MAPK+K_A<=>PK=>P_K_A+P_MAPK
% 5. P_K_A=>K_A
% 6. PKC+2AA<=>Act_PKC
% 7. AA<=>APC
% 8. P_MAPK+APC<=>PA=>AA+P_MAPK
% 9. MAPK+Act_PKC<=>MA=>P_MAPK+Act_PKC
% 10. Ca+pmca<=>pmcaCa=>pmca
% in [(Molar*s)^-1]
function setup_erk_parameters()
global ERK_PARAMS_RATES
ERK_PARAMS_RATES = [4e12 8.0 0.0 ;             %1
     0.025090663e8 40.0 10.0;   %2
     0.501831326e8 0.4 0.1;     %3
     0.050184337e8 40.0 10.0;   %4
     0.05 0.0 0.0;              %5
     1e12 2.0 0.0;              %6  
     0.2 0.01 0.0;              %7
     0.250918674e8 20.0 5.0;    %8
     0.050184337e8 4.0 1.0;     %9
     0.06e7 7.0 5.0];           %10
end