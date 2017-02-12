clear;
%% %% electrophysiology

%% SOMA compartment
EM_SOMA = -70.0; %mV
EK_SOMA = -90.0; %mV
ENA_SOMA = 50.0; %mV
G_PAS_SOMA = 1e-4; %S/cm2
GKBAR_INA_SOMA = 0.05; %S/cm2
GKBAR_IDR_SOMA = 0.005; % S/cm2
VTRAUB_SOMA = -63.0; %mV
GKBAR_IM_SOMA = 7e-5; % S/cm2
TAUMAX_IM_SOMA = 824.5; %ms
%capacitance
CM_SOMA = 1e-6; % F/cm2
%specific axial resistance
RA_SOMA = 35.4; %Ohm.cm

%current injection protocol for soma
%represented in the form
%[current amplitude [mA]] [end of stimulation [s]]
I_SOMA = [ 0.4*1e-6 1; 1.3*1e-6 6; 0.4*1e-6 45]; 
CAI_SOMA = 1.5-3; % M

%% DEND compartment
EM_DEND = EM_SOMA; %mV
G_PAS_DEND = 6e-4; %S/cm2
%capacitance
CM_DEND = CM_SOMA; % F/cm2
%specific axial resistance
RA_DEND = RA_SOMA; %Ohm.cm

%% SPINE compartment
EM_SPINE = EM_SOMA; %mV
G_PAS_SPINE = 1e-5; %S/cm2
%capacitance
CM_SPINE = CM_SOMA; % F/cm2
%specific axial resistance
RA_SPINE = RA_SOMA; %Ohm.cm

%ECA_SPINE = 132.457934164; %mV
EK_SPINE = EK_SOMA;
GKBAR_ICAL_SPINE = 3e-3; % S/cm2
GKBAR_IKA_SPINE = 0.00345; % S/cm2

%% %% common constants
CELSIUS	= 36; % degC
CAO 	= 2; %	mM
FARADAY = 9.6485309e4; %coulomb;
R =  8.31441; %Joules/Kelvin/mole
AVOGADRO = 6.02214e23; %mol-1
COLUMB = 6.242e18; %electrons/sec in 1 columb
%% %% morphology
L_SOMA = 96e-4; %cm
DIAM_SOMA = 96e-4; %cm

NSEG_DEND = 15;
L_DEND = 500e-4/NSEG_DEND;
DIAM_DEND = 1e-4;

L_SPINE = 1e-4;
DIAM_SPINE = 1e-4;
% cylinder 2*r*h*pi
AREA_SPINE = pi*L_SPINE*DIAM_SPINE;
%pi*(diam/2)^2*h
VOLL_SPINE = AREA_SPINE*DIAM_SPINE/4*1e-3; %L
%% %% save constants to a file
save('modelconst.mat'); 
clear;
