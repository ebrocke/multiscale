function setup_cell_parameters()
global CONST
global CELL_PARAMS_CAREV0 CELL_PARAMS_CAREV1
global CELL_PARAMS_RA CELL_PARAMS_I CELL_PARAMS_I2M CELL_PARAMS_N2M

%CONST = load('modelconst.mat');

%coefficients for the ca reversial potential calculate
CELL_PARAMS_CAREV0 = (1e3) * (CONST.R*(CONST.CELSIUS+273.15))/(2*CONST.FARADAY);
CELL_PARAMS_CAREV1 = CELL_PARAMS_CAREV0*log(CONST.CAO);

%calculate coefficient to convert number of molecules [#] to Molar [M]
CELL_PARAMS_N2M = 1/(CONST.AVOGADRO*CONST.VOLL_SPINE);

%coefficients to convert ca current [mA/cm2] to Molar [M]
CELL_PARAMS_I2M = -1e-3*CELL_PARAMS_N2M*CONST.COLUMB*CONST.AREA_SPINE/2;

%axial resistance Ra = RA*4*l/(pi*d^2)
ra1 = CONST.RA_SOMA*4.0*CONST.L_SOMA/CONST.DIAM_SOMA^2; % Ohm
ra2 = CONST.RA_DEND*4.0*CONST.L_DEND/CONST.DIAM_DEND^2;
ra3 = CONST.RA_SPINE*4.0*CONST.L_SPINE^2/CONST.DIAM_SPINE;
%soma->dendrite
ra12 = CONST.L_SOMA*CONST.DIAM_SOMA*(ra1+ra2)/2; % Ohm.cm2
%dendrite->soma
ra21 = CONST.L_DEND*CONST.DIAM_DEND*(ra1+ra2)/2;
%dendrite->dendrite
ra22 = CONST.RA_DEND*4.0*CONST.L_DEND^2/CONST.DIAM_DEND; 
%dendrite->spine
ra23 = CONST.L_DEND*CONST.DIAM_DEND*(ra2+ra3)/2;
%spine->dendrite
ra32 = CONST.L_SPINE*CONST.DIAM_SPINE*(ra2+ra3)/2;

CELL_PARAMS_RA = [ra12 ra21 ra22 ra23 ra32];
CELL_PARAMS_I = [CONST.I_SOMA(:,1)./(pi*CONST.L_SOMA*CONST.DIAM_SOMA), CONST.I_SOMA(:,2)];

end