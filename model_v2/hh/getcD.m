function [c D] = getcD(T, SOL, frac, cai)
% y' = c(x) + D(x)y
global CONST NGATES_SOMA
global CELL_PARAMS_RA  CELL_PARAMS_I
global CELL_PARAMS_CAREV0 CELL_PARAMS_CAREV1
inj = CELL_PARAMS_I;
ra = CELL_PARAMS_RA;

%% D matrix
% D-matrix for the soma
D(1,1) = -(CONST.GKBAR_INA_SOMA*SOL(1)^3*SOL(2)+...
    CONST.GKBAR_IDR_SOMA*SOL(3)^4+...
    CONST.GKBAR_IM_SOMA*SOL(4)+...
    CONST.G_PAS_SOMA+...
    1/ra(1))/CONST.CM_SOMA;

D(1,2) = 1/ra(1)/CONST.CM_SOMA;

% First dendrite
cm_d = CONST.CM_DEND;
gm_d = CONST.G_PAS_DEND;

D(2,1) = 1/ra(2)/cm_d;
D(2,2) = (-1/ra(2)-1/ra(3)-gm_d)/cm_d;
D(2,3) = 1/ra(3)/cm_d;

%for inner dendrites, all is the same...
a_d = -1/ra(3)-1/ra(3)-gm_d;
for k = 3:CONST.NSEG_DEND % fishy
    D(k,k-1) = 1/ra(3)/cm_d;
    D(k,k) = a_d/cm_d;
    D(k,k+1) = 1/ra(3)/cm_d;
end

% Last dendrite
D(end+1,end) = (-1/ra(3)-1/ra(4)-gm_d)/cm_d;
D(end,end-1) = 1/ra(3)/cm_d;
D(end,end+1) = 1/ra(4)/cm_d;

% Spine solver
D(end+1,end) = -(CONST.GKBAR_ICAL_SPINE*SOL(NGATES_SOMA+1)^2*SOL(NGATES_SOMA+2)+...
    CONST.GKBAR_IKA_SPINE*frac+...
    CONST.G_PAS_SPINE+...
    1/ra(5))/CONST.CM_SPINE;
D(end,end-1) = 1/ra(5)/CONST.CM_SPINE;

%% c-vector for the soma

if T(1) >= inj(1,2) && T(2) < inj(2,2)
    I = inj(2,1); %mA/cm2
else
    I = inj(1,1);
end

c = (CONST.ENA_SOMA*CONST.GKBAR_INA_SOMA*SOL(1)^3*SOL(2)+...
    CONST.EK_SOMA*CONST.GKBAR_IDR_SOMA*SOL(3)^4+...
    CONST.EK_SOMA*CONST.GKBAR_IM_SOMA*SOL(4)+...
    CONST.EM_SOMA*CONST.G_PAS_SOMA+...
    +I)/CONST.CM_SOMA;

% For dendrites, all is the same
c(1,2:CONST.NSEG_DEND+1) = CONST.EM_DEND*CONST.G_PAS_DEND/CONST.CM_DEND;

% Spine solver 
carev = CELL_PARAMS_CAREV1 - CELL_PARAMS_CAREV0*log(cai);
c(1,end+1) = (carev*CONST.GKBAR_ICAL_SPINE*SOL(NGATES_SOMA+1)^2*SOL(NGATES_SOMA+2)+...
            CONST.EK_SPINE*CONST.GKBAR_IKA_SPINE*frac+...
            CONST.EM_SPINE*CONST.G_PAS_SPINE)/CONST.CM_SPINE;
c = c';

end