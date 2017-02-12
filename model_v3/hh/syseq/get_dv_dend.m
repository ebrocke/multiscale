function dv_dend = get_dv_dend(in, ra)
global CONST NGATES

v=in(end-CONST.NSEG_DEND+1:end);
v_soma  = in(NGATES+1);
v_spine = in(NGATES+2);

cm = CONST.CM_DEND;
gm = CONST.G_PAS_DEND;
em = CONST.EM_DEND;

dv_dend = zeros(CONST.NSEG_DEND,1);

if CONST.NSEG_DEND > 1
    dv_dend(1) = -((v(1)-em)*gm+(v(1)-v_soma)/ra(2)+(v(1)-v(2))/ra(3))/cm;
    dv_dend(end) = -((v(end)-em)*gm+(v(end)-v(end-1))/ra(3)+(v(end)-v_spine)/ra(4))/cm;
else
    dv_dend(1) = -((v(1)-em)*gm+(v(1)-v_soma)/ra(2))/cm;
end

for i=2:CONST.NSEG_DEND-1
   dv_dend(i) = -((v(i)-em)*gm+(2*v(i)-v(i-1)-v(i+1))/ra(3))/cm;
end
end