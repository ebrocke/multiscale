%plots only Ca Ka and P_MAPK
function plot_results(t, y, ca_sti)
global CA_BASE T_CLAMP_START T_CLAMP_END
if ca_sti
    ca_plot = zeros(length(t),1);
    for i=1:length(t)
        if t(i) >=T_CLAMP_START && t(i) < T_CLAMP_END
            ca_plot(i)=ca_sti;
        else
            ca_plot(i)=CA_BASE;
        end
    end
    y = [ca_plot y];
end
plot(t,y(:,[1,6,9])*1e6), %in uM
legend('ca','p\_mapk','ka');
end