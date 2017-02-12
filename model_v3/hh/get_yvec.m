function y = get_yvec(SOL)
global NGATES
y = SOL([NGATES+1,NGATES+3:end,NGATES+2]); %voltages in the order: soma,dend,spine
end