function SOL = set_yvec(SOL, y)
global NGATES
SOL([NGATES+1,NGATES+3:end,NGATES+2]) = y; %voltages in the order: soma,dend,spine
end