function SOL = set_xvec(SOL, x)
global NGATES
SOL(1:NGATES) = x;
end