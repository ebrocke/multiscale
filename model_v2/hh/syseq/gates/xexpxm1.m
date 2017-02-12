function y = xexpxm1(x)
if x > -37
   if abs(x) < 1e-8
      y = 1-0.5*x;
   else
      y = x/(exp(x)-1);
   end
else
   y = -tanh(x+37)+37;
end