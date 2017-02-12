function y = expp(x)
if x <= 2
   y = exp(x);
else
   y = exp(2)*(tanh(x-2)+1);
end
end