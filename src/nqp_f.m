% calculate f
function f = nqp_f(x, param)
%  x must be vector  
H = param.H;
h = param.h;
f = 0.5*x'*H*x + h'*x + param.c;
end
