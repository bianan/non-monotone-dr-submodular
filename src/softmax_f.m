% calculate f
function f = softmax_f(x, param)
%  x must be vector  
L = param.L;
n = param.n;
I = eye(n);
f = log(det(diag(x)*(L-I) + I));
end
