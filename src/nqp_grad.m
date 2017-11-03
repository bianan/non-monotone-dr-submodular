% calculate g
function grad = nqp_grad(x, param)

% grad:  n*1
H = param.H;
h = param.h;
grad = x'*H + h';
grad = grad';

end
