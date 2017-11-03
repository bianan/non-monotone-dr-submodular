function  [opt_x, opt_f, fs,  runtime] ...
    = nonmonotone_fw_variant(f, grad, param,max_iter);
% fs:  function value in each iteration
% opt_x: returned solution
% opt_f: returned fun. value
fs = [];
gamma_cons = 1/max_iter;
% m = param.m;
n = param.n;
x = zeros(n, 1);
f_t = f(x, param);
fs = [fs f_t];
tic;
t=0;
% grad_t = zeros(n,1);
while t < 1
    % calculate gradient
    grad_t = grad(x, param); 
    vm = LMO_fw_variant(grad_t, param, x);   % 
    gamma = min(gamma_cons, 1 -t);
    x = x + gamma*vm;
    t = t+ gamma;
    f_t = f(x, param);
    fs = [fs f_t];
end
runtime = toc;
opt_x = x;
opt_f = fs(end);
end


function vm = LMO_fw_variant(grad_t, param, x)
% returned vm:  n*1
% x: current solution
% currently using 'interior-point'

% m = param.m;
n = param.n;
lb = param.lb;
ub = param.ub;
ub_new = ub - x;
b = param.b;
A = param.A;
Aeq = param.Aeq; beq = param.beq;
opts= param.opts;
x0 = [];
% solve lp
vm = linprog(-grad_t,A,b,Aeq,beq,lb,ub_new,x0, opts);
[s1, s2] =  size(vm);
if n ~= s1 || 1 ~= s2
  vm = zeros(n, 1); % in case of returning NaN solution
end
end

