function  [opt_x, opt_f, fs,  runtime ] ...
    = twophase_fw(f, grad, param,max_iter);
% fs:  function value in each iteration
% gamma_t:  stepsize

n = param.n;
tic;
%% first local
x_0 = zeros(n, 1);
epsilon = 1e-6;
[x_l, g_l, opt_fx, fxs] = non_convex_fw(f, grad, param,max_iter, x_0,...
    epsilon);

%% second local
new_ubar = param.ub - x_l;
z_0 = zeros(n, 1);
epsilon_p = 1e-6;
[z_lp, g_lp, opt_fz, fzs] = non_convex_fw(f, grad, param,max_iter, z_0,...
    epsilon_p,  new_ubar);

runtime = toc;
%%
if opt_fx >= opt_fz
    opt_x = x_l;
    opt_f = opt_fx;
    fs = fxs;
else
    opt_x = z_lp;
    opt_f = opt_fz;
    fs = fzs;
end
end



function  [opt_x, g_l, opt_f, fs] ...
    = non_convex_fw(f, grad, param, L, x_0, epsilon, new_ubar);

if nargin < 7
    new_ubar = param.ub;
end

fs = [];
gs = [];
% n = param.n;
% grad_0 = grad(x_0, param);

f_k = f(x_0, param);
fs = [fs, f_k];
% g_0 = dot(grad_0, d_0);

tic;
k=0;
x_k = x_0;
opt_x = x_k;
g_min = Inf;
while k <= L
    % calculate gradient
    grad_k = grad(x_k, param);
    
    v = LMO_non_convex(grad_k, param, new_ubar);   %
    d_k = v - x_k;
    g_k = dot(d_k, grad_k);
    if g_k < g_min
        g_min = g_k;
        opt_x = x_k;
        opt_f = fs(end);
    end
    
    gs = [gs, g_k];
    if g_k <= epsilon
        break
    end
    
    gamma_k = 2/(k+2);
    x_k = x_k + gamma_k*d_k;
    f_k = f(x_k, param);
    fs = [fs f_k];
    k = k +1;
    
end
g_l = g_min;
end



function vm = LMO_non_convex(grad_t, param, new_ubar)
% returned vm:
% currently using 'interior-point'
% x: current solution

% m = param.m;
n = param.n;
lb = param.lb;
ub = new_ubar;

b = param.b;
A = param.A;
Aeq = param.Aeq; beq = param.beq;
opts= param.opts;
x0 = [];
% solve lp
vm = linprog(-grad_t,A,b,Aeq,beq,lb,ub, x0, opts);
[s1, s2] =  size(vm);
if n ~= s1 || 1 ~= s2
  vm = zeros(n, 1); % in case of returning NaN solution, 
end
end
