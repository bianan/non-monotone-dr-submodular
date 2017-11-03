function  [opt_x, opt_f, fs,  runtime] ...
    = proj_grad(f, grad, param, max_iter);
% fs:  function value in each iteration
% stepsize: 1/(k +1) 
fs = [];

n = param.n;
x = zeros(n, 1);
f_t = f(x, param);
fs = [fs f_t];
tic;
k=0;
% grad_t = zeros(n,1);
while k <= max_iter
    % calculate gradient
    grad_t = grad(x, param); 
    gamma = 1/(1+k);
    x = x + gamma*grad_t; 
    % projection to poly 
    x = proj_polytope(x, param);
    k = k+1;
    f_t = f(x, param);
    fs = [fs f_t];
end
runtime = toc; 
opt_x = x; % possibly the max one 
opt_f = fs(end);
end


function y  = proj_polytope(x, param);

lb=param.lb;
ub=param.ub;
A = param.A;
b = param.b;
n = length(ub);

% formulate as QP 
H = eye(n);
h = -x;
opt_quad = optimoptions('quadprog','Display', 'off');
y = quadprog(H, h, A,b, [], [], lb, ub, [], opt_quad);
end
