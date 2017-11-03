function  [opt_x, opt_f, fs,  runtime] ...
    = quadprog_ip(f, grad, param,max_iter);
% fs:  function value in each iteration
% opt_x: returned solution
% opt_f: returned fun. value
% fs: empty
fs = [];
tic;

lb = param.lb;
ub = param.ub;
b = param.b;
A = param.A;
Aeq = param.Aeq; beq = param.beq;
opts= param.opts_quadprogIP;
H = -param.H;
h = -param.h;
% x0 = [];
[opt_x, opt_f, ~, ~] = quadprogIP(H,h, A,b,Aeq,beq,lb,ub, opts);
opt_f = -opt_f + param.c;
runtime = toc;
end
