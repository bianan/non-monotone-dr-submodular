function [f,grad, param] = gen_data(dataset,n, mn_type);

if 1 == mn_type
    m = floor(n/2);
elseif 2 == mn_type
    m = n;
elseif  3 == mn_type
    m = floor(1.5*n);
else
    disp('wrong value of mn_type');
end

switch dataset
    case 1, %  quadratic data, uniform distribution
        
        alpha = 0.1;
        nu = 0.01;
        H = -1*rand(n);  % each entry between [-1, 0]
        H = symmetry(H);
        ub = zeros(n, 1);
        lb = zeros(n, 1);
        A = nu+rand(m, n);
        b = 1*ones(m, 1);
        
        % set ub to be the tightest bound
        for j = 1:n
            ub(j) = min(b./A(:,j));
        end
        h = - 0.2*H'*ub;  %  non-monotone gradient
        
        param.m = m;
        param.n = n;
        param.A = A;
        param.b = b;
        param.H = H;
        param.h = h;
        param.ub = ub;     %
        param.lb = lb;
        
        param.Aeq = []; param.beq = [];
        param.opts = optimoptions('linprog','Algorithm','interior-point', 'Display', 'off');
        param.opts_quadprogIP = struct(...
            'max_time'            ,10000,...
            'fix_var'             ,1e-8 ,...
            'tol'                 ,1e-6 ,...
            'constant'            ,0    ,...
            'Diagnostics'         ,'off',...
            'TolXInteger'         ,1e-12 ,...
            'nodeselect'          ,1    ,...
            'BranchStrategy'      ,1    ,...
            'display'             , 2 );
        
        % find the global minimium by quadprogIP:
        [~, f_min, ~, ~] = quadprogIP(H,h, A,b,param.Aeq,param.beq,lb,ub, param.opts_quadprogIP);
        c = -f_min + alpha*abs(f_min);
        param.c = c;
        
        f = @nqp_f;
        grad = @nqp_grad;
        
    case 2, %  quadratic  data with  exp. distribution
        
        
        
        alpha = 0.1;
        nu = 0.01;
        lambdaH = 1;
        lambdaA = 0.25;
        
        H = -1*exprnd(1/lambdaH, n, n);  % each entry
        H = symmetry(H);
        ub = zeros(n, 1);
        lb = zeros(n, 1);
        %
        A = nu +exprnd(1/lambdaA ,m, n);
        b = 1*ones(m, 1);
        
        % set ub to be the tightest bound
        for j = 1:n
            ub(j) = min(b./A(:,j));
        end
        
        h = - 0.2*H'*ub;  % non-monotone gradient
        
        param.m = m;
        param.n = n;
        param.A = A;
        param.b = b;
        param.H = H;
        param.h = h;
        param.ub = ub;  %
        param.lb = lb;
        
        param.Aeq = []; param.beq = [];
        param.opts = optimoptions('linprog','Algorithm','interior-point', 'Display', 'off');
        param.opts_quadprogIP = struct(...
            'max_time'            ,10000,...
            'fix_var'             ,1e-8 ,...
            'tol'                 ,1e-6 ,...
            'constant'            ,0    ,...
            'Diagnostics'         ,'off',...
            'TolXInteger'         ,1e-12 ,...
            'nodeselect'          ,1    ,...
            'BranchStrategy'      ,1    ,...
            'display'             , 2 );
        
        % find the global minimium by quadprogIP:
        % x0 = [];
        [~, f_min, ~, ~] = quadprogIP(H,h,A,b,param.Aeq,param.beq,lb,ub,param.opts_quadprogIP);
        c = -f_min + alpha*abs(f_min);
        param.c = c;
        
        f = @nqp_f;
        grad = @nqp_grad;
        
        
    case 3, % softmax,  uniform polytope constraints
        
        
        
        nu = 0.01;
        largest_eigen = 1.5;
        eigens = rand(n,1)*largest_eigen;
        % specified eigens
        D = diag(eigens); % Random negative eigenvalues
        V = orth(rand(n)); % Random unitary matrix
        L = V*D*V'; %
        
        lb = zeros(n,1);
        ub = zeros(n,1);
        
        A = nu + rand(m, n);
        b = 2*ones(m, 1);
        
        for i = 1:n
            ub(i) = min(b./A(:,i));
        end
        
        param.m = m;
        param.n = n;
        param.A = A;
        param.b = b;
        param.L = L;
        param.ub = ub;     %
        param.lb = lb;
        
        param.Aeq = []; param.beq = [];
        param.opts = optimoptions('linprog','Algorithm','interior-point', 'Display', 'off');
        
        f = @softmax_f;
        grad = @softmax_grad;
        
        
    case 4, % softmax objective,   polytope constraints with exponential distribution
        
        
        nu = 0.01;
        lambdaA = 1;
        largest_eigen = 1.5;
        %         eigens = linspace(0, largest_eigen, n);
        eigens = rand(n,1)*largest_eigen;
        % specified eigens
        D = diag(eigens); % Random negative eigenvalues
        V = orth(rand(n)); % Random unitary matrix
        L = V*D*V'; %
        
        
        lb = zeros(n,1);
        ub = zeros(n,1);
        
        A = nu +exprnd(1/lambdaA ,m, n);
        b = 2*ones(m, 1);
        
        for i = 1:n
            ub(i) = min(b./A(:,i));
        end
        
        param.m = m;
        param.n = n;
        param.A = A;
        param.b = b;
        param.L = L;
        param.ub = ub;     %
        param.lb = lb;
        
        param.Aeq = []; param.beq = [];
        param.opts = optimoptions('linprog','Algorithm','interior-point', 'Display', 'off');
        
        f = @softmax_f;
        grad = @softmax_grad;
        
end
end
