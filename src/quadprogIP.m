function [x_sol, fval_sol, time_sol, stats]= quadprogIP(H,f,A,b,Aeq,beq,LB,UB,options)
%% [x_sol,fval_sol,time_sol,stats] = quadprogIP(H,f,A,b,Aeq,beq,LB,UB,options)
%
% Authors: Wei Xia, Luis Zuluaga, Juan C. Vera
% quadprogIP is a non-convex quadratic program solver which solves problem with 
% the following form:
%
%    min      1/2*x'*H*x + f'*x
%    s.t.       A * x <= b
%             Aeq * x == beq
%             LB <= x <= UB
%
% --------------------------------------------------------------
% --> This code requires the Matlab interface to CPLEX 12.2  <--
% --> or later!                                              <--
% --------------------------------------------------------------
%
% Syntax:
%
%   x_sol = quadprogIP(H,f)
%   x_sol = quadprogIP(H,f,A,b)
%   x_sol = quadprogIP(H,f,A,b,Aeq,beq)
%   x_sol = quadprogIP(H,f,A,b,Aeq,beq,LB,UB)
%   x_sol = quadprogIP(H,f,A,b,Aeq,beq,LB,UB)
%   x_sol = quadprogIP(H,f,A,b,Aeq,beq,LB,UB,options)
%   [x_sol,fval_sol] = quadprogIP(H,f,...)
%   [x_sol,fval_sol,time_sol] = quadprogIP(H,f,...)
%   [x_sol,fval_sol,time_sol,stats] = quadprogIP(H,f,...)
%
% Input arguments:
%
% * H,f,A,b,Aeq,beq,LB,UB: identical to the corresponding input
%   arguments for MATLAB's QUADPROG function; see also the QP
%   formulation above
%
% * options: a structure with the following possible fields (defaults in
%   parentheses):
%
%   1) max_time (10000): the maximum amount of time QUADPROGBB
%      is allowed to run, in seconds
%
%   2) TolXInteger (1e-8): Specifies the amount by which an integer vari-
%      able can be different from an integer and still be considered 
%      feasible.
%
%   3) tol (1e-6): all-purpose numerical tolerance. For example, when
%      |LB(i) - UB(i)| < tol, we treat LB(i) as equal to UB(i), that is,
%      x(i) is fixed.
%
%   4) display (off): has the following different display levels
%         iter : display information at each iteration
%         final : display information of final iteration 
%         notify : display the final status
%
%   5) diagnostics (off):display the diagnostics information of Cplexmilp
%      solver's process.
%
%   6) max_iter (1000000): The maximum number of iterations allowed 
%
%   7) BranchStrategy 
%          (-1)Branch on variable with minimum infeasibility
%          (0) Automatic: let CPLEX choose variable to branch on; default
%          (1) Branch on variable with maximum infeasibility
%          (2) Branch based on pseudo costs
%          (3) Strong branching
%          (4) Branch based on pseudo reduced costs
%   8) Nodeselect
%          (0) Depth-first search
%          (1) Best-bound search; default
%          (2) Best-estimate search
%          (3) Alternative best-estimate search
%
% Output arguments:
%
% * x_sol,fval_sol: the solution and objective value of the QP; check
%   stat.status for the solution status, i.e., whether it is optimal
%
% * time_sol: time used by the branch-and-bound algorithm, in seconds
%
% * stats: a structure with more information:
%
%   1) time_pre: time spent on preprocessing
%
%   2) time_PB:  time spent on calculating primal bounds 
%
%   3) time_DB: time spent on calculating dual bounds
%
%   4) time_IP:  time spent on Integer branch-and-bound
%
%   5) nodes:    total number of nodes solved
%
%   6) status: final status of the solution
%
%      'opt_soln'  : optimal solution found
%      'time_limit': time limit specified by options.max_time was
%                    excedeeded
%      'inf_or_unb': the problem is infeasible or unbounded
%      
% 
% 


% start recording time 
tic;

% construct default option parameter object
defaultopt = struct(...
  'max_time'            ,10000,...
  'fix_var'             ,1e-8 ,...
  'tol'                 ,1e-6 ,...
  'constant'            ,0    ,...
  'Diagnostics'         ,'off',...
  'TolXInteger'         ,1e-12 ,...
  'nodeselect'          ,1    ,...
  'BranchStrategy'      ,1    ,...
  'display'             ,0    ,...
  'max_iter'            ,1000000);


% Argument checking
if nargin < 2 
  fprintf('Usage: \n');
  fprintf('[fval,x,time,stat] = QPIP(H,f,A,b,Aeq,beq,LB,UB)\n');
  if nargin > 0
     error('QPIP requires at least 2 input arguments.');
  end
  return
end

% Check if H is symmetric
[n1,n_vars] = size(H);
if n1 ~= n_vars
  error('H must be a square matrix!');
end

% if H not symmetric
H = .5*(H + H');
n = size(f,1);
old_H = H;
old_A = A;
old_Aeq = Aeq;

% Check dimension of H and f
if n ~= n1
  error('Dimensions of H and f are not consistent!');
end

% Set fixed variable indicator to 0
fix = 0;

% Assgin default option parameters if none given
if nargin < 9
  options = defaultopt;
else
  if isstruct(options)
    if ~isfield(options,'max_time')
      options.max_time = defaultopt.max_time;
    end
    if ~isfield(options,'tol')
      options.tol = defaultopt.tol;
    end
    if ~isfield(options,'Diagnostics')
      options.Diagnostics = defaultopt.Diagnostics;
    end
    if ~isfield(options,'max_iter')
      options.max_iter = defaultopt.max_iter;
    end
    if ~isfield(options,'TolXInteger')
      options.TolXInteger = defaultopt.TolXInteger;
    end
    if ~isfield(options,'BranchStrategy')
      options.BranchStrategy = defaultopt.BranchStrategy;
    end   
    if ~isfield(options,'display')
      options.display = defaultopt.display;
    end
    if ~isfield(options,'nodeselect')
      options.nodeselect = defaultopt.nodeselect;
    end
  else
    fprintf('The input argument options is not a struct!\n');
    fprintf('Overwritten with default options.\n\n');
    options = defaultopt;
  end
end

% Assign lower bound to negative infinity if none given
try
    if isempty(LB)
        LB = -inf*ones(n_vars,1);
    end
catch
    LB = -inf*ones(n_vars,1);
end

% Assign upper bound to infinity if none given
try
    if isempty(UB)
        UB = inf*ones(n_vars,1);
    end
catch
    UB = inf*ones(n_vars,1);
end

% Allow input without A and b
try
    A = A; b = b;
catch
    fprintf('No input of A or b detected. Using A = [] and b = [].\n\n');
    A = []; b = [];
end

% Allow input without Aeq and beq
try
    Aeq = Aeq; beq = beq;
catch
    fprintf('No input of Aeq or beq detected. Using Aeq = [] and beq = [].\n\n');
    Aeq = []; beq = [];
end

% Initiate problem status
stats.status = 'SOL NOT FOUND';

% End recording preprocessing
time_prep = toc;


% Calculate explicit primal bounds
x_var = size(H,1);
% Save original bounds for futre transformation
LB_o = LB;


% Convert problem to standard form
[LB,UB,time_PB] = prepbound(H,f,A,b,Aeq,beq,LB,UB,options);

% Find the indices of fixed variables
Fx_ind = find(abs(UB-LB)<options.fix_var);

% If all variables fixed, return solution
if length(Fx_ind) == n_vars
    x_sol = UB;
    fval_sol = 0.5*UB'*H*UB + f'*UB;
    time_sol = time_prep+time_PB;
    stats.total_time = time_sol;
    display(x_sol);
    display(fval_sol);
else
    % If some variables are fixed, reduce the problem
    if ~isempty(Fx_ind) 
        Fx = zeros(n_vars,1);
        Fx(Fx_ind) = UB(Fx_ind); 
        index = true(size(A,2),1);
        index(Fx_ind) = false;
        b = b - A(:,~index)*Fx(~index);
        A = A(:,index);
        beq = beq - Aeq(:,~index)*Fx(~index);
        Aeq = Aeq(:,index);
        temp_H = H(:,index);
        H_n = temp_H(index,:);
        index = +index;
        add_f = H*Fx;
        f = f + 2*add_f;
        index = logical(index);
        f = f(index);
        constant = sum(sum(H.*(Fx*Fx')));
        LB = LB(index);
        UB = UB(index);
        H = H_n;
        reduce = size(H,2);
        fix = 1;
    end
    
    
    [H,f,Aeq,beq,cons,LB,UB,time_refm] = standardform(H,f,A,b,Aeq,beq,LB,UB);



    % Append the box constraints to A
    %A = [A; eye(n_vars); -eye(n_vars)];
    %b = [b; UB; -LB];
    %num_vars = size(H,1);

    if issimplex(old_A,old_Aeq,n,LB)
        BigM = (norm(H,'fro')+norm(f,2))*n_vars*ones(n_vars,1);
        time_DB = 0;
    else
        % Compute the upper bounds for dual variables
        if isempty(old_A) & isempty(old_Aeq)
            [BigM,time_DB] = analytic_dual_bounds(H,f,Aeq,beq,LB,UB,x_var);
        else
            [BigM,time_DB] = dualbounds(H,f,Aeq,beq,LB,UB,options);
        end
    end

    % Prepare the problem formulation for integer program
    [f_IP,A_IP,b_IP,Aeq_IP,beq_IP,LB_IP,UB_IP,ctype_IP,time_PrepIP] = preIP(H,f,Aeq,beq,LB,UB,BigM);

    tic;
    lhs = [-inf * ones(size(A_IP,1),1);beq_IP];
    rhs = [b_IP;beq_IP];

    p = Cplex();                                                                    
    p.Model.sense = 'minimize';
    p.Model.obj   = f_IP;
    p.Model.lb    = LB_IP;
    p.Model.ub    = UB_IP;
    p.Model.ctype = ctype_IP;
    p.Model.A     = [A_IP;Aeq_IP];
    p.Model.lhs   = lhs;
    p.Model.rhs   = rhs;

    % Set options
    c_options = set_cplex_options(p, options);

    p.solve;
    b = toc;

    % Record calculation time
    time_IP = b;

    % Record integer branch and bound time
    tic;
    x = p.Solution.x;
    fval = p.Solution.objval;
    % Scale the solution to get the solution of original problem
    x_sol = x(1:n_vars)+LB_o;
    fval_sol = (1/2)*(fval)+cons;

    if fix == 1
        fval_sol = fval_sol + constant;
        Fx(index) = x_sol;
        x_sol = Fx;
    end
    
    %fval_sol = 0.5*x_sol'*old_H*x_sol + old_f'*x_sol;
    % Finish recording the of post calculation time
    time_post = toc;

    % Calculate the time for preprocessing and solving integer program
    time_Pre = time_refm+time_PB+time_DB+time_prep;
    time_sol = time_PrepIP+time_IP+time_post;

    stats.time_Pre = time_Pre;
    stats.time_IP = time_sol;
    stats.total_time = time_Pre+time_sol;

    % Update problem status
    if time_sol > options.max_time
        stats.status = 'time_limit';
    end
end
end



%% Auxillary functions

function [simplex] = issimplex(A,Aeq,n_vars,LB)
%% Check if the problem has form of a Simplex problem
    simplex = 0;
    if isempty(A) & ~isempty(Aeq) & size(Aeq,1) == 1 & ...
    sum(abs(Aeq - ones(1,n_vars))) == 0 & sum(abs(LB)) == 0
        simplex = 1;
    end
end

function [H,f,Aeq,beq,cons,LB,UB,time_refm] = standardform(H,f,A,b,Aeq,beq,LB,UB)
%% Transform problem into standard form and scale the varaibles to be between 0 and 1

tic;

if issimplex(A,Aeq,length(f),LB)
    cons = 0;
else
    m_ineq = size(A,1);
    n_var = size(H,1);

    n = 2*n_var+m_ineq;

    % Scale the coefficient matrices and bounds

    U = UB - LB;
    fn = [H*LB+ f; zeros(m_ineq+n_var,1)];
    cons = 0.5*LB'*H*LB + f'*LB;
    f = fn;

    if ~isempty(Aeq)
        beq = beq - Aeq*LB;    
    end

    if ~isempty(A)
        b = b - A*LB;    
        Us = b+abs(A)*U;
        Aeq = [A eye(m_ineq) zeros(m_ineq,n_var); eye(n_var) zeros(n_var,m_ineq) eye(n_var);Aeq zeros(size(Aeq,1),m_ineq+n_var)];
        beq = [b;U;beq];
        UB = [U;Us;U];
    else
        Aeq = [eye(n_var) zeros(n_var,m_ineq) eye(n_var);Aeq zeros(size(Aeq,1),m_ineq+n_var)];
        beq = [U;beq];
        UB = [U;U];
    end

    LB = zeros(n,1);


    H = [H zeros(n_var,n_var+m_ineq);zeros(n_var+m_ineq,2*n_var+m_ineq)];
end
time_refm = toc;

end


function [LB,UB,time_PB] = prepbound(H,f,A,b,Aeq,beq,LB,UB,options)
% Computes bounds for primal variables
tic;

n_vars = size(H,1);
ctype(1:n_vars) = 'C';
f_aux = zeros(n_vars,1);

c_options = cplex_options(options);

% Find Lower Bounds on original variables
I_lo = find(~isfinite(LB));
x0 = [];
for i=1:length(I_lo)
    f_aux(I_lo(i)) = 1;
    [x, fval, exitflag,output] = cplexmilp(f_aux,A,b,Aeq,beq,[],[],[],LB,UB,ctype,x0,c_options);
    if output.cplexstatus >= 103
        error('PROBLEM DOES NOT SATISFY BOUNDED ASSUMPTIONS');
    else
        LB(I_lo(i)) = fval;
    end;
    f_aux(I_lo(i)) = 0;
    x0 = x;
end;

%Find Upper Bounds on original variables
I_up = find(~isfinite(UB));
x0 = [];
for i=1:length(I_up)
    f_aux(I_up(i)) = -1;
    [x, fval, exitflag,output] = cplexmilp(f_aux,A,b,Aeq,beq,[],[],[],LB,UB,ctype,x0,c_options);    
    if output.cplexstatus >= 103
        error('PROBLEM DOES NOT SATISFY BOUNDED ASSUMPTIONS');
    else
        UB(I_up(i)) = -fval;
    end;
    f_aux(I_up(i)) = 0;
    x0 = x;
end;

time_PB = toc;

end


function [BigM,time_DB] = analytic_dual_bounds(H,f,A,b,LB,UB,n)
m_ineq = size(H,1);

BigM = (0.5*sum(max(abs(H)))*(2*max(UB)+1)+norm(f,1))*ones(m_ineq,1);

time_DB = toc;
end


function [BigM,time_DB] = dualbounds(H,f,Aeq,beq,LB,UB,options)
%Find dual variables bounds

BigM = zeros(size(H,1),1);

% Add lower and upper bounds to inequalities

n_vars = size(H,1);
%m_ineq = size(A,1);
m_eq = size(Aeq,1);

% Calculate bounds for dual variables using optimization

% variable order [Up_vec(X),x,y,lambda,rho]
% x are the original variables

r_vars = n_vars * n_vars;

if m_eq > 0
    Aeq_BD = [zeros(m_eq, r_vars) Aeq zeros(m_eq, m_eq + 2*n_vars)];
    beq_BD = beq;
else
    Aeq_BD = [];
    beq_BD = [];
end;

% Add Normal KKT


Aeq_BD = [Aeq_BD; zeros(n_vars, r_vars) H Aeq' -eye(n_vars) eye(n_vars)];
beq_BD = [beq_BD; -f];


% Add linearized KKT
H_lin = H(:);

Aeq_BD = [Aeq_BD; H_lin' f' beq' zeros(1,n_vars) UB'];

beq_BD = [beq_BD; 0];


% Construct Upper and Lower Bounds
S_UB = (UB*UB');
vec_S_UB = S_UB(:);

S_LB = (LB*LB');
vec_S_LB = S_LB(:);

LB_BD = [vec_S_LB; LB; -Inf*ones(m_eq,1); zeros(n_vars,1); zeros(n_vars,1)];
UB_BD = [vec_S_UB; UB; Inf*ones(m_eq,1);Inf*ones(n_vars,1);Inf*ones(n_vars,1)];

% Solve for lambda bounds
for i=1:n_vars
    f_BD = zeros(1,r_vars+n_vars+m_eq+2*n_vars);
    f_BD(n_vars + r_vars + m_eq + i) = -1;
    
    [x1, fval1, exitflag1,output1] = cplexlp(f_BD,[],[],Aeq_BD,beq_BD,LB_BD,UB_BD,[],options);
    if output1.cplexstatus >= 103 
        error('UUPS DUAL VARIABLES ARE UNBOUNDED');
    else
        if isempty(fval1)
            BigM(i) = abs(UB(i));
        else
            BigM(i) = -fval1;
        end
    end;
end;

time_DB = toc;
end



function [f_IP,A_IP,b_IP,Aeq_IP,beq_IP,LB_IP,UB_IP,ctype_IP,time_PrepIP] = preIP(H,f,Aeq,beq,LB,UB,BigM)
% Set up the integer program formulation for the problem

tic;    

% Save primal and dual vairable sizes
m_eq = size(Aeq,1);
n_vars = size(H,1);

%% variable order [x lambda mu z]
%% objective vector
f_IP = [f; zeros(n_vars,1); -beq; zeros(n_vars,1)];


%% Constraint Matrix

%Aeq = beq
Aeq_IP = [Aeq, sparse(m_eq,n_vars + m_eq + n_vars)];  beq_IP = beq;

%H x - lambda + Aeq' nu = -f 
Aeq_IP = [Aeq_IP; H -eye(n_vars) Aeq' zeros(n_vars, n_vars)]; beq_IP = [beq_IP; -f];

% x - zU <= 0
A_IP = [eye(n_vars), sparse(n_vars,n_vars + m_eq), -diag(UB)];  b_IP = zeros(n_vars,1);

% lambda + M z <= Me
A_IP = [A_IP; sparse(n_vars,n_vars) speye(n_vars) sparse(n_vars,m_eq) diag(BigM)]; b_IP = [b_IP; BigM];

%if issimplex(A,Aeq,n_vars,LB)
% -z_u - z_L <= -1
%    A_IP = [A_IP; sparse(n_vars,n_vars + m_ineq + m_eq) -speye(n_var) -speye(n_vars)]; b_IP = [b_IP; -ones(n_vars,1)];
%end

% Variable Upper and Lower bounds
LB_IP = [LB; zeros(n_vars,1); -Inf*ones(m_eq,1); zeros(n_vars,1)];
UB_IP = [UB; BigM; Inf*ones(m_eq,1); ones(n_vars,1)];


%% Integer variables
ctype_IP(1:n_vars + n_vars + m_eq) = 'C';
ctype_IP(n_vars + n_vars + m_eq + 1 : n_vars + n_vars + m_eq + n_vars) = 'B';

time_PrepIP = toc;
end



function [c_options] = cplex_options(options)
% Set options according to user's specification

c_options = cplexoptimset('cplex'); 
c_options.display = 'off';
c_options.diagnostics = 'off';
c_options.mip.strategy.variableselect = options.BranchStrategy;
c_options.timelimit = options.max_time;
c_options.TolFun = options.tol;
c_options.MaxTime = options.max_time;
c_options.mip.strategy.nodeselect = options.nodeselect;
c_options.mip.tolerence.mipgap = options.tol;
c_options.mip.tolerances.integrality = options.TolXInteger;
end

function [options] = set_cplex_options(p,options)
% Set options according to user's specification

%p.Param.display.Cur = 'off';
%p.Param.diagnostics.Cur = 'off';
p.Param.mip.strategy.variableselect.Cur = options.BranchStrategy;
p.Param.timelimit.Cur = options.max_time;
%p.Param.simplex.tolerances.optimality.Cur = options.tol;
p.Param.mip.strategy.nodeselect.Cur = options.nodeselect;
p.Param.mip.tolerances.mipgap.Cur = options.tol;
p.Param.mip.tolerances.integrality.Cur = options.TolXInteger;
end
