close all;
clear all;
seed = 0; % fix seed to have consistent outcome 
rng(seed,'twister');
addpath('src/');
% add path of CPLEX, e.g., on my machine, it is
% '/opt/ibm/ILOG/CPLEX_Studio127/cplex/matlab/x86-64_linux/'
% replace the following line with the path on your machine 
addpath('/opt/ibm/ILOG/CPLEX_Studio127/cplex/matlab/x86-64_linux/');
%% setup
% different way of generating data
data_id = 1
% 1: quadratic: uniform distribution
% 2: quadratic: exponential distribution
% 3: softmax: polytope constraints with uniform distribution
% 4: softmax: polytope constraints with exponential distribution
data_name= {'quad_uniform', 'quad_exp' ,'softmax_uniform', 'softmax_exp'}; %
%
solver_names = {...
    'non-monotone Frank-Wolfe variant',...      % 1
    'quadprogIP', ...                           % 2
    'two-phase Frank-Wolfe',...                 % 3
    'ProjGrad (1/(k+1))',...                    % 4
    };
nm_names=length(solver_names);
result_path = 'results/';
if ~ exist(result_path)
    mkdir(result_path);
end
%
if 1 == data_id | 2 == data_id % for quadratic functions, use global-optimal solver quadprogIP
    solver_list = [1 2 3 4]
else
    solver_list = [1 3 4]
end
% three different ways to set #constraints:
% m = floor(0.5n), m = n, m = floor(1.5n)
mn_type_names = {'m-half-n', 'm-n', 'm-onehalf-n'};
nm_exps = 20;   % set # repeated experiments
K = 5;  % number of trials for different dimensionalities
nbase = 8;
n_stepsize = 2;
ns = [];
for i = 1:K
    n_tmp = nbase + (i-1)*n_stepsize;
    ns = [ns  n_tmp];
end

nm_solver = length(solver_list);
max_iter = 100;

%% run 
for mn_type = 1:3
    
    % generate subfix of the results name
    subfix = [mn_type_names{mn_type} '-n_exp' int2str(nm_exps) '-seed' int2str(seed)];
    
    
    for idx_exp = 1:nm_exps % repeated experiments
        
        for i=1:K  % i:  param idx
            n_tmp = ns(i);
            % generate random data
            [f, grad, param] = gen_data(data_id, n_tmp, mn_type);
            fprintf('experiment_id %d: data_name-%s, m: %d, n: %d\n', ...
                idx_exp, data_name{data_id}, param.m, n_tmp);
            
            for t = solver_list
                [opt_x, opt_f, fs,  runtime] = launch_solver(f, grad,  param, t, max_iter);
                results{idx_exp, i, t}.opt_f = opt_f; % optimal function value
                results{idx_exp, i, t}.fs = fs;  % function values in the iterations 
                results{idx_exp, i, t}.opt_x  = opt_x; % solution returned 
                results{idx_exp, i, t}.runtime = runtime;
            end
            %
        end
    end
    
    file_name = [result_path data_name{data_id} '-' subfix];
    save(file_name, 'results');
    
    %% plot figures
    fig_scale = 1;
    fWidth=400*fig_scale;
    fHeight=240*fig_scale;
    plot_opt = {'--db', ':^r', '--sm',...
        '-.vk'};  % line config. for the solvers
    close all
    
    hFig = figure;
    set(hFig, 'Units', 'points');
    set( hFig, 'Position', [0 0 fWidth fHeight]);
    set(hFig,'PaperPositionMode','auto');
    set(hFig, 'PaperUnits','points', 'PaperSize', [fWidth fHeight],...
        'PaperPosition', [0 0 fWidth fHeight]);
    set(hFig, 'Name', data_name{data_id});
    
    %  plot objective w.r.t #dimensionalities for softmax
    if 3 == data_id | 4 == data_id
        opt_fs = zeros(nm_exps, K,nm_names);
        for id = 1:nm_exps
            for i = 1:K
                for t = solver_list
                    opt_fs(id, i,t) =results{id, i,t}.opt_f;
                end
            end
        end
        opt_mean =squeeze( mean(opt_fs, 1));
        opt_std = squeeze( std(opt_fs, 1, 1));
        
        opt_max = max(max(opt_mean(:, solver_list)));
        opt_min = min(min(opt_mean(:, solver_list)));
        
        
        % ----------------------------------------------
        hands = [];
        for t = solver_list(2:end) % plot
            hi = plot(ns, opt_mean(:, t), plot_opt{t},'linewidth',2);
            hands = [hands hi];
            hold on;
            errorbar(ns, opt_mean(:, t), opt_std(:, t), plot_opt{t},'linewidth',1)
            hold on;
        end
        t = solver_list(1);
        hi = plot(ns, opt_mean(:, t),  plot_opt{t},'linewidth',2*fig_scale);
        hands = [hands hi];
        hold on;
        errorbar(ns, opt_mean(:, t), opt_std(:, t), plot_opt{t},'linewidth',1*fig_scale)
        hold off;
        
        legend(hands, solver_names{[solver_list(2:end), solver_list(1)]}, ...
            'Location','northoutside');
        legend('boxoff');
        
        set(gca,'fontsize',14*fig_scale)
        % set the proper y limit
        %         axis([ns(1)-1 ns(end)+1   opt_min opt_max])
        axis([ns(1)-1 ns(end)+1   -Inf 0.21])
        %         axis([ns(1)-1 ns(end)+1   opt_min - 0.01 opt_max + 0.015])
        xlabel('Dimensionality');
        ylabel('Function value');
        fig_name = [result_path, data_name{data_id}, '_', subfix];
        saveas(hFig, fig_name, 'pdf')
        
        % plot the approx. ratio for quadratic functions
    elseif 1 == data_id | 2 == data_id
        opt_solver_idx = 2; % solver index of the quadprogIP
        opt_fs = zeros(nm_exps, K,nm_names);
        for id = 1:nm_exps
            for i = 1:K
                for t = solver_list
                    opt_fs(id, i,t) =results{id, i,t}.opt_f;
                end
            end
        end
        
        opt_fs_simple = opt_fs(:,:,solver_list);
        solvers_no_opt = solver_list([1:opt_solver_idx-1, opt_solver_idx+1:end]);
        opt_fs_simple_no_opt = opt_fs_simple(:,:,solvers_no_opt);
        
        for id = 1:nm_exps
            for i = 1:K
                opt_fs_simple_no_opt(id,i,:) =   ...
                    opt_fs_simple_no_opt(id,i,:)./opt_fs_simple(id,i,opt_solver_idx);
            end
        end
        
        
        opt_mean =squeeze(mean(opt_fs_simple_no_opt, 1));
        opt_std = squeeze( std(opt_fs_simple_no_opt, 1, 1));
        opt_max = max(max(opt_mean));
        opt_min = min(min(opt_mean));
        
        hands = [];
        for t = 2:length(solvers_no_opt) % plot
            hi = plot(ns, opt_mean(:, t), plot_opt{t+1},'linewidth',2);
            hands = [hands hi];
            hold on;
            errorbar(ns, opt_mean(:, t), opt_std(:, t), plot_opt{t+1},'linewidth',1)
            hold on;
        end
        t = 1;
        hi = plot(ns, opt_mean(:, t),  plot_opt{t},'linewidth',2*fig_scale);
        hands = [hands hi];
        hold on;
        errorbar(ns, opt_mean(:, t), opt_std(:, t), plot_opt{t},'linewidth',1*fig_scale)
        hold off;
        
        legend(hands, solver_names{[solvers_no_opt(2:end), solvers_no_opt(1)]}, ...
            'Location','northoutside');
        legend('boxoff');
        % set(gca,'Interpreter','latex');
        set(gca,'fontsize',14*fig_scale)
        % axis([ns(1)-1 ns(end)+1   opt_min-0.01 1])
        axis([ns(1)-1 ns(end)+1   opt_min 1])
        % axis([ns(1)-1 ns(end)+1   -Inf Inf])
        xlabel('Dimensionality');
        ylabel('Approx. ratio');
        fig_name = [result_path, 'approx-ratio_', data_name{data_id}, '_', subfix];
        saveas(hFig, fig_name, 'pdf')
        
    end
end % mn_type
