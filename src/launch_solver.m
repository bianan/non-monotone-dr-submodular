function [x_opt, opt_f, fs, runtime] = launch_solver(f, grad, param,method, max_iter);

switch method
    
    case 1, % nonmonotone_fw_variant
        [x_opt, opt_f, fs, runtime] ...
            = nonmonotone_fw_variant(f, grad, param,max_iter);
        
    case 2, % quadprogIP
        [x_opt, opt_f, fs, runtime] ...
            = quadprog_ip(f, grad, param,max_iter);
        
    case 3, % twophase FW
        [x_opt, opt_f, fs, runtime] ...
            = twophase_fw(f, grad, param, max_iter);
        
    case 4, % proj grad
        %   1/(k+1)
        [x_opt, opt_f, fs, runtime] ...
            = proj_grad(f, grad, param, max_iter);
        
end
end
