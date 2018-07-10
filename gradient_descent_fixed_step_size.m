function [new_eps, number_calls_obj, number_calls_adj] = ...
    gradient_descent_fixed_step_size(eps, last_eps, min_eps, max_eps, ...
    compute_objective, compute_gradient)

    step_size = 20;
    
    number_calls_obj = 0;
    [gradient_adj, number_calls_adj] = compute_gradient(eps);
    
    new_eps = update_eps(eps, eps - step_size * gradient_adj, min_eps, max_eps);  
end
