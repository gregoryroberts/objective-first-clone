function [performance_plot, number_calls_obj, number_calls_adj, final_eps] = evaluate_descent_method(...
    descent_method, ...
    c, omega, input, boundary_condition, eps_initial, min_eps_array, max_eps_array, ...
    g_x_partial, compute_objective, compute_gradient )

    clear objective_function;
    clear adjoint_gradient;
        
    number_calls_obj = 0;
    number_calls_adj = 0;
    
    num_iterations = 40;
    last_eps = eps_initial;
    
    eps = eps_initial;
    
    performance_plot = zeros(num_iterations, 1);
    
    calculate_gradient = @(eps) compute_gradient(omega, eps, input, g_x_partial, boundary_condition);
    calculate_objective = @(eps) compute_objective(c, omega, eps, input, boundary_condition);
    
    for n = 1 : 1 : num_iterations
        
        % For now, we enforce a projective method and enforce epsilon to
        % stay within its minimum and maximum values.
        
        [performance, number_calls_obj] = objective_function(c, omega, eps, input, boundary_condition, 'ignore_fn_call');
        performance_plot(n) = performance;
        
        [new_eps, number_calls_obj, number_calls_adj] = ...
            descent_method(eps, last_eps, min_eps_array, max_eps_array, ...
            calculate_objective, calculate_gradient);
        eps = new_eps;
        
        [get_gradient, number_calls_adj] = compute_gradient(omega, eps, input, g_x_partial, boundary_condition, 'ignore_fn_call');
        get_norm = norm(get_gradient(:), 2);
        
        fprintf('Iteration = %d, Current performance = %f, Gradient norm = %f!\n\n', n, performance, get_norm);

        last_eps = eps;
    end
    
    final_eps = eps;
    performance_plot = performance_plot(1:n);
end

