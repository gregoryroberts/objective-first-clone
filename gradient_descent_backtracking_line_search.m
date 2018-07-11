function [new_eps, number_calls_obj, number_calls_adj] = ...
    gradient_descent_backtracking_line_search(eps, last_eps, min_eps, max_eps, ...
    compute_objective, compute_gradient)

    beta = 0.5;
    alpha = 0.5;
    t = 100.0;
    min_step_size = 2.0;
        
    [row_dims, col_dims] = eps_dims_to_change(eps);

    [gradient_adj, number_calls_adj] = compute_gradient(eps);
    delta_x_full = -gradient_adj;

    [initial_objective, number_calls_obj] = compute_objective(eps);
    shift_eps = update_eps(eps, eps + t * delta_x_full, min_eps, max_eps);
    [objective, number_calls_obj] = compute_objective(shift_eps);
    
    select_vector_gradient = gradient_adj(row_dims, col_dims);
    vector_gradient = select_vector_gradient(:);

    delta_x_matrix = ...
        -(min(t * gradient_adj, max_eps - eps) .* (gradient_adj > 0) + ...
        max(t * gradient_adj, min_eps - eps) .* (gradient_adj <= 0));
    delta_x_matrix = delta_x_matrix(row_dims, col_dims);
    delta_x = delta_x_matrix(:);
    
%     h = 1e-8;

%     shift_eps = update_eps(eps, eps - h * gradient_adj, min_eps, max_eps);
%     [objective_test, number_calls_obj] = compute_objective(shift_eps);
% 
%     found_change = objective_test - initial_objective;
%     exp_change = (vector_gradient' * delta_x);
%     
%     fprintf('CHANGE = %e\nEXPECTED CHANGE = %e [%f]\n', ...
%         found_change, exp_change, exp_change / found_change);    
%     fprintf('Step size = %f!\n', t);

        
    while (((objective - initial_objective) > alpha*vector_gradient'*delta_x) && (t > min_step_size))
        t = beta * t;
        shift_eps = update_eps(eps, eps + t * delta_x_full, min_eps, max_eps);
        [objective, number_calls_obj] = compute_objective(shift_eps);

        delta_x_matrix = ...
            -(min(t * gradient_adj, max_eps - eps) .* (gradient_adj > 0) + ...
            max(t * gradient_adj, min_eps - eps) .* (gradient_adj <= 0));
        delta_x_matrix = delta_x_matrix(row_dims, col_dims);
        delta_x = delta_x_matrix(:);

        fprintf('%e versus %e (step_size = %e)\n', objective - initial_objective, alpha*t*vector_gradient'*delta_x, t)
    end

    t = max(t, min_step_size);

    new_eps = update_eps(eps, eps + t * delta_x_full, min_eps, max_eps);
end
