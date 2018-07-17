function [val, gradient, obj_fn_iter, grad_norm_iter, number_calls] = value_and_gradient_nlopt(omega, z, input, g_x_partial, bc, compute_objective, change_rows, change_cols, full_z, deps_dz, get_eps_from_z)
    persistent num_calls;
    persistent obj_fn;
    persistent grad_norm;
    if isempty(num_calls)
        num_calls = 0;
        obj_fn = [];
        grad_norm = [];
    end
    num_calls = num_calls + 1;
    
    number_calls = num_calls;
    
    reshape_z = reshape(z, [length(change_rows), length(change_cols)]);
    full_z(change_rows, change_cols) = reshape_z;
    
    full_eps = get_eps_from_z(full_z);
    
    tic;
    [~, Ey, Hz, gradient_full] = ob1_fdfd_adj(omega, full_eps, input, g_x_partial, bc);
    adj_time = toc;

    gradient_full_z = deps_dz(full_z) .* gradient_full;
    
    reduce_gradient = gradient_full_z(change_rows, change_cols);
    gradient = reduce_gradient(:);
    
    val = real(compute_objective(Hz(end,:).', Ey(end,:).'));
    
    obj_fn(num_calls) = val;
    grad_norm(num_calls) = norm(gradient_full_z(:), 2);
    
    obj_fn_iter = obj_fn;
    grad_norm_iter = grad_norm;
      
    % sometimes imag part of val can be really small so we use real part of
    % it even though it should just be a magnitude (order of ops must
    % affect this).
    fprintf('Current objective on call number %d is %f.  Adj took %f seconds \n', ...
        number_calls, val, adj_time);
end

