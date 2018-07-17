clc; clf; close all;

% Create the initial structure, a silicon waveguide in air.
eps_rows = 40;
eps_cols = 80;
num_eps_pixels = eps_rows * eps_cols;
eps0 = ones([eps_rows eps_cols]);

% Final simulation space size
simulation_width = 160;
simulation_height = 120;

% The maximum and minimum values of epsilon allowed in the structure
max_eps_value = 12.25;
min_eps_value = 1.0;
max_eps_array = 12.25 * ones(size(eps0));
min_eps_array = 1.0 * ones(size(eps0));

% The input mode coupling from and the output mode coupling to
input_mode = 1;
output_mode = 2;
design_frequency = 0.15;

% Run the rectangular waveguide completely through the middle of the
% coupler to begin with
eps0(:,30:50) = max_eps_value;

% Fill everything, except for a 2-cell boundary layer, with epsilon = 9.
% The 2-cell boundary layer is left unchanged because the boundary conditions
% use these cells.
eps0(3:end-2,3:end-2) = 9.0;

% Setup the design. Specifically, we want
% *   a design frequency of 0.15,
% *   to limit epsilon from 1 to 12.25, and
% *   to use the fundamental mode as an input on the left and
%     the second-order mode as the output mode on the right.
spec = setup(design_frequency, eps0, [min_eps_value max_eps_value], [input_mode output_mode]);

% Set the initial values of epsilon in the structure to the initial values.
% These will evolve in time over the course of the optimization.
eps = eps0;

upper_eps = max_eps_value * ones(num_eps_pixels, 1);
lower_eps = min_eps_value * ones(num_eps_pixels, 1);

% The desired output for Hz
c = spec.out.Hz;

% Compute the objective function for the optimization.  We are maximizing
% the power flowing through the end of the waveguide coupler in the desired
% output mode (unnormalized).  This requires knowledge of both Hz and Ey at
% the output to construct the Poynting vector.  Since we are looking for a
% minimization problem, we are actually computing the negative of this
% power flow.
compute_objective = @(input_Hz, input_Ey) objective_power(spec.out.Hz, spec.out.Ey, input_Hz, input_Ey);

% Partial derivative of the power-based objective function, g, with respect to
% the fields, x (Hz).  This will be the adjoint source.
g_x_partial = @(input_Hz, Ey_matrix, S_filter) x_partial_objective_power(spec.out.Hz, spec.out.Ey, input_Hz, Ey_matrix, S_filter);

min_z_value_unscaled = 1.0 / max_eps_value;
max_z_value_unscaled = 1.0 / min_eps_value;
z_offset = min_z_value_unscaled;
z_scale = (max_z_value_unscaled - min_z_value_unscaled);
z_scale_inverse = 1.0 / z_scale;
min_z_array_scaled = zeros(size(eps));
max_z_array_scaled = ones(size(eps));

get_z_from_eps = @(eps) ((1.0 ./ eps) - z_offset) * z_scale_inverse;
get_eps_from_z = @(z) 1.0 ./ ((z * z_scale) + z_offset);
deps_dz = @(z) (-1.0 ./ power(((z * z_scale) + z_offset), 2)) * z_scale;

% Setup gradient descent parameters.  Run for a fixed number of iterations
% with a fixed step size.
% This step size is likely far too big
% for the convergence test will want to do a running average to make sure convergence is stable
% or check something like not letting the objective value change too much
% between iterations as a way to check for convergence
step_size_max_unnormalized = 0.05;
step_size_min_unnormalized = 5e-5;
num_iter = 50;

% Should we also be moving a variable scaled to be between 0 and 1 (min/max
% epsilon values)

% We will follow the same strategy as the objective-first approach and 
z = get_z_from_eps(eps);

dims = size(eps);

num_objective_fn_calls = 0;
num_gradient_calls = 0;
gradient_norm = zeros(num_iter, 1);
avg_abs_grad = zeros(num_iter, 1);
obj_fn = zeros(num_iter, 1);

num_parameters = prod(size(eps(change_rows, change_cols)));

direction = zeros(num_parameters, num_iter);
reduced_gradient = zeros(num_parameters, num_iter);

type = 'pr';

for n = 1 : 1 : num_iter
    tic;
    eps = get_eps_from_z(z);
    ob1_plot_positive(dims, {'eps', eps}, {'z', z});
    [Ex_adj, Ey_adj, Hz_adj, gradient_adj] = ob1_fdfd_adj(spec.omega, eps, spec.in, g_x_partial, spec.bc);
    num_gradient_calls = num_gradient_calls + 1;
    adjoint_time = toc;
    
    gradient_adj_z = deps_dz(z) .* gradient_adj;
    
    current_objective = compute_objective(Hz_adj(end,:).', Ey_adj(end,:).');
    num_objective_fn_calls = num_objective_fn_calls + 1;
    
    reduced_gradient_z = gradient_adj_z(change_rows, change_cols);
    reduced_gradient(:, n) = reduced_gradient_z(:);
    
    last_idx = max(1, n - 1);
    
    if (strcmp(type, 'pr'))
    
        beta = ((reduced_gradient(:, n) - reduced_gradient(:, last_idx))' * reduced_gradient(:, n)) / ...
            (reduced_gradient(:, last_idx)' * reduced_gradient(:, last_idx));
        
    else
        
        beta = (reduced_gradient(:, n)' * reduced_gradient(:, n)) / ...
            (reduced_gradient(:, last_idx)' * reduced_gradient(:, last_idx));

    end
    
    direction(:, n) = -reduced_gradient(:, n) + beta * direction(:, last_idx);

    cur_direction = reshape(direction(:, n), [length(change_rows), length(change_cols)]);
    cur_direction_full = zeros(size(z));
    cur_direction_full(change_rows, change_cols) = cur_direction;
    get_max_direction_value = max(abs(cur_direction(:)));
    
    step_size = step_size_max_unnormalized / get_max_direction_value;
    step_size_min = step_size_min_unnormalized / get_max_direction_value;
    cur_z = update_z(z, z + step_size * gradient_adj_z, min_z_array_scaled, max_z_array_scaled);
    cur_eps = get_eps_from_z(cur_z);
    
    [~, Ey_check, Hz_check] = ob1_fdfd(spec.omega, cur_eps, spec.in, spec.bc);
    check_objective = compute_objective(Hz_check(end,:).', Ey_check(end,:).');

    while ((check_objective >= current_objective) && (0.5 * step_size > step_size_min))
        num_objective_fn_calls = num_objective_fn_calls + 1;
        step_size = 0.5 * step_size;
        cur_z = update_z(z, z + step_size * cur_direction_full, min_z_array_scaled, max_z_array_scaled);
        cur_eps = get_eps_from_z(cur_z);
        
        [~, Ey_check, Hz_check] = ob1_fdfd(spec.omega, cur_eps, spec.in, spec.bc);
        check_objective = compute_objective(Hz_check(end,:).', Ey_check(end,:).');
    end
       
    z = update_z(z, z + step_size * cur_direction_full, min_z_array_scaled, max_z_array_scaled);
    fprintf('Iteration #%d, Next Step Size = %f, Objective Function Value = %f, Computation Time (sec) = %f!\n', ...
        n, step_size, check_objective, adjoint_time);
    
    obj_fn(n) = check_objective;
    gradient_norm(n) = norm(gradient_adj_z(:), 2);
    avg_abs_grad(n) = mean(abs(gradient_adj_z(:)));
end

eps = get_eps_from_z(z);

% Simulate the result and analyze the success of the design process.
simulate(spec, eps, [simulation_width simulation_height]);

save_epsilon_filename = sprintf('adjoint_eps_%d_%d.csv', input_mode, output_mode);
csvwrite(save_epsilon_filename, eps);g

save_obj_fn_filename = sprintf('obj_fn_conj_grad_%s_%d_%d.csv', type, input_mode, output_mode);
save_grad_norm_filename = sprintf('grad_norm_conj_grad_%s_%d_%d.csv', type, input_mode, output_mode);
save_avg_abs_grad_filename = sprintf('obj_fn_avg_abs_grad_conj_grad_%s_%d_%d.csv', type, input_mode, output_mode);
save_num_obj_calls_filename = sprintf('num_obj_calls_conj_grad_%s_%d_%d.csv', type, input_mode, output_mode);
save_num_grad_calls_filename = sprintf('num_grad_calls_conj_grad_%s_%d_%d.csv', type, input_mode, output_mode);

csvwrite(save_obj_fn_filename, obj_fn);
csvwrite(save_grad_norm_filename, gradient_norm);
csvwrite(save_avg_abs_grad_filename, avg_abs_grad);
csvwrite(save_num_obj_calls_filename, num_objective_fn_calls);
csvwrite(save_num_grad_calls_filename, num_gradient_calls);

