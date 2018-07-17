clc; clf; close all;

% Try out the nlopt MMA

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
max_eps_array = max_eps_value * ones(size(eps0));
min_eps_array = min_eps_value * ones(size(eps0));

% The input mode coupling from and the output mode coupling to
input_mode = 1;
output_mode = 3;
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
% the overlap of the output mode from the structure with the desired output
% waveguide mode specified.  As a minimization problem, this is equivalent
% ot minimizing the negative of that same overlap integral.
% compute_objective = @(input_Hz) (-(c'*input_Hz)*conj(c'*input_Hz));
compute_objective = @(input_Hz, input_Ey) objective_power(spec.out.Hz, spec.out.Ey, input_Hz, input_Ey);

% Partial derivative of the objective function, g, with respect to the
% field variable input_Hz
% g_x_partial = @(input_Hz, S_filter) -c'*S_filter*conj(c'*(S_filter*input_Hz));
g_x_partial = @(input_Hz, Ey_matrix, S_filter) x_partial_objective_power(spec.out.Hz, spec.out.Ey, input_Hz, Ey_matrix, S_filter);

% Setup gradient descent parameters.  Run for a fixed number of iterations
% with a fixed step size.
% This step size is likely far too big
% for the convergence test will want to do a running average to make sure convergence is stable
% or check something like not letting the objective value change too much
% between iterations as a way to check for convergence


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

z = get_z_from_eps(eps);
z0 = get_z_from_eps(eps0);

[change_rows, change_cols] = eps_dims_to_change(eps);

my_func = @(input_z) value_and_gradient_nlopt(spec.omega, input_z, spec.in, ...
    g_x_partial, spec.bc, compute_objective, change_rows, change_cols, ...
    z0, deps_dz, get_eps_from_z);

clear value_and_gradient_nlopt;


opt.maxeval = 50;
opt.stopval = stop_obj_val;
opt.algorithm = NLOPT_LD_MMA;
opt.min_objective = my_func;
trim_z = z(change_rows, change_cols);
opt.lower_bounds = zeros(size(trim_z(:)));
opt.upper_bounds = ones(size(trim_z(:)));
[z_opt, fmin, retcode] = nlopt_optimize(opt, trim_z(:));

eps_opt = get_eps_from_z(z_opt);

[~, ~, obj_fn, grad_norm, final_number_calls] = my_func(z_opt);
final_number_calls = final_number_calls - 1;

fprintf('The number of calls needed to get to %f was %d!\n', stop_obj_val, final_number_calls);

reshape_eps = reshape(eps_opt, [length(change_rows), length(change_cols)]);
full_eps = eps;
full_eps(change_rows, change_cols) = reshape_eps;

% Simulate the result and analyze the success of the design process.
simulate(spec, full_eps, [simulation_width simulation_height]);

save_epsilon_filename = sprintf('adjoint_eps_%d_%d.csv', input_mode, output_mode);
csvwrite(save_epsilon_filename, eps);

save_obj_fn_filename = sprintf('obj_fn_mma_%d_%d.csv', input_mode, output_mode);
save_grad_norm_filename = sprintf('grad_norm_mma_%d_%d.csv', input_mode, output_mode);
save_num_obj_calls_filename = sprintf('num_obj_calls_mma_%d_%d.csv', input_mode, output_mode);
save_num_grad_calls_filename = sprintf('num_grad_calls_mma_%d_%d.csv', input_mode, output_mode);

csvwrite(save_obj_fn_filename, obj_fn(1:end-1));
csvwrite(save_grad_norm_filename, grad_norm(1:end-1));
csvwrite(save_num_obj_calls_filename, final_number_calls);
csvwrite(save_num_grad_calls_filename, final_number_calls);

