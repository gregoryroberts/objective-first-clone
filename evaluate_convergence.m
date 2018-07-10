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
% the overlap of the output mode from the structure with the desired output
% waveguide mode specified.  As a minimization problem, this is equivalent
% ot minimizing the negative of that same overlap integral.
% compute_objective = @(omega, eps, input, boundary_condition) objective_function(c, omega, eps, input, boundary_condition);
% Partial derivative of the objective function, g, with respect to the
% field variable input_Hz
g_x_partial = @(input_Hz, S_filter) -c'*S_filter*conj(c'*(S_filter*input_Hz));
% g_x_partial = @(input_Hz, S_filter) g_x_partial_normalized_overlap(c, input_Hz, S_filter);

%% Evaluate Fixed Step Size
clear objective_function;
clear adjoint_gradient;

[fixed_step_performance_plot, fixed_step_obj_calls, fixed_step_adj_calls, fixed_final_eps] = ...
    evaluate_descent_method(@gradient_descent_fixed_step_size, ...
    c, spec.omega, spec.in, spec.bc, eps, ...
    min_eps_array, max_eps_array, g_x_partial, @objective_function, ...
    @adjoint_gradient);

fprintf('Fixed step size\nFinal Performance: %f\nNumber objective calls: %d\nNumber adjoint calls: %d\n\n', ...
    fixed_step_performance_plot(end), fixed_step_obj_calls, fixed_step_adj_calls);

% save_epsilon_filename = sprintf('adjoint_eps_bench_%d_%d.csv', input_mode, output_mode);
% csvwrite(save_epsilon_filename, fixed_final_eps);

figure;
hold on;
plot(fixed_step_performance_plot, 'g');

% simulate(spec, fixed_final_eps, [simulation_width simulation_height]);

%% Evaluate Backtracking Line Search
clear objective_function;
clear adjoint_gradient;

[backtrack_performance_plot, backtrack_step_obj_calls, backtrack_step_adj_calls, backtrack_final_eps] = ...
    evaluate_descent_method(@gradient_descent_backtracking_line_search, ...
    c, spec.omega, spec.in, spec.bc, eps, ...
    min_eps_array, max_eps_array, g_x_partial, @objective_function, ...
    @adjoint_gradient);

fprintf('Backtracking Line Search\nFinal Performance: %f\nNumber objective calls: %d\nNumber adjoint calls: %d\n\n', ...
    backtrack_performance_plot(end), backtrack_step_obj_calls, backtrack_step_adj_calls);

% save_epsilon_filename = sprintf('adjoint_eps_bench_%d_%d.csv', input_mode, output_mode);
% csvwrite(save_epsilon_filename, backtrack_final_eps);

plot(backtrack_performance_plot, 'b');

% simulate(spec, backtrack_final_eps, [simulation_width simulation_height]);
