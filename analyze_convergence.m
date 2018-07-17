close all; clf; clc;

input_modes = [1, 1];
output_modes = [2, 3];

methods = {'backtrack', 'fixed', 'mma', 'conj_grad_pr', 'conj_grad_fr'};
title_methods = {'Backtrack', 'Fixed', 'MMA', 'Conj. Grad. (Polak, Ribiere)', 'Conj. Grad. (Fletcher, Reeves)'};
colors = {'g', 'b', 'r', 'm', 'c'};
plot_types = {'obj_fn', 'grad_norm'};
title_plot_types = {'Objective Function', 'Gradient Norm'};

for m = 1 : 1 : length(plot_types)
    for k = 1 : 1 : length(input_modes)
    
        figure;
        hold on;
        
        make_title = sprintf('%s (%d -> %d)', title_plot_types{m}, input_modes(k), output_modes(k));
        title(make_title);
        
        for n = 1 : 1 : length(methods)
            load_filename = sprintf('%s_%s_%d_%d.csv', ...
                plot_types{m}, methods{n}, input_modes(k), output_modes(k));
            load_obj_fn_calls = sprintf('num_obj_calls_%s_%d_%d.csv', ...
                methods{n}, input_modes(k), output_modes(k));
            load_grad_calls = sprintf('num_grad_calls_%s_%d_%d.csv', ...
                methods{n}, input_modes(k), output_modes(k));
            get_color = colors{n};
            data = csvread(load_filename);
            
            get_obj_fn_calls = csvread(load_obj_fn_calls);
            get_grad_calls = csvread(load_grad_calls);
            
            legends{n} = sprintf('%s: %d Obj Fn Calls, %d Grad Calls', ...
                title_methods{n}, get_obj_fn_calls, get_grad_calls);
            
            plot(1:length(data), data, get_color);
        end
        
        legend(legends);

        hold off;
    end
end
