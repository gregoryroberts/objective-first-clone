function [gradient, number_of_calls] = adjoint_gradient(omega, eps, input, g_x_partial, boundary_condition, varargin)
    persistent num_calls;
    
    if (isempty(num_calls))
        num_calls = 0;
    end
    
    if (isempty(varargin))
        num_calls = num_calls + 1;
    end
            
    number_of_calls = num_calls;
    [unused_Ex, unused_Ey, unused_Hz_adj, gradient] = ob1_fdfd_adj(omega, eps, input, g_x_partial, boundary_condition);
end

