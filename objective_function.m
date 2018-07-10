function [value, number_of_calls] = objective_function(c, omega, eps, input, boundary_condition, varargin)
    persistent num_calls;
    
    if (isempty(num_calls))
        num_calls = 0;
    end
    
    if (isempty(varargin))
        num_calls = num_calls + 1;
    end

    number_of_calls = num_calls;

    [unused_Ex, unused_Ey, Hz] = ob1_fdfd(omega, eps, input, boundary_condition);
    extract_Hz = Hz(end,:).';
%     value = (-(c'*extract_Hz)*conj(c'*extract_Hz)) / ((extract_Hz'*extract_Hz) * (c' * c));
    value = (-(c'*extract_Hz)*conj(c'*extract_Hz));
end

