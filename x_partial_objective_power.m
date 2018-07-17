function [partial] = x_partial_objective_power(spec_Hz, spec_Ey, Hz, Ey_hat, S_filter)
    projected_Hz = spec_Hz * (spec_Hz' * S_filter * Hz);
    projected_Ey = spec_Ey * (spec_Ey' * S_filter * Ey_hat * Hz);
    
    first_part_derivative = projected_Hz' * spec_Ey * (spec_Ey' * S_filter * Ey_hat);
    second_part_derivative = (spec_Hz * (spec_Hz' * S_filter)).' * conj(projected_Ey);
    
    partial = -0.5 * (first_part_derivative.' + second_part_derivative);
end

