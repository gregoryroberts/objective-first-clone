function [g_x_partial] = g_x_partial_normalized_overlap(c, x, S_filter)
    S_x = S_filter * x;
    h = (c'*S_x)*conj(c'*S_x);
    dh_dx = c'*S_filter*conj(c'*S_x);

    top = (S_x' * S_x) * dh_dx - h * S_x' * S_filter;
    bottom = (c' * c) * power(S_x' * S_x, 2);
    
    g_x_partial = -top ./ bottom;
end

