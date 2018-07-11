function [new_eps] = update_eps(cur_eps, next_eps, min_eps, max_eps)
    new_eps = cur_eps;
    
    [row_dims, col_dims] = eps_dims_to_change(cur_eps);
    
    new_eps(row_dims,col_dims) = next_eps(row_dims,col_dims);
    new_eps(row_dims,col_dims) = min(new_eps(row_dims,col_dims), max_eps(row_dims,col_dims));
    new_eps(row_dims,col_dims) = max(new_eps(row_dims,col_dims), min_eps(row_dims,col_dims));
end

