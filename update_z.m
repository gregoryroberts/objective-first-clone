function [new_z] = update_z(cur_z, next_z, min_z, max_z)
    new_z = cur_z;
    
    [row_dims, col_dims] = eps_dims_to_change(cur_z);
    
    new_z(row_dims,col_dims) = next_z(row_dims,col_dims);
    new_z(row_dims,col_dims) = min(new_z(row_dims,col_dims), max_z(row_dims,col_dims));
    new_z(row_dims,col_dims) = max(new_z(row_dims,col_dims), min_z(row_dims,col_dims));
end

