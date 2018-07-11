function [row_dims, col_dims] = eps_dims_to_change(eps)

    start_row_idx = 4;
    start_col_idx = 3;

    [eps_rows, eps_cols] = size(eps);

    row_dims = start_row_idx:(eps_rows - 2);
    col_dims = start_col_idx:(eps_cols - 2);

end

