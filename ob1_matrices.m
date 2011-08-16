function [A, S] = ob1_matrices(dims, active_box)

path(path, '~/wave-tools/helper');

    %
    % Form A matrices.
    %

% Shortcut to form a 2D derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

A{1} = [-(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))]; % Curl for E-field.
A{2} = 0.5 * [S_(0,0)+S_(1,0); S_(0,0)+S_(0,1)]; % Spread for epsilon.
A{3} = [(S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; % Curl for H-field.

    
    %
    % Form selection matrices, used to select active elements of x and p.
    %

S.x = my_selection(dims, [2 2]);
S.p = my_selection(dims, ceil((dims - active_box) / 2));
S.r = my_selection(dims, [1 1]);



function [S] = my_selection(dims, border)
S = zeros(dims);
S(border(1)+1:end-border(1), border(2)+1:end-border(2)) = 1;
ind = find(S);
n = length(ind);
S =  sparse(ind, 1:n, ones(n,1), prod(dims), n);
