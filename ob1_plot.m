function ob1_plot(x, p, dims, option)

N = prod(dims);

Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:2*N), dims);
% Hz = reshape(x(2*N+1:3*N), dims);
% Hz = reshape(x, dims);

e = reshape(p, dims);
% ex = reshape(epsilon(1:N), dims);
% ey = reshape(epsilon(N+1:2*N), dims);

energy = 0.5 * (e .* abs(Ex).^2 + e .* abs(Ey).^2);
switch option
    case 'quick'
        plot_fields(dims, {'\epsilon', e}, {'energy', energy});
        
    case 'full'
        figure(1);
        plot_fields(dims, {'\epsilon', e});

        figure(2);
        plot_fields(dims, ...
            {'Re(Ex)', real(Ex)}, {'Im(Ex)', imag(Ex)}, {'|Ex|', abs(Ex)}, ...
            {'Re(Ey)', real(Ey)}, {'Im(Ey)', imag(Ey)}, {'|Ey|', abs(Ey)});
end

drawnow;
