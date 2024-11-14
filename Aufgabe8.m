close all;
clear vars;

% l >= 7
l = 12;

dim_inp = 6:l;
dim_out = 3;

% time taken for each dimension 2^6 - 2^l
time_Linvb = zeros(l - 5, 1);
time_Linvb2 = zeros(l - 5, 1);
time_Linvb_Vec = zeros(l - 5, 1);
time_Linvb2_Vec = zeros(l - 5, 1);

% time per function for each dimension
for t = dim_inp
    dim = 2 ^ t;
    L = non_sing_tril_gen(dim);
    Z = rand(dim, dim_out);

    tic;
    sol_Linvb = Linvb(L, Z);
    time_Linvb(t - 5) = toc;

    tic;
    sol_Linvb2 = Linvb2(L, Z);
    time_Linvb2(t - 5) = toc;

    tic;
    sol_Linvb_Vec = Linvb_Vec(L, Z);
    time_Linvb_Vec(t - 5) = toc;

    tic;
    sol_Linvb2_Vec = Linvb2_Vec(L, Z);
    time_Linvb2_Vec(t - 5) = toc;
end

loglog(dim_inp, time_Linvb, 'DisplayName', 'Linvb'); hold on;
loglog(dim_inp, time_Linvb2, 'DisplayName', 'Linvb2'); hold on;
loglog(dim_inp, time_Linvb_Vec, 'DisplayName', 'Linvb\_Vec'); hold on;
loglog(dim_inp, time_Linvb2_Vec, 'DisplayName', 'Linvb2\_Vec');
xlabel('Dimension der Matrix 2^t');
ylabel('Zeit in Sekunden');
legend show;
grid on;

% generate a regular lower triangular matrix with dimension n
function A = non_sing_tril_gen(n)
    A = tril(randn(n));

    % regenerate if matrix is singular
    if det(A) == 0
        A = non_sing_tril_gen(n);
    end

end

% compute L*x=Z with lower triangular L for each column in Z
function X = Linvb(L, Z)
    n = size(L, 1);
    m = size(Z, 2);
    X = zeros(n, m);

    for i = 1:m
        z = Z(:, i);
        x = zeros(n, 1);

        for j = 1:n

            for k = 1:j - 1
                z(j) = z(j) - L(j, k) * x(k);
            end

            x(j) = z(j) / L(j, j);
        end

        X(:, i) = x;
    end

end

% compute L*x=z with lower triangular L for each column in Z
function X = Linvb2(L, Z)
    [n, m] = size(Z);
    X = zeros(n, m);

    for i = 1:m
        z = Z(:, i);
        x = zeros(n, 1);

        for j = 1:n
            x(j) = z(j) / L(j, j);

            for k = j + 1:n
                z(k) = z(k) - L(k, j) * x(j);
            end

        end

        X(:, i) = x;
    end

end

% compute L*x=Z with lower triangular L for each column in Z (vectorized)
function X = Linvb_Vec(L, Z)
    [n, m] = size(Z);
    X = zeros(n, m);

    for i = 1:m
        z = Z(:, i);
        x = zeros(n, 1);

        for j = 1:n
            z(j) = z(j) - L(j, 1:j - 1) * x(1:j - 1);
            x(j) = z(j) / L(j, j);
        end

        X(:, i) = x;
    end

end

% compute L*x=z with lower triangular L for each column in Z (vectorized)
function X = Linvb2_Vec(L, Z)
    n = size(L, 1);
    m = size(Z, 2);
    X = zeros(n, m);

    for i = 1:m
        z = Z(:, i);
        x = zeros(n, 1);

        for j = 1:n
            x(j) = z(j) / L(j, j);
            z(j + 1:n) = z(j + 1:n) - L(j + 1:n, j) * x(j);
        end

        X(:, i) = x;
    end

end
