close all;
clear vars;

n = 10;
b = zeros(n, 1);
A = hilb(n);
x = ones(n, 1);

for j = 1:n
    b(j) = sum(1 ./ ((1:n) + j - 1));
end

e = ones(n, 1);
B = spdiags([-e 2 * e -e], -1:1, n, n);
B = full(B);

%% keine Störung
[condition, abs_err, residuum] = calc_parameters(b, A, x);
disp(['Hilbert Matrix' newline '(w/o error) Kondition: ' num2str(condition) ' Absoluter Fehler: ' num2str(abs_err) ' Residuum: ' num2str(residuum)]);

[b_condition, b_abs_err, b_residuum] = calc_parameters(b, B, x);
disp(['Saiten Matrix' newline '(w/o error) Kondition: ' num2str(b_condition) ' Absoluter Fehler: ' num2str(b_abs_err) ' Residuum: ' num2str(b_residuum)]);

%% Störung
eps = 1e-3;
b_s = b .* (1 + eps * rand(n, 1));
err_amount = norm(b - b_s, 2);
[condition_s, abs_err_s, residuum_s] = calc_parameters(b_s, A, x);
disp([newline 'Hilbert Matrix' newline '(w error) Eingabefehler: ' num2str(err_amount) ' Absoluter Fehler: ' num2str(abs_err_s) ' Residuum: ' num2str(residuum_s)]);

[b_condition_s, b_abs_err_s, b_residuum_s] = calc_parameters(b_s, B, x);
disp(['Saiten Matrix' newline '(w error) Eingabefehler: ' num2str(err_amount) ' Absoluter Fehler: ' num2str(b_abs_err_s) ' Residuum: ' num2str(b_residuum_s)]);

%% Abschätzung relativer Fehler
rel_err_b = rel_error_bound(b, A, residuum);
rel_err_b_s = rel_error_bound(b, A, residuum_s);
disp([newline 'Hilbert Matrix' newline 'Abschätzung relativer Fehler: (w/o error) ' num2str(rel_err_b) ' (w error) ' num2str(rel_err_b_s)]);

b_rel_err_b = rel_error_bound(b, A, b_residuum);
b_rel_err_b_s = rel_error_bound(b, B, b_residuum_s);
disp(['Saiten Matrix' newline 'Abschätzung relativer Fehler: (w/o error) ' num2str(b_rel_err_b) ' (w error) ' num2str(b_rel_err_b_s)]);

%% Abschätzung Fehler
err_b = err_bound(condition, b, b_s);
disp([newline 'Hilbert Matrix' newline 'Abschätzung Fehler: ' num2str(err_b)]);

b_err_b = err_bound(b_condition, b, b_s);
disp(['Saiten Matrix' newline 'Abschätzung Fehler: ' num2str(b_err_b)]);

%%
function [c, a, r] = calc_parameters(b, A, x)
    xn = A \ b;
    c = cond(A, 2);
    a = norm(x - xn, 2);
    r = norm(A * xn - b, 2);
end

function rb = rel_error_bound(b, A, r)
    rb = norm(inv(A), 2) * norm(A, 2) * r / norm(b, 2);
end

% k(A) * (norm(err_b)/norm(b)
function ef = err_bound(cond, b, b_s)
    ef = cond * (norm(b_s) / norm(b));
end
