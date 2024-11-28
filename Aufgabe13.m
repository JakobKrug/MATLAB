close all;
clear vars;

%% 13d)
n = 10;
diagonals = {rand(n-1,1), rand(n,1), rand(n-1,1)};
b = rand(n,1);

test_lr_tridiag(diagonals{1}, diagonals{2}, diagonals{3}, b);
    
diagonals{2}(3) = 0;
test_lr_tridiag(diagonals{1}, diagonals{2}, diagonals{3}, b);


function [lund, rhd, rond] = my_lr_tridiag(und, hd, ond)
    n = length(hd);        
    lund = zeros(n-1, 1); rhd = hd; rond = ond; 

    if any(rhd == 0)
        error('Error: 0 on main diag encountered, cannot perform LU-decomposition w/o pivoting');
    end

    for i = 1:n-1
        lund(i) = und(i) / rhd(i);
        rhd(i+1) = rhd(i+1) - lund(i) * rond(i);
    end
end

function test_lr_tridiag(und, hd, ond, b)
    try
        [lund, rhd, rond] = my_lr_tridiag(und, hd, ond);
        disp('Decomposition of:'); disp(diag(und,-1)+diag(hd)+diag(ond,1)); disp(diag(lund,-1)+diag(rhd)+diag(rond,1));
        
        x = solve_lr_tridiag(lund, rhd, rond, b);
        A = diag(lund,-1)+diag(rhd)+diag(rond,1);
        disp(x); disp((diag(lund,-1)+diag(rhd)+diag(rond,1))\b);
    catch e
        disp(e.message);
    end
end

function x = solve_lr_tridiag(lund, rhd, rond, b)
    n = size(rhd, 1);
    
    % L*y=b
    y = zeros(n,1); y(1) = b(1);
    for i = 2:n
         y(i) = b(i) - lund(i-1) * y(i-1);
    end

    % R*x=y
    x = zeros(n,1); x(n) = y(n) / rhd(n);
    for i = n-1:-1:1
        x(i) = (y(i)-x(i+1)*rond(i)) / rhd(i);
    end
end