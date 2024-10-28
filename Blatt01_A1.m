% Define a function
function Blatt01_A1()
    % Call the chatgpt function
    chatgpt();

    % Call the quadpoly function
    % Vectors are defined using square brackets
    %            1
    %[1; 2; 3] = 2 and [1, 2, 3] = (1, 2, 3)
    %            3
    y = quadpoly([1; 2; 3], [1, 2, 3; 1, 2, 3; 1, 2, 3], [1; 2; 3], 1);
    disp(['The result is: ' num2str(y)]);

    % Call the test_quadpoly function
    test_quadpoly1(5);
    test_quadpoly2(5);
end

%Aufgabe 1a
function chatgpt()
    % Generate a random 3 x 3 matrix
    matrix = rand(3);
    % Call the mean function
    meanValue = mean(matrix, 'all');
    % Set a threshold value
    threshold = 0.5;
    % Check if the mean is greater than the threshold
    if meanValue > threshold
        disp('The mean is greater than the threshold.');
    else
        disp('The mean is not greater than the threshold.');
    end

    % Display the matrix and the calculated mean
    disp('Matrix:');
    disp(matrix);
    disp(['Mean: ' num2str(meanValue)]);
end

%Aufgabe 1b
function y = quadpoly(x, A, b, c)
    % Calculate the quadratic polynomial
    y = 0.5 * transpose(x) * A * x + transpose(b) * x + c;
end

%Aufgabe 1c
function test_quadpoly1(cycles)
    intervall = linspace(-2, 2);
    result = zeros(length(intervall), cycles); % Initialize result matrix

    for j = 1:cycles
        % Initialize temp_result
        temp_result = zeros(size(intervall));

        A = (rand(1) - 0.5) * 4;
        b = (rand(1) - 0.5) * 4;
        c = (rand(1) - 0.5) * 4;

        for i = 1:length(intervall)
            x = intervall(i);
            temp_result(i) = quadpoly(x, A, b, c);
        end

        result(:, j) = temp_result; % Store the results of this cycle
    end
    % Create a meshgrid for plotting
    [X, Y] = meshgrid(1:cycles, intervall);
    surf(X, Y, result);

    xlabel('Cycle');
    ylabel('Interval');
    zlabel('Result of quadpoly');
    title('Surface plot of Interval vs. Result of quadpoly1 over Cycles');

end

%Aufgabe 1d
function test_quadpoly2(cycles)
    intervall = linspace(-2, 2);
    result = zeros(length(intervall), cycles); % Initialize result matrix

    for j = 1:cycles
        % Initialize temp_result
        temp_result = zeros(size(intervall));

        A = (rand(2) - 0.5) * 4;
        b = (rand(2, 1) - 0.5) * 4;
        c = (rand(1) - 0.5) * 4;

        for i = 1:length(intervall)
            x = [intervall(i); intervall(i)]; % Ensure x is a vector
            temp_result(i) = quadpoly(x, A, b, c); % Ensure quadpoly returns a scalar
        end

        result(:, j) = temp_result; % Store the results of this cycle
    end

    % Create a meshgrid for plotting
    [X, Y] = meshgrid(1:cycles, intervall);
    surf(X, Y, result);

    xlabel('Cycle');
    ylabel('Interval');
    zlabel('Result of quadpoly');
    title('Surface plot of Interval vs. Result of quadpoly2 over Cycles');
end
