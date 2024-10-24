% Define a function
function Blatt01()
    % Generate a random 3 x 3 matrix
    matrix = rand(3, 3);
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

    % Call the quadpoly function
    % Vectors are defined using square brackets
    %            1
    %[1; 2; 3] = 2 and [1, 2, 3] = (1, 2, 3)
    %            3 
    quadpoly([1; 2; 3], [1, 2, 3; 1, 2, 3; 1, 2, 3], [1; 2; 3], 1);

    test_quadpoly();
end

function quadpoly(x, A, b, c)
    % Calculate the quadratic polynomial
    y = 0.5 * transpose(x) * A * x + transpose(b) * x + c;
    % Display the result
    disp(['The result is: ' num2str(y)]);
end

function test_quadpoly()
    intervall = linspace(-2, 2);
    for i = 1:length(intervall)
        x = intervall(i);
        quadpoly(x, 1, 1, 1);
    end
end