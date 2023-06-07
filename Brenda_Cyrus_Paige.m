%{
Group Members:
1. Brenda Cheptoo - ENE211-0004/2018
2. Cyrus Muthui - ENE211-0010/2018
3. Paige Muva - ENE211-0017/2018
%}
% Coefficients of the equations
A = [5 1 2; 1 4 -2; 2 3 8];
% Constants on the right-hand side
B = [19; -2; 39];
% Initial solution guess
X = [1; 1; 1];
% Maximum number of iterations
maxIterations = 100;
% Error tolerance
tolerance = 1e-6;

% Size of the system
n = length(B);

% Matrix to store values of x, y, z at each iteration
valuesMatrix = zeros(maxIterations, n);

% Gauss-Seidel iteration
for k = 1:maxIterations
    X_old = X; % Store the previous iteration solution
    
    for i = 1:n
        sigma = 0;
        
        for j = 1:n
            if j ~= i
                sigma = sigma + A(i, j) * X(j);
            end
        end
        
        X(i) = (B(i) - sigma) / A(i, i);
    end
    
    % Store values of x, y, z at each iteration
    valuesMatrix(k, :) = X';
    
    % Check convergence
    if norm(X - X_old, inf) < tolerance
        disp('Solution converged:');
        disp(X);
        break;
    end
end

% Trim the valuesMatrix to remove unused rows
valuesMatrix = valuesMatrix(1:k, :);

% Display the matrix of values after each iteration
disp('Values of x, y, z after each iteration:');
disp(valuesMatrix);

% Plot the values of x, y, z over the number of iterations
iterations = 1:k;
plot(iterations, valuesMatrix(:, 1), 'r-', iterations, valuesMatrix(:, 2), 'g-', iterations, valuesMatrix(:, 3), 'b-');
xlabel('Number of Iterations');
ylabel('Values of x, y, z');
legend('x', 'y', 'z');
