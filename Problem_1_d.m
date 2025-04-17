clear; clc; close all;

f = @f0; % objective function
grad_f = @df0; % gradient function 

gammas = [0.005,0.12,0.05];
gamma_labels = {'Small step size', 'Too large step size', 'Just right step size'};
colors = ['r', 'g', 'b'];

x1 = linspace(-3,3, 100);
x2 = linspace(-3,3, 100);
[X1,X2] = meshgrid(x1, x2);
f_x1_x2 = f(X1, X2);
tol = 1e-6;

for idx = 1:length(gammas)
    gamma = gammas(idx);
    color = colors(idx);
    
    x = [0;0]; % starting x=(0,0)
    lambda = [0;0]; % starting lambda=(0,0)
    x_matrix = x;
    lambda_matrix = lambda;
    
    max_iter = 10000;
    fprintf('Step size = %.3f\n',gamma);
    for k = 1:max_iter
        % Primal update: minimization of L(x,lambda)
        grad_L = grad_f(x(1),x(2))' - [2;1]*lambda(1) - [1;2]*lambda(2);% gradient of L
        x = x - gamma * grad_L;
        x = max(x, 0); % projection
        
        % Dual update: maximization of g(lambda)
        g1 = 3 - (2*x(1)+x(2)); % constraint
        g2 = 3 - (x(1)+2*x(2)); 
        lambda(1) = lambda(1) + gamma*g1; % gradient ascent
        lambda(2) = lambda(2) + gamma*g2;
        lambda = max(lambda, 0); % project to lambda >= 0

        x_matrix = [x_matrix, x];
        lambda_matrix = [lambda_matrix, lambda];
        [constraint_1,constraint_2,lambda,comp1,comp2,grad_L,optimal_flag] = verifyKKT(x, lambda, tol, grad_f);
        if optimal_flag
            disp(['Converged at iteration ', num2str(k)]);
            break;
        end
    end
    x_final = x_matrix(:, end);
    lambda_final = lambda_matrix(:, end);
    [constraint_1,constraint_2,lambda,comp1,comp2,grad_L,optimal_flag] = verifyKKT(x_final, lambda_final, tol, grad_f);
    
    fprintf('constraint_1 = %e, constraint_2 = %e\n', constraint_1, constraint_2); % primal feasibility
    fprintf('lambda = [%e, %e]\n', lambda(1), lambda(2));% dual feasibility
    fprintf('comp1 = %e, comp2 = %e\n', comp1, comp2);% Complementary Slackness
    fprintf('grad L = %e %e\n', grad_L(1),grad_L(2));% gradient
    if optimal_flag
        fprintf('All KKT conditions are satisfied. The solution is optimal.\n\n\n');
    else
        fprintf('KKT conditions are NOT fully satisfied. The solution may not be optimal.\n\n\n');
    end

    %% Plot 1: Primal variable trajectory on contour plot
    figure;
    contour(X1, X2, f_x1_x2, 50, 'LineWidth', 2); hold on;
    plot(x_matrix(1,:), x_matrix(2,:), 'o-', ...
         'Color', color, 'LineWidth', 2, ...
         'DisplayName', [gamma_labels{idx}]);
    xlabel('x1'); ylabel('x2');
    title(gamma_labels{idx});
    grid on; axis tight;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'best');
    hold off;

    %% Plot 2: Evolution of x1(k) and x2(k)
    figure;
    plot(x_matrix(1,:), '-', ...
         'Color', color, 'LineWidth', 2, ...
         'DisplayName', ['x1(k), ' gamma_labels{idx}]); hold on;
    plot(x_matrix(2,:), '--', ...
         'Color', color, 'LineWidth', 2, ...
         'DisplayName', ['x2(k), ' gamma_labels{idx}]);
    xlabel('Iteration k'); ylabel('Value');
    title(['x1(k), x2(k): ', gamma_labels{idx}]);
    grid on;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'best');
    hold off;

    %% Plot 3: Evolution of lambda1(k) and lambda2(k)
    figure;
    plot(lambda_matrix(1,:), '-', ...
         'Color', color, 'LineWidth', 2, ...
         'DisplayName', ['\lambda_1(k), ' gamma_labels{idx}]); hold on;
    plot(lambda_matrix(2,:), '--', ...
         'Color', color, 'LineWidth', 2, ...
         'DisplayName', ['\lambda_2(k), ' gamma_labels{idx}]);
    xlabel('Iteration k'); ylabel('Value');
    title(['\lambda_1(k), \lambda_2(k): ', gamma_labels{idx}]);
    grid on;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'Location', 'best');
    hold off;
end