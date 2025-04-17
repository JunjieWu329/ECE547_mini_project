%% x1(k) and x2(k) 
clear;
clc;
close all;

max_iter = 1000;
x0 = [0;0]; 
gammas = [0.01, 0.12, 0.05]; % step size
gamma_labels = {'Small step size','Too large step size','Just right step size'};
colors = ['r','g','b'];
tol = 1e-4;

for i = 1:length(gammas)
    figure;  
    hold on;
    gamma = gammas(i);
    x = x0;
    x_matrix = x;
    
    for k = 1:max_iter
        gradient = df0(x(1),x(2));
        x = x - gamma * gradient';
        x_matrix = [x_matrix,x];
        
        new_gradient = df0(x(1), x(2));
        if norm(new_gradient) < tol
            break;
        end
    end
    
    k_vals = 0:(size(x_matrix,2)-1);
    plot(k_vals, x_matrix(1,:), '--', 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', ['x1(k), ' gamma_labels{i}]);
    plot(k_vals, x_matrix(2,:), '-', 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', ['x2(k), ' gamma_labels{i}]);
    
    xlabel('Iteration k', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Value', 'FontSize', 14, 'FontWeight', 'bold');
    legend('show', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
    xlim([0 max_iter])
end

%% contour plot 
x1_vals = linspace(-3,3, 100);
x2_vals = linspace(-3,3, 100);
[X1, X2] = meshgrid(x1_vals, x2_vals);
f_x1_x2 = f0(X1, X2);

for i = 1:length(gammas)
    figure;
    contour(X1, X2, f_x1_x2, 50, 'LineWidth', 2);
    hold on;
    xlabel('x1', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('x2', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    
    gamma = gammas(i);
    x = x0;
    trajectory = x';
    for k = 1:max_iter
        gradient = df0(x(1), x(2));
        x = x - gamma * gradient';
        trajectory = [trajectory; x'];
        
        if norm(gradient') < 1e-6
            break;
        end
    end
    
    plot(trajectory(:,1), trajectory(:,2), 'o-', ...
        'Color', colors(i), 'LineWidth', 2, 'DisplayName', gamma_labels{i});
    
    % verify
    x_final = trajectory(end, :)';
    grad_f = df0(x_final(1), x_final(2))';
    if norm(grad_f) < tol
        fprintf('(Step size = %.3f) Optimality condition satisfied: gradient = [%e %e]\n', ...
            gamma, grad_f(1), grad_f(2));
    else
        fprintf('(Step size = %.3f) Optimality condition not satisfied: gradient = [%e %e]\n', ...
            gamma, grad_f(1), grad_f(2));
    end
    legend('show', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
end