clear;
clc;
close all;
max_iter = 1000;
x0 = [0; 0]; 
gammas = [0.01, 0.12, 0.05]; % step size
gamma_labels = {'Small step size', 'Too large step size', 'Just right step size'};
colors = ['r', 'g', 'b'];
tol = 1e-6;
%% x1k x2k
for i = 1:length(gammas)
    figure; hold on;
    gamma = gammas(i);
    x = x0;
    x_matrix = x;

    for k = 1:max_iter
        gradient = df0(x(1), x(2));
        y = x - gamma * gradient'; % gradient
        x = project_feasible(y); % project to feasible set
        x_matrix = [x_matrix, x];

        residual = norm(x - project_feasible(x - gamma*df0(x(1),x(2))'));
        residual = norm(x - project_feasible(x - gamma*gradient'));
        if residual < tol
            break;
        end
    end

    k_vals = 0:(size(x_matrix,2)-1);
    plot(k_vals, x_matrix(1,:),'--', 'Color',colors(i),'LineWidth', 2, ...
        'DisplayName', ['x1(k), ' gamma_labels{i}]);
    plot(k_vals, x_matrix(2,:),'-', 'Color', colors(i),'LineWidth', 2, ...
        'DisplayName', ['x2(k), ' gamma_labels{i}]);
    xlabel('Iteration k','FontSize',14,'FontWeight','bold');
    ylabel('Value','FontSize',14,'FontWeight','bold');
    title('Evolution of x1(k) and x2(k)','FontSize', 16,'FontWeight','bold');
    legend('show','FontSize',14,'FontWeight', 'bold');
    grid on;
    set(gca,'FontSize',14,'FontWeight','bold');
    hold off;
end

x1 = linspace(-3,3, 100);
x2 = linspace(-3,3, 100);
[x1,x2] = meshgrid(x1, x2);
f_x1_x2= f0(x1, x2);

for i = 1:length(gammas)
    figure;
    contour(x1, x2, f_x1_x2, 50, 'LineWidth', 2); hold on;
    xlabel('x1', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('x2', 'FontSize', 14, 'FontWeight', 'bold');
    title('Gradient Projection Trajectories', 'FontSize', 16, 'FontWeight', 'bold');
    grid on; set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    tol = 1e-6;

    gamma = gammas(i);
    x = x0;
    trajectory = x'; 
    for k = 1:max_iter
        gradient = df0(x(1), x(2));
        y = x - gamma * gradient';
        x = project_feasible(y);
        trajectory = [trajectory; x'];

        residual = norm(x - project_feasible(x - gamma*df0(x(1),x(2))'));
        if residual < tol
            break;
        end
    end

    plot(trajectory(:,1), trajectory(:,2), 'o-', ...
        'Color', colors(i), 'LineWidth', 2, 'DisplayName', gamma_labels{i});

    %% verify
    x_final = trajectory(end, :)';
    grad = df0(x_final(1), x_final(2));  
    x_star = x_final - gamma * grad';
    x_proj = project_feasible(x_star);
    residual = norm(x - x_proj);

    if residual < tol
        fprintf('(Step size = %.3f) Optimality condition satisfied:residual = %e\n',gamma,residual);
    else
        fprintf('(Step size = %.3f) Optimality condition not satisfied: residual = %e\n',gamma, residual);
    end
    legend('show', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
end

