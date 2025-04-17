clc;clear all;close all
Num_Links=12;
Num_Flows=7;
Max_Links_On_Path=4;

Link_Capacity=[1.0 1.0 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 2.0 2.0];
Flow_Weight=[1.0 1.0 1.0 1.0 2.0 2.0 2.0];

Flow_Path=[ [3 9 -1 -1]; ...
    [4 9 -1 -1]; ...
    [3 10 -1 -1]; ...
    [4 10 -1 -1]; ...
    [5 8 -1 -1]; ...
    [5 9 -1 -1]; ...
    [1 6 11 -1]];

A = zeros(Num_Links, Num_Flows); % 12*7
for i = 1:Num_Flows % ith flow
    for k = 1:Max_Links_On_Path % 1:1:4
        j = Flow_Path(i, k);
        if j ~= -1
            A(j, i) = 1; % flow i passes through link j
        end
    end
end

lambda = ones(Num_Links, 1)*0.1; % initialize for Lagrange multipliers
gamma = 0.25;% step size
Max_Iter = 200; 

lambda_matrix = zeros(Num_Links, Max_Iter);
x_matrix = zeros(Num_Flows, Max_Iter); 

for iter = 1:Max_Iter
    x = zeros(Num_Flows, 1); % rate flow matrix x
    for i = 1:Num_Flows
        idx = find(A(:,i) == 1); % find j:Aji=1
        if isempty(idx)
            continue;
        end
        denominator = sum(lambda(idx));
        if denominator < 1e-10
            x(i) = 0; 
        else
            x(i) = Flow_Weight(i) / denominator;
        end
    end

    gradient_g = Link_Capacity' - A*x;
    lambda = lambda-gamma*gradient_g; % update lambda
    lambda = max(lambda,0); % lambda >= 0

    lambda_matrix(:, iter) = lambda;
    x_matrix(:, iter) = x;
end

figure;
plot(x_matrix','LineWidth',2);
xlabel('Iteration');
ylabel('Flow Rates');
title('Convergence of Flow Rates');
legend('Flow 1', 'Flow 2', 'Flow 3', 'Flow 4', 'Flow 5', 'Flow 6', 'Flow 7');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure;
colors = jet(12);  hold on;
for i = 1:12
    plot(lambda_matrix(i, :), 'LineWidth', 2, 'Color', colors(i, :));
end
hold off;
legend('Link 1','Link 2','Link 3','Link 4','Link 5','Link 6', ...
       'Link 7','Link 8','Link 9','Link 10','Link 11','Link 12'); 
xlabel('Iteration');
ylabel('Dual Variables (\lambda)');
title('Convergence of Dual Variables');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'bold');


%% Verify final results from (b)
final_flow_rate = x_matrix(:, end);
final_lambda    = lambda_matrix(:, end);
tol = 1e-6;
% 1) Primal: Ax <= c
Ax = A * final_flow_rate;
max_Ax_minus_c = max(Ax - Link_Capacity');
fprintf('=== Primal feasibility: Ax <= c ===\n');
fprintf('Ax values:'); fprintf('%g ', Ax); fprintf('\n');
fprintf('Link capacities (c):'); fprintf('%g ', Link_Capacity); fprintf('\n');
fprintf('Max(Ax-c) = %e\n', max_Ax_minus_c);
if max_Ax_minus_c <= tol
    fprintf('Result: PASS\n');
else
    fprintf('Result: FAIL\n');
end

% 2) Primal: x >= 0
min_flow = min(final_flow_rate);
fprintf('\n=== Primal feasibility:x >= 0 ===\n');
fprintf('Final flow rates:'); fprintf('%g ', final_flow_rate); fprintf('\n');
fprintf('Min(x) = %e\n', min_flow);
if min_flow >= -tol
    fprintf('Result: PASS\n');
else
    fprintf('Result: FAIL\n');
end

% 3) Dual: lambda >= 0
min_lambda = min(final_lambda);
fprintf('\n=== Dual feasibility: lambda >= 0 ===\n');
fprintf('Final lambda values:'); fprintf('%g ', final_lambda); fprintf('\n');
fprintf('Min(lambda) = %e\n', min_lambda);
if min_lambda >= -tol
    fprintf('Result: PASS\n');
else
    fprintf('Result: FAIL\n');
end

% 4) Complementary slackness: lambda_i * (Ax_i - c_i) = 0
comp_slack = final_lambda .* (Ax - Link_Capacity');
max_comp_violation = max(abs(comp_slack));
fprintf('\n=== Complementary slackness ===\n');
fprintf('lambda.* (Ax-c):'); fprintf('%g ', comp_slack); fprintf('\n');
fprintf('Max|lambda.*(Ax-c)| = %e\n', max_comp_violation);
if max_comp_violation <= tol
    fprintf('Result: PASS\n');
else
    fprintf('Result: FAIL\n');
end

% 5) Stationarity: w_i/x_i - sum_j lambda_j = 0
stationarity = zeros(Num_Flows,1);
for i = 1:Num_Flows
    idx = find(A(:,i));
    stationarity(i) = Flow_Weight(i)/final_flow_rate(i) - sum(final_lambda(idx));
end
max_stat_violation = max(abs(stationarity));
fprintf('\n=== Stationarity ===\n');
fprintf('w_i/x_i-sum lambda:'); fprintf('%g ', stationarity); fprintf('\n');
fprintf('Max|stationarity| = %e\n', max_stat_violation);
if max_stat_violation <= tol
    fprintf('Result: PASS\n');
else
    fprintf('Result: FAIL\n');
end




