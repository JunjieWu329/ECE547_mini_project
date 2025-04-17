%% 3D plot
x1 = linspace(-5,5,100);
x2 = linspace(-5,5,100);
[x1,x2] = meshgrid(x1, x2);
f_x1_x2 = f0(x1, x2);
figure;
mesh(x1, x2, f_x1_x2);
xlabel('x1');
ylabel('x2');
zlabel('f(x1,x2)');
title('3D Plot of f(x1,x2)');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;

%% 2D plot
rng(0);
directions = -1+2*rand(3,2); % [-1,1]
t = linspace(-5, 5, 100);
figure;
hold on;
x1_start = 0;
x2_start = 0;
for i = 1:3
    d = directions(i,:); 
    x1_line = x1_start+d(1)*t;
    x2_line = x2_start+d(2)*t;

    f_x1_x2 = f0(x1_line, x2_line);
    plot(t, f_x1_x2, 'DisplayName', sprintf('Direction [%.2f, %.2f]', d(1), d(2)),'LineWidth',2);
end
xlabel('t');
ylabel('f(t)');
title('Function Restricted to A Line');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
legend show;
grid on;
hold off;