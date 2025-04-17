% function x_proj = project_feasible(x)
% a1 = [-2; -1]; b1 = -3;
% a2 = [-1; -2]; b2 = -3;
% 
% while (2*x(1) + x(2) < 3) || (x(1) + 2*x(2) < 3)
%     if (2*x(1) + x(2) < 3)
%         x = x + (b1 - a1'*x) / (norm(a1)^2) * a1;
%     end
%     if (x(1) + 2*x(2) < 3)
%         x = x + (b2 - a2'*x) / (norm(a2)^2) * a2;
%     end
%     x(1) = max(x(1), 0);
%     x(2) = max(x(2), 0);
% end
% x_proj = x;
% end



function x_proj = project_feasible(x)
    A = [-2 -1; -1 -2];  
    b = [-3; -3];
    lb = [0; 0];
    options = optimoptions('quadprog', 'Display', 'none');
    x_proj = quadprog(eye(2), -x, A, b, [], [], lb, [], [], options);
end