function [constraint_1,constraint_2,lambda,comp1,comp2,grad_L,optimal_flag] = verifyKKT(x, lambda, tol, grad_f)
    constraint_1 = 2*x(1)+x(2)-3; % >=0
    constraint_2 = x(1)+2*x(2)-3; % >=0

    comp1 = lambda(1) * constraint_1;
    comp2 = lambda(2) * constraint_2;

    grad_L = grad_f(x(1), x(2))' - [2; 1]*lambda(1) - [1; 2]*lambda(2);

    % judge
    optimal_flag = true;
    if any(x<-tol) || constraint_1<-tol || constraint_2<-tol
        optimal_flag = false;
    end

    if ~all(lambda >= 0)
        optimal_flag = false;
    end

    if abs(comp1) >= tol || abs(comp2) >= tol
        optimal_flag = false;
    end

    if norm(grad_L) >= tol
        optimal_flag = false;
    end
end