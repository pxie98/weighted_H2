function H = Hessian(func, x0)% x0为求Hessian矩阵的基准点
    h = 0.01;
    n = size(x0, 1);
    H = zeros(n, n);
    
    for i = 1:n
        e_i = zeros(n,1);
        e_i(i,1) = h;
        for j = i:n
            e_j = zeros(n,1);
            e_j(j,1) = h;
            % 计算\nablda^2 f(xk)
            H(i,j) = (func(x0 + e_i + e_j) - func(x0 + e_i) - func(x0 + e_j) + func(x0)) / (h * h);
            H(j,i) = H(i,j);
        end
    end
end