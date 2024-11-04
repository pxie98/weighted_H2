function df = nabla_f(func, x0)% x0为求Hessian矩阵的基准点
    h = 0.1;
    n = size(x0, 1);
    df = zeros(n, 1);
    
    for i = 1:n
        e_i = zeros(n,1);
        e_i(i,1) = h;
        % 计算\nablda f(xk)
        df(i,1) = (func(x0 + e_i) - func(x0 - e_i)) / (2*h);
        % disp(func(x0 + e_i))
        % disp(func(x0 - e_i))
    end
end
