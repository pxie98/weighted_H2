
function [H, g] = Cquad_frob(X_value, F_values, delta, X_center, clist) %, clist

    X_value=X_value';
X_center=X_center';

    eps = 2.220446049250313e-16;
    tol_svd = eps.^5;

    




    [n, m] = size(X_value);

    H = zeros(n, n);
    g = zeros(n, 1);

   
    

     Y = X_value - diag(X_center) * ones(n, m);

      Y = (X_value - diag(X_center) * ones(n, m))/delta;

    if (m < (n+1)*(n+2)/2)
        c1 = clist(1);
        c2 = clist(2);
        c3 = clist(3);
       



        r=1;

        
        omega1 = (c1 * r ^ 4 / (2 * (n + 4) * (n + 2)) + c2 * r ^ 2 / (n + 2) + c3);
        omega2 = (c1 * r ^ 2 / (n + 2) + c2);
        omega3 = c1 * r ^ 4 / (4 * (n + 4) * (n + 2));
        omega4 = c1 * r ^ 2 / (n + 2);
        omega5 = c1;
   
        b = [F_values; zeros(n + 1, 1)];

        A = (1 / (8 * omega1)) * (Y' * Y).^2;
        Y_sums = sum(Y.^2, 1);
        Bxpc = Y_sums' * Y_sums;
        A = A - (omega3 / (8 * omega1 * (n * omega3 + omega1))) * Bxpc;

        J = zeros(m, 1);
        for i = 1:m
            J(i) = 1 - omega4 / (4 * omega1 + 4 * n * omega3) * (Y(1:end, i)' * Y(1:end, i));
        end
        

        line1 = [A, J, Y'];  
        line2 = [J', ((n * omega4 ^ 2) / (2 * n * omega3 + 2 * omega1) - 2 * omega5) * ones(1, 1), zeros(1,n)];
        line3 = [Y, zeros(n, 1), -2 * omega2 * eye(n)];

        W = [line1; line2; line3];

        lambda_0 = pinv(W) * b;

        g = lambda_0(m+2:end);
        c = lambda_0(m+1);    

        inner_sum = 0;
        for j = 1:m
            inner_sum = inner_sum + lambda_0(j) * (reshape(Y(1:end, j), 1, n) * reshape(Y(1:end, j), n, 1));
        end

        H = H - (1 / (2 * omega1)) * (2 * omega3 * ((1 / (2 * (2 * n * omega3 + 2 * omega1))) * inner_sum - n * omega4 * c / (2 * n * omega3 + 2 * omega1)) + omega4 * c) * eye(n);

        Lambdaxpc =zeros(m,m);
        Lambdaxpc = diag(1 / (4 * omega1) * (lambda_0(1:m)));

        H = H + Y * Lambdaxpc * Y';
        
        H = (H + H') / 2;


   


H = H / delta^2;
g = g / delta;

end

