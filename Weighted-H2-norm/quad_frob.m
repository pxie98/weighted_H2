%DFO-H2-norm-a-matlab-version
%Copyright: Pengcheng Xie
%Email: xpc@lsec.cc.ac.cn

function [H, g_hat, c_hat] = quad_frob(X_value, F_values, clist)
    eps = 2.220446049250313e-16;
    tol_svd = eps.^5;
    % n = number of variables m = number of points
    [n, m] = size(X_value);

    H = zeros(n, n);
    g = zeros(n, 1);

    % Shift the points to the origin
    Y = X_value - diag(X_value(1:end, 1)) * ones(n, m);

    if (m < (n+1)*(n+2)/2)
        c1 = clist(1);
        c2 = clist(2);
        c3 = clist(3);
        r = 1;
        % V2 = (np.pi ** (n / 2) / math.gamma(n / 2 + 1))
        
        omega1 = (c1 * r ^ 4 / (2 * (n + 4) * (n + 2)) + c2 * r ^ 2 / (n + 2) + c3);
        omega2 = (c1 * r ^ 2 / (n + 2) + c2);
        omega3 = c1 * r ^ 4 / (4 * (n + 4) * (n + 2));
        omega4 = c1 * r ^ 2 / (n + 2);
        omega5 = c1;
        
        b = [F_values; zeros(n + 1, 1)];
       
        A = 1 / (8 * omega1) * (Y' * Y).^2;
        len_A = size(A,1);
        for i = 1:len_A
            for j = 1:len_A
                A(i,j) = A(i,j) - (omega3 / (8 * omega1 * (n * omega3 + omega1))) * (Y(1:end, i)' * Y(1:end, i)) * (Y(1:end, j)' * Y(1:end, j));
            end
        end

        J = zeros(m, 1);
        
        for i = 1:m
            J(i) = 1 - omega4 / (4 * omega1 + 4 * n * omega3) * (Y(1:end, i)' * Y(1:end, i));
        end
        

        line1 = [A, J, Y'];  
        line2 = [J', ((n * omega4 ^ 2) / (2 * n * omega3 + 2 * omega1) - 2 * omega5) * ones(1, 1), zeros(1,n)];
        line3 = [Y, zeros(n, 1), -2 * omega2 * eye(n)];

        W = [line1; line2; line3];
        lambda_0 = quad_Frob_compute_coeffs(W, tol_svd, b, 'partial');

        % Grab the coeffs of linear terms (g) and the ones of quadratic terms
        % (H) for g.T s + s.T H s
        g = lambda_0(m+2:end);
        
        c = lambda_0(m+1);    % new c

        inner_sum = 0;
        for j = 1:m
            inner_sum = inner_sum + lambda_0(j) * (reshape(Y(1:end, j), 1, n) * reshape(Y(1:end, j), n, 1));
        end

        H = H - (1 / (2 * omega1)) * (2 * omega3 * ((1 / (2 * (2 * n * omega3 + 2 * omega1))) * inner_sum - n * omega4 * c / (2 * n * omega3 + 2 * omega1)) + omega4 * c) * eye(n);
        for j = 1:m
            H = H + 1 / (4 * omega1) * (lambda_0(j) * (reshape(Y(1:end, j), n, 1) * reshape(Y(1:end, j), 1, n)));
        end
        
        g_hat = g - H * X_value(1:end, 1);
        c_hat = c + 0.5 * X_value(1:end, 1)' * H * X_value(1:end, 1) - g' * X_value(1:end, 1);
        
    else  % Construct a full model
        % Here enough points. Solve the sys of equations.
        b = F_values;
        phi_Q = [];
        for i = 1:m
            y = Y(1:end, i);
            y = y(newaxis);    
            yguodu = (y.^2);
            aux_H = y * y' - 0.5 * diag(yguodu(0));
            aux = [];
            for j = 1:n
                aux = [aux, aux_H(j:n, j)];
            end

            phi_Q = [phi_Q; aux];
        end

        W = [ones(m, 1), Y.T];
        W = [W, phi_Q];

        lambda_0 = quad_Frob_compute_coeffs(W, tol_svd, b, option='full');

        % Retrieve the model coeffs (g) and (H)
        g = lambda_0(1:n+1, 1:end);  % 1 n
        cont = n+1;
        H = zeros(n, n);

        for j = 1:n   
            H(j:n, j) = lambda_0(cont:cont + n - j, 1:end);  % n+1+n- 0  n-1 2n+1 n+2
            cont = cont + n - j;
        end

        H = H + H' - diag(diag(H));
        
        
        c = [];
        g_hat = [];
        c_hat = [];
    end
end

