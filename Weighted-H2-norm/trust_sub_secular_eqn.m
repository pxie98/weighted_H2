%DFO-H2-norm-a-matlab-version
%Copyright: Pengcheng Xie
%Email: xpc@lsec.cc.ac.cn

function [value] = trust_sub_secular_eqn(lambda_0, eigval, alpha, delta)
    

    m = size(lambda_0,1);
    n = size(eigval,1);
    unn = ones(n, 1);
    unm = ones(m, 1);
    M = eigval * unm' + unn * lambda_0';
    MC = M;
    MM = alpha * unm';
    M(M ~= 0.0) = MM(M ~= 0.0) ./ M(M ~= 0.0);
    M(MC == 0.0) = inf;
    M = M.*M;
    value = sqrt(unm / (M' * unn));

    if size(value(value == inf),1)
        inf_arg = (value == inf);
        value(inf_arg) = zeros(size(value(inf_arg),1), 1);
    end

    value = (1.0./delta) * unm - value;
end

