clear all
clc
warning('off', 'all')

% Add paths to BenDFO repository
addpath('./m/')
addpath('./data/')
addpath('./m2/')
addpath('./minq5/')

load('dfo.dat');

% Parameters
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", ...
             "relnormal", "reluniform", "relwild", "smooth", "wild3"];
probtypes = probtypes(1:6); % Ensure only first 10 types are used

max_iter = 500;
% C = [1, 0, 0; 0, 1, 0; 0, 0, 1; 1/3, 1/3, 1/3; 1/2, 1/2, 0; 0, 1/2, 1/2; 1/2, 0, 1/2];
% C = [1/2, 0, 1/2]; prob 16 nf=9
 % C = [1/3, 1/3, 1/3]
% C = [0.26, 0, 0.74; 0.2, 0, 0.8; 0.1, 0, 0.9; 0.3, 0, 0.7; 0.4, 0, 0.6; 0.6,0,0.4; 0.7,0,0.3; 0.8, 0,0.2; 0.9,0,0.1]
 
 C=[     0         0    1.0000
         0    0.1000    0.9000
         0    0.2000    0.8000
         0    0.3000    0.7000
         0    0.4000    0.6000
         0    0.5000    0.5000
         0    0.6000    0.4000
         0    0.7000    0.3000
         0    0.8000    0.2000
         0    0.9000    0.1000
         0    1.0000         0
    0.1000         0    0.9000
    0.1000    0.1000    0.8000
    0.1000    0.2000    0.7000
    0.1000    0.3000    0.6000
    0.1000    0.4000    0.5000
    0.1000    0.5000    0.4000
    0.1000    0.6000    0.3000
    0.1000    0.7000    0.2000
    0.1000    0.8000    0.1000
    0.1000    0.9000         0
    0.2000         0    0.8000
    0.2000    0.1000    0.7000
    0.2000    0.2000    0.6000
    0.2000    0.3000    0.5000
    0.2000    0.4000    0.4000
    0.2000    0.5000    0.3000
    0.2000    0.6000    0.2000
    0.2000    0.7000    0.1000
    0.2000    0.8000         0
    0.3000         0    0.7000
    0.3000    0.1000    0.6000
    0.3000    0.2000    0.5000
    0.3000    0.3000    0.4000
    0.3000    0.4000    0.3000
    0.3000    0.5000    0.2000
    0.3000    0.6000    0.1000
    0.3000    0.7000   0.0000
    0.4000         0    0.6000
    0.4000    0.1000    0.5000
    0.4000    0.2000    0.4000
    0.4000    0.3000    0.3000
    0.4000    0.4000    0.2000
    0.4000    0.5000    0.1000
    0.4000    0.6000   0.0000
    0.5000         0    0.5000
    0.5000    0.1000    0.4000
    0.5000    0.2000    0.3000
    0.5000    0.3000    0.2000
    0.5000    0.4000    0.1000
    0.5000    0.5000         0
    0.6000         0    0.4000
    0.6000    0.1000    0.3000
    0.6000    0.2000    0.2000
    0.6000    0.3000    0.1000
    0.6000    0.4000   0.0000
    0.7000         0    0.3000
    0.7000    0.1000    0.2000
    0.7000    0.2000    0.1000
    0.7000    0.3000   0.0000
    0.8000         0    0.2000
    0.8000    0.1000    0.1000
    0.8000    0.2000   0.0000
    0.9000         0    0.1000
    0.9000    0.1000   0.0000
    1.0000         0         0];


% C = [1, 0, 0];

c_len = size(C, 1);

% Initialize storage for results
max_prob = length(probtypes) * size(dfo, 1);
history = zeros(max_iter, max_prob, c_len);
yhist = zeros(max_iter, max_prob, c_len);
iter_hist = zeros(max_iter, max_prob, c_len, 9);

% Disable figure display
set(0, 'DefaultFigureVisible', 'off');  

ip = 0; % Global counter for problem number

% Main Loop for Solving Problems
for p = 1:length(probtypes)
    for row = 1:size(dfo, 1)
        ip = ip + 1; % Increment problem counter

        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp(strcat(int2str(ip), ':'));
        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');

        % Problem setup
        BenDFO.nprob = dfo(row, 1);
        BenDFO.n = dfo(row, 2);
        BenDFO.m = dfo(row, 3);
        BenDFO.factor_power = dfo(row, 4);
        BenDFO.sigma = 1e-2;

        X0 = dfoxs(BenDFO.n, BenDFO.nprob, 10^BenDFO.factor_power);
        func = @(x) calfun(x, BenDFO, probtypes(p));
        dfunc = @(x) out3(x(:)', BenDFO, probtypes(p));
        
        % % Run solver for each coefficient combination
        % for ic = 1:c_len
        %     [hist, ~, ~, ~, ~] = dfo1_4test(max_iter, func, X0, C(ic, :));
        %     history(:, ip, ic) = revise_hist(max_iter, hist, func(X0));
        %     yhist(:, ip, ic) = [reshape(history(:, ip, ic), 1, size(history, 1))];
        % end
        
        % Run solver for each coefficient combination

        n=dfo(row, 2);
        rhobeg=max(1,norm(X0,'inf'));
        
        for ic = 1:c_len
            fprintf('%d, ', ic);
            if mod(ic,15) == 0
                fprintf('\n');
            end
            % [hist, ~, ~, ~, ~] = dfo1_4test(max_iter, func, X0, C(ic, :));
            [X, F, flag, xkin, ihist] = poundersC(func, dfunc, X0', n, (2*n+1), max_iter, 1e-12, C(ic, :), rhobeg, 0, 1, [], 1, -Inf(n,1)', Inf(n,1)', 0,2);
            history(:, ip, ic) = revise_hist(max_iter, F, func(X0));
            yhist(:, ip, ic) = [reshape(history(:, ip, ic), 1, size(history, 1))];
            iter_hist(:, ip, ic, :) = reshape(ihist, size(ihist,1), 1, 1, size(ihist,2));
        end

        fprintf('%d\n', c_len+1);
        [X, F, flag, xkin, ihist] = pounders(func, dfunc, X0', n, 2*n+1, max_iter, 1e-12, [], rhobeg, 0, 1, [], 1, -Inf(n,1)', Inf(n,1)', 0,2);
        history(:, ip, 67) = revise_hist(max_iter, F, func(X0));
        yhist(:, ip, 67) = [reshape(history(:, ip, 67), 1, size(history, 1))];
        iter_hist(:, ip, 67, :) = reshape(ihist, size(ihist,1), 1, 1, size(ihist,2));


        % for ic = 1:c_len
        %     fprintf('%d, ', ic);
        %     if mod(ic,15) == 0
        %         fprintf('\n');
        %     end
        %     % [hist, ~, ~, ~, ~] = dfo1_4test(max_iter, func, X0, C(ic, :));
        %     [X, F, flag, xkin, ihist] = poundersC(func, dfunc, X0', n, (n+3), max_iter, 1e-12, C(ic, :), rhobeg, 0, 1, [], 1, -Inf(n,1)', Inf(n,1)', 0,2);
        %     history(:, ip, ic) = revise_hist(max_iter, F, func(X0));
        %     yhist(:, ip, ic) = [reshape(history(:, ip, ic), 1, size(history, 1))];
        %     iter_hist(:, ip, ic, :) = reshape(ihist, size(ihist,1), 1, 1, size(ihist,2));
        % end
        % 
        % fprintf('%d\n', c_len+1);
        % [X, F, flag, xkin, ihist] = pounders(func, dfunc, X0', n, n+3, max_iter, 1e-12, [], rhobeg, 0, 1, [], 1, -Inf(n,1)', Inf(n,1)', 0,2);
        % history(:, ip, 67) = revise_hist(max_iter, F, func(X0));
        % yhist(:, ip, 67) = [reshape(history(:, ip, 67), 1, size(history, 1))];
        % iter_hist(:, ip, 67, :) = reshape(ihist, size(ihist,1), 1, 1, size(ihist,2));


        
        % for k = 1:
        %     F_best(k) = min(F(1:k)); % 取 F 的前 k 个元素的最小值
        % end
        % yhist(:, ip, 8)=F_best;
    end
end

save history_2n+1
% save history_n+3

% Helper functions
function hist = revise_hist(max_iter, hist, f0)
    % Initialize history if empty
    if isempty(hist)
        hist = f0 * ones(1, max_iter);
    else
        hist(1) = min(f0, hist(1));
        for i = 2:min(max_iter, numel(hist))
            hist(i) = min(hist(i - 1), hist(i));
            % hist(i) = hist(i);
        end
        % Extend or trim hist to match max_iter
        if numel(hist) < max_iter
            hist = [hist; hist(end) * ones(max_iter - numel(hist), 1)];
        else
            hist = hist(1:max_iter);
        end
    end
end

function G=out3(x, BenDFO, ptypes)
	[~, ~, G, ~] = calfun(x, BenDFO, ptypes);
end

