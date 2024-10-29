%clear all
%clc
%warning('off','all')

% These 2 lines will need to change depending on your location of the BenDFO repo
addpath('/home/wild/oldcomp/repos/poptus/BenDFO/m/')
addpath('/home/wild/oldcomp/repos/poptus/BenDFO/data/')


load('dfo.dat');

probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"];


max_iter = 500;
M = 1;

ifnoise = 1;
noise = RandStream('mt19937ar', 'Seed', 1);
noise.NormalTransform = 'Ziggurat';

solvers = ["dfo1"];


C = [1,0,0;...
    0,1,0;...
    0,0,1;...
    1/3,1/3,1/3;...
    1/2,1/2,0;...
    0,1/2,1/2];
c_len = size(C,1);


probtypes = probtypes(9); % For now, only run on the smooth version of the problems


max_prob = length(probtypes)*size(dfo, 1);
history = zeros(max_prob, c_len, M, max_iter);
yhist = zeros(max_prob, c_len, M, max_iter+1);
iter_hist = zeros(max_prob, c_len, M, 6, max_iter);



set(0, 'DefaultFigureVisible', 'off');  
t1 = clock;





ip=0; %Global counter of problem number

Results = cell(length(probtypes), size(dfo, 1));
for p = 1:length(probtypes)
    for row = 1:size(dfo, 1) % 53 problems
    ip = ip+1; % Problem counter updated
    
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(strcat(int2str(row), '. ', dfo(row, 1), ':'));
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            BenDFO.nprob = dfo(row, 1);
            BenDFO.n = dfo(row, 2);
            BenDFO.m = dfo(row, 3);
            BenDFO.factor_power = dfo(row, 4);
            BenDFO.sigma = 1e-2; % If you want to add noise

            X0 = dfoxs(BenDFO.n, BenDFO.nprob, 10^BenDFO.factor_power);
            rhobeg =  max(1, norm(X0,'inf'));

   	    func = @(x) calfun(x, BenDFO, probtypes(p));
   	    


        for ic = 1:c_len
            hist = [];
            [hist, x_dfo1, fval_dfo1, ~, ~, ihist] = dfo1_4test(max_iter, func, X0, C(ic,:));
            history(ip, ic, 1, :) = revise_hist(max_iter, hist, func(X0));
            %iter_hist(ip, ic, 1, :, :) = revise_ihist(max_iter, ihist);
        end
    end
end    









function hist = revise_hist(max_iter, hist, f0)
    if size(hist,1) == 0
        hist = f0 * ones(1, max_iter);
    else
        hist(1) = min(f0, hist(1));
        for i = 2:min(max_iter, size(hist))
            hist(i) = min(hist(i-1), hist(i));
        end

        if size(hist,1) >= max_iter
            hist = hist(1:max_iter);
        else
            hist = [hist; hist(end) * ones(max_iter - size(hist,1),1)];
        end
    end
end

function hist = revise_ihist(max_iter, ihist)
    if size(ihist,2) >= max_iter
        hist = ihist(:,1:max_iter);
    else
        hist = [ihist; ihist(:,end) .* ones(6, max_iter - size(ihist,2))];
    end
end
