clear all
clc
warning('off','all')

filename = 'problems1';
fileID = fopen(filename);
C = textscan(fileID, '%s %f');
prob = C{1};
tstn = C{2};
fclose(fileID);

max_iter = 500;
M = 1;

ifnoise = 1;
noise = RandStream('mt19937ar', 'Seed', 1);
noise.NormalTransform = 'Ziggurat';

% solvers = ["TRDS", "TRDSB", "fminsearch", "fminunc", "dfom", "dfoc", "dfo1", "para-dfo1", "para-dfo2", "dfo-plus", "subdfo", "new", "newuoam", "newuoas", "CMA-ES"];
solvers = ["dfo1"];


C = [1,0,0;...
    0,1,0;...
    0,0,1;...
    1/3,1/3,1/3;...
    1/2,1/2,0;...
    0,1/2,1/2];
c_len = size(C,1);

max_prob = size(prob, 1);
history = zeros(max_prob, c_len, M, max_iter);
yhist = zeros(max_prob, c_len, M, max_iter+1);
iter_hist = zeros(max_prob, c_len, M, 6, max_iter);


set(0, 'DefaultFigureVisible', 'off');  
t1 = clock;
for ip = 1:max_prob
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(strcat(int2str(ip), '. ', prob{ip}, '_', int2str(tstn(ip)), ':'));
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    
    func = @(x) evalfun(prob{ip}, x, tstn(ip), ifnoise * 1e-8 * randn(noise));
    [x_0_set, rhobeg, ~, ~] = setuptest (prob{ip}, tstn(ip));
    for m = 1:M
        x_0 = x_0_set * (1 + (-1)^m * 0.1 * m) + 0.05 * (-1)^(m+1) * m;
        % disp('-------------------------------------------------------');
        % disp(strcat('m=', int2str(m), ', x_0=', num2str(x_0(1)), ':'));
        % disp('-------------------------------------------------------');
        
        for ic = 1:c_len
            hist = [];
            [hist, x_dfo1, fval_dfo1, ~, ~, ihist] = dfo1_4test(max_iter, func, x_0, C(ic,:));
            history(ip, ic, m, :) = revise_hist(max_iter, hist, func(x_0));
            iter_hist(ip, ic, m, :, :) = revise_ihist(max_iter, ihist);
        end
    end
    

    PlotStyle={
        struct('Color',[1,0,0],'LineStyle','-'),...
        struct('Color',[0,1,0],'LineStyle','-'),...
        struct('Color',[0,0,1],'LineStyle','-'),...
        struct('Color',[0,0,0],'LineStyle','-'),...%    
        struct('Color',[1,1,0],'LineStyle','-'),...%yellow
        struct('Color',[1,0,1],'LineStyle','-'),...%pink
        struct('Color',[0,1,1],'LineStyle','-'),...
        struct('Color',[0.5,0.5,0.5],'LineStyle','-'),...%gray
        struct('Color',[136,0,21]/255,'LineStyle','-'),...%dark red
        struct('Color',[255,127,39]/255,'LineStyle','-'),...%orange
        struct('Color',[0,162,232]/255,'LineStyle','-'),...%Turquoise
        struct('Color',[163,73,164]/255,'LineStyle','-'),...%purple    
        struct('Color',[1,0,0],'LineStyle','--'),...
        struct('Color',[0,1,0],'LineStyle','--'),...
        struct('Color',[0,0,1],'LineStyle','--'),...
        struct('Color',[0,0,0],'LineStyle','--'),...%    
        struct('Color',[1,1,0],'LineStyle','--'),...%yellow
        struct('Color',[1,0,1],'LineStyle','--'),...%pink
        struct('Color',[0,1,1],'LineStyle','--'),...
        struct('Color',[0.5,0.5,0.5],'LineStyle','--'),...%gray
        struct('Color',[136,0,21]/255,'LineStyle','--'),...%dark red
        struct('Color',[255,127,39]/255,'LineStyle','--'),...%orange
        struct('Color',[0,162,232]/255,'LineStyle','--'),...%Turquoise
        struct('Color',[163,73,164]/255,'LineStyle','--'),...%purple    
        struct('Color',[1,0,0],'LineStyle','-.'),...
        struct('Color',[0,1,0],'LineStyle','-.'),...
        struct('Color',[0,0,1],'LineStyle','-.'),...
        struct('Color',[0,0,0],'LineStyle','-.'),...%    
        struct('Color',[1,1,0],'LineStyle',':'),...%yellow
        struct('Color',[1,0,1],'LineStyle','-.'),...%pink
        struct('Color',[0,1,1],'LineStyle','-.'),...
        struct('Color',[0.5,0.5,0.5],'LineStyle','-.'),...%gray
        struct('Color',[136,0,21]/255,'LineStyle','-.'),...%dark red
        struct('Color',[255,127,39]/255,'LineStyle','-.'),...%orange
        struct('Color',[0,162,232]/255,'LineStyle','-.'),...%Turquoise
        struct('Color',[163,73,164]/255,'LineStyle','-.'),...%purple
    };


    figstr = [prob{ip} '-' int2str(tstn(ip))];

    % yhist
    for m = 1:M
        for ic = 1:c_len
            hist_tmp = reshape(history(ip, ic, m, :), 1, size(history, 4));
            yhist(ip, ic, m, :) = [1, hist_tmp / func(x_0)];
        end
    end


    disp('--------------------------------------------------------------------');
    disp('--------------------------------------------------------------------');
    disp(['Solving problem' num2str(ip) ':' prob{ip} '_' int2str(tstn(ip)) '  Process：' num2str(ip / max_prob * 100) '%']);
    t2 = clock;
    time_consump = etime(t2, t1);
    fprintf("Total time is %.3f seconds.\n",time_consump);
    fprintf("Expected total time is %.3f seconds.\n",time_consump / (ip / max_prob * 100) * 100);
    disp('--------------------------------------------------------------------');
    disp('--------------------------------------------------------------------');
end

save yhist
cd('..')

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