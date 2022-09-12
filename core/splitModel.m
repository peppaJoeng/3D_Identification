function result = splitModel(X, Y, opt, sigma2, result_dir)
%% preprocess
paths.rpca_path = 0;
if opt.debug
    paths.off_path = [result_dir, '/off'];    % store off file
    paths.pic_path = [result_dir, '/pic'];    % store pictures
    paths.X_path = [paths.off_path, '/X'];          % store off file of blocks of X
    paths.Y_path = [paths.off_path, '/Y'];          % store off file of blocks of Y
    paths.rpca_path = [paths.pic_path, '/inexact']; % store results of rpca 
    paths.block_path = [paths.pic_path, '/block'];  % store results of block
    checkAndMkdir(paths.off_path);
    checkAndMkdir(paths.pic_path);
    checkAndMkdir(paths.X_path);
    checkAndMkdir(paths.Y_path);
    checkAndMkdir(paths.rpca_path);
    checkAndMkdir(paths.block_path);
end

%% split the point cloud along the logest axis
% 选取X，Y，Z轴中点覆盖长度最长的轴并以该轴为中心排序
[X, Y] = sort_By_Longest_Axis(X, Y, opt);

% 将模型均匀分为num份，其中每一份点数最大不超过thred
[N, ~] = size(X);
[M, ~] = size(Y);
% 是否需要切割
if opt.segment
    num = ceil(max(N, M) / opt.thred);
else
    num = 1;
end

n = floor(N / num);
m = floor(M / num);
if opt.debug
    disp(['These models are splitted to ',num2str(num), ' blocks.']);
    disp(['Each block of X has ', num2str(n), ' points.'])
    disp(['Each block of Y has ', num2str(m), ' points.'])
end

%% compute rpca result each block
result = zeros(1, num);
for i = 1 : num - 1
    x = X((i - 1) * n + 1 :  n * i, :, :);
    y = Y((i - 1) * m + 1 :  m * i, :, :);
    result(i) = compute_result(x, y, sigma2, opt, paths, i);
end

x = X((num - 1) * n + 1 : end, :, :);
y = Y((num - 1) * m + 1 : end, :, :);
result(num) = compute_result(x, y, sigma2, opt, paths, num);
end

%% compute result
function dis = compute_result(x, y, sigma2, opt, paths, num)
    switch opt.metric
        case "LR"
            dis = LR_result(x, y, sigma2, opt, paths, num);
            disp(['using LR : ', num2str(dis)]);
        case "KURT"
            dis = kurto_result(x, y, sigma2, opt);
            disp(['using KURT : ', num2str(dis)]);
        case "CORR"
            dis = corr_result(x, y, sigma2, opt);
            disp(['using CORR : ', num2str(dis)]);
        case "ALL"
            dis_LR = LR_result(x, y, sigma2, opt, paths, num);
            dis_KURT = kurto_result(x, y, sigma2, opt);
            dis_CORR = corr_result(x, y, sigma2, opt);
            
            disp(['using LR : ', num2str(dis_LR)]);
            disp(['using KURT : ', num2str(dis_KURT)]);
            disp(['using CORR : ', num2str(dis_CORR)]);
            save(opt.savename, 'dis_LR', 'dis_KURT', 'dis_CORR');
            
            dis = -789456123;
        case "KC"
            dis_KURT = kurto_result(x, y, sigma2, opt);
            dis_CORR = corr_result(x, y, sigma2, opt);
            
            disp(['using KURT : ', num2str(dis_KURT)]);
            disp(['using CORR : ', num2str(dis_CORR)]);
            save(opt.savename, 'dis_KURT', 'dis_CORR');
            
            dis = -987654321;
        otherwise
            warning('Unexpected metric type.');
            return;
    end
end


%% low rank function
function dis = LR_result(x, y, sigma2, opt, paths, num)
    if opt.debug
        % 将xy分块比较结果保存为jpg图片格式
        % savePic(num2str(k), [block_path, '/part_', num2str(k)], x, y);
        % 将xy分块比较结果保存到fig文件中
        %block_path = path.block_path;
        saveFig(num2str(num), [paths.block_path, '/part_', num2str(num)], x, y);
        write_off([paths.X_path, '/part_', num2str(num), '.off'], x);
        write_off([paths.Y_path, '/part_', num2str(num), '.off'], y);
    end
    p = resbonsibility(x, y, sigma2, opt);
    A= rpca(p, num, paths.rpca_path, opt);
    dis = compute_dis(A);
end

function A = rpca(p, i, rpca_path, opt)
    disp('inexact_alm_rpca');
    tic;
    [A, ~, ~] = inexact_alm_rpca(p');
    if opt.debug
        inexact = figure('name',['A_', num2str(i)]);
        plot(A(:,1));
        saveas(inexact,[rpca_path,'/part_', num2str(i)],'jpg');
    end
    disptime(toc);
end

%% kurtosis
function dis = kurto_result(x, y, sigma2, opt)
    p = resbonsibility(x, y, sigma2, opt);
    [~, N] = size(p);
    kur = kurtosis(p);
    dis = sum(kur) / N;
end


%% corresponse
function dis = corr_result(x, y, sigma2, opt)
    C = correspondence(x, y, sigma2, opt);
    [M, ~] = size(y);
    
    yy = zeros(M,3);
    for i = 1 : M
        yy(i, :) = x(C(i), :);
    end
    dis = corr2(y, yy);
end

%% save fucntion
function savePic(pic_name, store_name, x, y)
    pic = figure('Name',pic_name);
    disp(['===========The ', pic_name, 'th part==========']);
    disp(['x : ', num2str(size(x)), ' y : ', num2str(size(y))]);
    cpd_plot_iter(x, y); 
    saveas(pic,store_name,'jpg');
end

function saveFig(pic_name, store_name, x, y)
    pic = figure('Name',pic_name);
    disp(['===========The ', pic_name, 'th part==========']);
    disp(['x : ', num2str(size(x)), ' y : ', num2str(size(y))]);
    cpd_plot_iter(x, y); 
    saveas(pic,store_name,'fig');
end

%% prepare function
function checkAndMkdir(dirname)
    if exist(dirname,'dir')
        rmdir(dirname, 's');
    end
    mkdir(dirname);
end

%% spilt fucntion
function [X_new, Y_new] = sort_By_Longest_Axis(X, Y, opt)
    if opt.debug
        disp(['max value of X coordiantes : ', num2str(max(X))]);
        disp(['min value of X coordiantes : ', num2str(min(X))]);
        disp(['length of X coordiantes in 3 axes: ', num2str(max(X) - min(X))]);
    
        disp(['max value of Y coordiantes : ', num2str(max(Y))]);
        disp(['min value of Y coordiantes : ', num2str(min(Y))]);
        disp(['length of Y coordiantes in 3 axes: ', num2str(max(Y) - min(Y))]);
    end
    [max_length, axis] = max(max(X) - min(X));
    [y_length, y_axis] = max(max(Y) - min(Y));
    if max_length < y_length
        max_length = y_length;
        axis = y_axis;
    end
%     disp(['the longest axis is ', num2str(axis), ' whose length is ', num2str(max_length)]);
    if max_length == 1
        X_new = sortrows(X);
        Y_new = sortrows(Y);
    else
        X_new = sortrows(X, [axis, 1]);
        Y_new = sortrows(Y, [axis, 1]);
    end
    
end

%% compute probability
function P = resbonsibility(X, Y, sigma2, opt)
    % X and Y must normalize to zero mean and unit variance
    [X, Y, ~]=cpd_normalize(X,Y); 
    % comupte resbonsibility
    P = postProbablity(X, Y, sigma2, opt.outliers);
end

%% compute correspondence
function C = correspondence(X, Y, sigma2, opt)
    % X and Y must normalize to zero mean and unit variance
    [X, Y, ~]=cpd_normalize(X,Y); 
    % comupte resbonsibility
    C=cpd_Pcorrespondence(X, Y, sigma2,opt.outliers);
end
