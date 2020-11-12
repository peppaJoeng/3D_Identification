clear all; close all; clc;
% 
addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));
% store_name = 'xxxx';
% metric = "dsfdg";
% block_num = 1;
% dpath = strcat('distance_', store_name, '_', metric, '_', num2str(block_num), '.csv');
% disp(dpath);
% y = [1,0,0,0,0];
% writematrix(y,dpath);


source_path = '../model/frequent/bun_zipper_res4.ply';
des_path = '../model/stanford 3D/dragon_recon/dragon_vrip_res4.ply';
% des_path = '../model/frequent/bun_zipper_res3.ply';
% des_path = '../model/frequent/bun_zipper_res4.ply';
X = read_mesh(source_path);
Y = read_mesh(des_path);
% Y = X;

% % b453 b 453
% [N, D] = size(X);
% % generate target point set
% R = cpd_R(rand(1),rand(1),rand(1));
% t = rand(D, 1);
% X = rand(1) * X * R' + repmat(t', N, 1);


% % b378 b227
% X = downsample(X, 0.5);
% Y = downsample(Y, 0.2);

% b453 d521
Y = downsample(Y, 0.1);

[N, ~] = size(X);
[M, ~] = size(Y);

% 
% % option
opt.method = "rigid";
opt.viz=1;          
opt.outliers=0;     
opt.normalize=1;    
opt.scale=1;       
opt.rot=1;          
opt.corresp=1;      
opt.max_it=30;
opt.tol=1e-8;

[Transform,correspondence] = cpd_register(X, Y, opt);

%xxx = [255/255, 255/255, 240/255;0,0.447000000000000,0.741000000000000;0.850000000000000,0.325000000000000,0.0980000000000000;255/255, 20/255,147/255;0.929000000000000,0.694000000000000,0.125000000000000;0.494000000000000,0.184000000000000,0.556000000000000;0.466000000000000,0.674000000000000,0.188000000000000;0.301000000000000,0.745000000000000,0.933000000000000;139/255, 71/255, 38/255;0.635000000000000,0.0780000000000000,0.184000000000000; 84/255, 255/255, 159/255;255/255,193/255,193/255;0,0,0];
opt.debug=0;
[X, Y] = sort_By_Longest_Axis(X, Transform.Y, opt);
% write_off('X.off',X);
% write_off('Y.off',Y);

[XN, YN, ~]=cpd_normalize(X, Y); 
P = postProbablity(XN, YN, Transform.sigma2, opt.outliers);

plot(P(200,:))
% picname = strcat('bun',num2str(M),' x bun', num2str(N));
% plot_heatmap(P, picname);
% h = heatmap(P(1:200,1:200),'Colormap',parula(5), 'GridVisible','off');
% h = heatmap(P,'Colormap',parula(5), 'GridVisible','off');
%load('MyColormaps.mat', 'mymap');
% h = heatmap(P,'Colormap',xxx, 'GridVisible','off');
% %h.Colormap = repmat(linspace(0, 1, 255).', 1, 3);
% mmm = h.Colormap;
% % %h.XData = {};
% %h.YData = {};
% % h.Title = 'Posterior Probability Map of 2 coincident  PS';
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

%mymap = get(gcf,'Colormap');%gcf是get current figure的缩写
%save('MyColormaps','mymap');%把mymap变量保存为MyColormaps.mat，位置在matlab当前目录
% heatmap(P,'Colormap',parula(5));
% sss=sum(P(:,100));
% figure('Name','1');
% plot(P(:,100));
% PT = P';
% lll = zeros(1, N);
% for i = 1 : N
%      lll(i) = kurtosis(PT(i,:));
% end
% bbb = kurtosis(P);
% x = PT(1,:);
% u = mean(x);
% d = var(x);
% disp(u);
% disp(d);
% disp(1/N);
% disp((N-1)/(N*N));
% disp(lll(1));
% NN = (N - 1) * (N - 1);
% disp((NN + 1 / NN) / N);
% if opt.corresp, C=cpd_Pcorrespondence(XN, YN,Transform.sigma2,opt.outliers); else C=0; end;
% 
% %CC = repmat(C,1,3);
% C=cpd_Pcorrespondence(XN, YN,Transform.sigma2,opt.outliers);
% CC = repmat(C,1,3);
% % 
% 
% YY = zeros(M,3); 
% for i = 1 : M
%     %disp(C(i));
%     YY(i, :) = X(C(i), :);
%     %disp(X(C(i), :));
%     %disp(Y(C(i)));
%     %disp(YY(i));
% end
% R = corr2(Y, YY);
% figure('Name', 'hhh');
% cpd_plot_iter(X, Y);
% figure('Name', 'myname');
% cpd_plot_iter(YY, Y);
% write_off('X.off',X);
% write_off('Y.off',Y);

% 
% figure('Name','Before');
% cpd_plot_iter(X, Y); %title('Before');
% figure('Name','After registering Y to X');
% cpd_plot_iter(X, Transform.Y); % title('After registering Y to X');
% 
% tic;
% [AA,EE,ITER] = inexact_alm_rpca(P');
% figure('name','A');
% plot(AA(:,1));
% disptime(toc);

% 
% 
% % aa = PT(P, C);
% % bb = PF(P);
% % 
% % rho = corr(aa(:), bb(:), 'type','pearson');
picx = strcat('bun',num2str(M),' x bun', num2str(N));
plot_heatmap(P, picx);

tic;
[AA1,EE1,ITER1] = inexact_alm_rpca(P);
figure('name','AA');
plot(AA1(:,1));
disptime(toc);

picname = strcat('A bun',num2str(M),' x bun', num2str(N));
plot_heatmap(AA1, picname);
% 
% 
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
    disp(['the longest axis is ', num2str(axis), ' whose length is ', num2str(max_length)]);
    if max_length == 1
        X_new = sortrows(X);
        Y_new = sortrows(Y);
    else
        X_new = sortrows(X, [axis, 1]);
        Y_new = sortrows(Y, [axis, 1]);
    end
    
end
% 
% 
% function res = PT(P, corres)
%     [M, N] = size(P);
%     res = zeros(M, N);
%     for i = 1 : M
%         res(i, corres(i)) = 1;
%     end
% end
% 
% function res = PF(P)
%     res = P >= 0.5;
% end

