clear all; close all; clc;

addpath(genpath('../../utils'));
addpath(genpath('../../core'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));

modelnet10path = '../../model/modelnet.mat';
load(modelnet10path);
d_result_path = ['dis_result_m10_3.csv'];
disp(d_result_path);
p_result_path = 'path_result_m10_3.csv';

d = zeros(3, 3);
for i = 1:3
    Y = read_mesh(meshes.path(i));
    clax = meshes.cla(i);
    
    dist = zeros(1, 3);
    for j = 1:3
%         if j == i
%             dist(j) =  0;
%             continue;
%         end
        X = read_mesh(meshes.path(j));
        opt.debug = 0;
        opt.max_it = 130;
        opt.viz = 0;
        opt.segment = 0;
        opt.metric = "CORR";
%         try
            %dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
            dist(j) = Identification(X, Y, opt);
%         catch
%             dist(j) = 10000000000;
%         end
    end
    d(i, :) = dist;
end

writematrix(d, d_result_path); 

if opt.metric == "LR"
    [d, I] = sort(d, 2);
else
    [d, I] = sort(d, 2, 'descend');
end

[M, n] = size(d);

%M = 3;
%n = 3;
pathmat = strings(n, M + 1);

for i = 1 : M
    pathmat(i, 1) = meshes.path(i);
    pathmat(i, 2 : M + 1) = meshes.path(I(i, 1 : M));
end
writematrix(pathmat, p_result_path); 
