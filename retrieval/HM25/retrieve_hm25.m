clear all; close all; clc;

addpath(genpath('../../utils'));
addpath(genpath('../../core'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));

hm25path = '../../model/mHM25.mat';
load(hm25path);

[~, meshes_num] = size(meshes.path);
d_result_path = ['dis_result_hm25_', num2str(meshes_num), '.csv'];
p_result_path = ['path_result_hm25_', num2str(meshes_num), '.csv'];


d = zeros(meshes_num, meshes_num);
for i = 1 : meshes_num
    Y = read_mesh(meshes.path(i));
    clax = meshes.cla(i);
    
    dist = zeros(1, meshes_num);
    for j = 1 : meshes_num
        if j == i
            continue;
        end
        X = read_mesh(meshes.path(j));
        opt.debug = 0;
        opt.max_it = 150;
        opt.viz = 0;
        opt.split = 0;
        try
            dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
        catch
            dist(j) = 100000000;
        end
    end
    d(i, :) = dist;
end

writematrix(d, d_result_path); 

[d, I] = sort(d,2);
pathmat = strings(meshes_num, meshes_num);

for i = 1 : meshes_num
    pathmat(i ,:) = meshes.path(I(i,:));
end

writematrix(pathd, p_result_path); 
