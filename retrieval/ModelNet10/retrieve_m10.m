clear all; close all; clc;

addpath(genpath('../../utils'));
addpath(genpath('../../core'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));

modelnet10path = '../../model/modelnet.mat';
load(modelnet10path);

d_result_path = 'dis_result_m10_3.csv';
p_result_path = 'path_result_m10_3.csv';

d = zeros(3, 3);
for i = 1:3
    X = read_mesh(meshes.path(i));
    clax = meshes.cla(i);
    
    dist = zeros(1, 3);
    for j = 1:3
        if j == i
            dist(j) =  0;
            continue;
        end
        Y = read_mesh(meshes.path(j));
        opt.debug = 0;
        opt.max_it = 130;
        opt.viz = 0;
        opt.split = 0;
        % only if set split flag
        opt.thred = 2000;
        try
            dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
        catch
            dist(j) = 10000000000;
        end
    end
    d(i, :) = dist;
end

writematrix(d, d_result_path); 

[d, I] = sort(d,2);
pathd = [];
path = [];
[DM, DN] = size(d);
for i = 1 : DM
    II = I(i,:);
    path = meshes.path(II);
    pathd = [pathd; path];
end
writematrix(pathd, p_result_path); 
