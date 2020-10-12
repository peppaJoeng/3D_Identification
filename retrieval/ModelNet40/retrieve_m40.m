clear all; close all; clc;

addpath(genpath('../../utils'));
addpath(genpath('../../core'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));

modelnet10path = '../../model/modelnet10.mat';
load(modelnet10path);

result_path = 'result_m40.csv';

d = [];
for i = 1:3
    X = read_mesh(meshes(i).path);
    clax = meshes(i).cla;
    
    dist = [];
    for j = 1:3
        if j == i
            dist(end + 1) =  0;
            continue;
        end
        Y = read_mesh(meshes(j).path);
        opt.debug = 0;
        opt.max_it = 10;
        opt.viz = 0;
        try
            dist(end+1) = Identification(X, Y, opt, '');
        catch
            dist(end+1) = 10000000000;
        end
    end
    d = [d; dist];
end

disp(size(d));
writematrix(d, result_path); 
