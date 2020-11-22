clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

% N > M
source_path = '../model/frequent/Merlion2.obj';
des_path = '../model/frequent/Merlion2.obj';
result_dir = preprocess(source_path);

diary([result_dir, '/log_base.txt']);
diary on;

X = read_mesh(source_path);
Y = read_mesh(des_path);

disp(size(X));
disp(size(Y));

opt.max_it = 200;
opt.debug = 1;
opt.viz = 0; 
opt.segment = 0;
opt.metric = "LR";

distance = Identification(X, Y, opt, result_dir);

diary off;
