clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));

% N > M
source_path = '../model/xyzrgb_dragon-100k.ply';
result_dir = preprocess(source_path);

diary([result_dir, '/log.txt']);
diary on;

X = read_mesh(source_path);
Y = X;

R = cpd_R(rand(1),rand(1),rand(1));
[N, D] = size(X); 
t = rand(D, 1);
X = rand(1) * X * R' + repmat(t', N, 1);

X = HGMM(X, 20,[result_dir,'/X'], 19);
Y = HGMM(Y, 20,[result_dir,'/Y'], 25);

disp(size(X));
disp(size(Y));

opt.max_it = 200;
opt.debug = 1;
opt.viz = 0; 
distance = Identification(X, Y, opt, result_dir);

diary off;
