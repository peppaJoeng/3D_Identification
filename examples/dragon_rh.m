clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));

source_path = '../model/xyzrgb_dragon-100k.ply';
result_dir = preprocess(source_path);

diary([result_dir, '/log.txt']);
diary on;

X = read_mesh(source_path);
X = downsample(X,0.01);
Y = X;

[N, D] = size(X); 
t = [0.1;0.1;0.1];
s = 0.3;
X = s * X + repmat(t', N, 1);

opt.max_it = 20;
opt.debug = 1;
opt.viz = 0; 
distance = Identification(X, Y, opt, result_dir);

diary off;
