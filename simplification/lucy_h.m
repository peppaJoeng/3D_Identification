clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

% N > M
source_path = '../model/frequent/lucy_million_rst.ply';
des_path = '../model/frequent/lucy_million.ply';
result_dir = preprocess(source_path);

diary([result_dir, '/log_hgmm.txt']);
diary on;

X = read_mesh(source_path);
Y = read_mesh(des_path);

X = HGMM(X, 20,[result_dir,'/X'], 80 );
Y = HGMM(Y, 20,[result_dir,'/Y'], 90);

disp(size(X));
disp(size(Y));

opt.max_it = 200;
opt.debug = 1;
opt.viz = 0; 
opt.segment = 0;
opt.metric = "ALL";
opt.savename = "lucy_hgmm.mat";

distance = Identification(X, Y, opt, result_dir);

diary off;
