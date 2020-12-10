clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

source_path = '../data/Merlion/Merlion_similaritytransform_number1.obj';
des_path = '../data/Merlion/Merlion.obj';

result_dir = 's5000';
mkdir(result_dir);
diary('log_5000.txt');
diary on;

X = read_mesh(source_path);
Y = read_mesh(des_path);

disp(size(X));
disp(size(Y));

opt.max_it = 200;
opt.debug = 1;
opt.viz = 0; 
opt.segment = 1;
opt.metric = "LR";
opt.thred = 5000;

distance = Identification(X, Y, opt, result_dir);

diary off;
