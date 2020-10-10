clear all; close all; clc;

addpath(genpath('../../../core'));
addpath(genpath('../../../ModelNet10'));
addpath(genpath('../../../../inexact_alm_rpca'));

% N > M
source_path = '../model/Merlion2.obj';
result_dir = ['result_s/', source_path(10:end-4)];
disp(result_dir);

if ~exist('./result_s','dir')
    mkdir('./result_s')
end
if ~exist(result_dir,'dir')
    mkdir(result_dir)
end
if exist([result_dir, '/log.txt'], 'file')
    delete([result_dir, '/log.txt']);
end
diary([result_dir, '/log.txt']);
diary on;

[X, ~] = read_obj(source_path);
Y = X;
disp(X);

R = cpd_R(rand(1),rand(1),rand(1));
[N, D] = size(X); 
t = rand(D, 1);
X = rand(1) * X * R' + repmat(t', N, 1);

X = HGMM(X, 3, [result_dir,'/X'], true, 0.9, 2.5, 3);
Y = HGMM(Y, 3, [result_dir,'/Y'], true, 0.9, 2.5, 3);
disp(size(X));
disp(size(Y));

opt.method = 'rigid';
opt.viz = 1;
opt.outliers = 0;
opt.normalize = 1;
opt.scale = 1;
opt.rot =1;
opt.corresp = 1;
opt.max_it = 200;
opt.tol = 1e-8;
% transformation
[Transform, ~] = cpd_register(X,Y,opt);

before = figure('Name','Before');
cpd_plot_iter(X, Y); 
saveas(before,[result_dir,'/before'],'jpg');

after = figure('Name','After registering Y to X');
cpd_plot_iter(X, Transform.Y);
saveas(after,[result_dir,'/after'],'jpg');

thred = 4000;
splitModel(X, Transform.Y, opt, result_dir, Transform.sigma2, thred);

diary off;