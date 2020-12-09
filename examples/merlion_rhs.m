clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../mex'));

% N > M
source_path = '../model/frequent/Merlion2.obj';
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

X = read_mesh(source_path);
Y = X;

R = cpd_R(rand(1),rand(1),rand(1));
[N, D] = size(X); 
t = rand(D, 1);
X = rand(1) * X * R' + repmat(t', N, 1);

% X = HGMM(X, 3,[result_dir,'/X'], 3);
Y = HGMM(Y, 3,[result_dir,'/Y'], 3.5);
%write_off('Merlion_hX.off',X);
write_off('Merlion_hY.off',Y);
disp(size(X));
disp(size(Y));
%{
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
%}
diary off;