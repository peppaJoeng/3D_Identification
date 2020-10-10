clear all; close all; clc;

addpath(genpath('../../../core'));
addpath(genpath('../../../ModelNet10'));
addpath(genpath('../../../../inexact_alm_rpca'));

% N > M
source_path = '../model/Armadillo.ply';
result_dir = ['result/', source_path(10:end-4)];
disp(result_dir);

if ~exist('./result','dir')
    mkdir('./result')
end
if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if exist([result_dir, '/log.txt'], 'file')
    delete([result_dir, '/log.txt']);
end

diary([result_dir, '/log.txt']);
diary on;

[X, ~] = read_ply(source_path);
Y = X;

R = cpd_R(rand(1),rand(1),rand(1));
[N, D] = size(X); 
t = rand(D, 1);
X = rand(1) * X * R' + repmat(t', N, 1);

% disp(size(X));
X = HGMM(X, 15, [result_dir,'/X'], false);
Y = HGMM(Y, 30, [result_dir,'/Y'], false);
disp(size(X));
disp(size(Y));

opt.method = 'rigid';
opt.viz = 1;
opt.outliers = 0;
opt.normalize = 1;
opt.scale = 1;
opt.rot =1;
opt.corresp = 1;
opt.max_it = 10;
opt.tol = 1e-8;
% transformation
[Transform, Correspondence] = cpd_register(X,Y,opt);

P = Transform.p;

before = figure('Name','Before');
cpd_plot_iter(X, Y); 
saveas(before,[result_dir,'/before'],'jpg');

after = figure('Name','After registering Y to X');
cpd_plot_iter(X, Transform.Y);
saveas(after,[result_dir,'/after'],'jpg');

disp('inexact_alm_rpca');
tic;
[A12, E12, ~] = inexact_alm_rpca(P');
inexact = figure('name','A12');
plot(A12(:,1));
saveas(inexact,[result_dir,'/inexact'],'jpg');
disptime(toc);

diary off;
