clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

% N > M
source_path = '../data/Horse/horse_noise_intensity0.0050_keepboundary_number3.off';
des_path = '../data/Horse/horse.off';

result_dir = 'horse_h2';
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end

diary([result_dir,'/simp_horse_h2.log']);
diary on;

opt.max_it = 200;
opt.debug = 0;
opt.viz = 0; 
opt.segment = 0;
opt.metric = "ALL";

X = read_mesh(source_path);
Y = read_mesh(des_path);


disp('=======HGMM2 models=========');
tic;
X = HGMM(X, 20,[result_dir,'/X'], 40);
disp(['HGMM X timeval : ', num2str(toc)]);
tic;
Y = HGMM(Y, 20,[result_dir,'/Y'], 40);
disp(['HGMM Y timeval : ', num2str(toc)]);

disp(size(X));
disp(size(Y));
opt.savename = "horse_h2.mat";
distance = Identification(X, Y, opt, result_dir);

diary off;
