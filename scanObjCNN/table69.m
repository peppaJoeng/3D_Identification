clear all; close all; clc;

addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

% N > M
source_path = './real/object_dataset_table_069_00004.off';
des_path = './real/PB_T50_RS_table_069_00004_4.off';

[path,name,ext]=fileparts(source_path);
nameList=strsplit(name, '_');

name=strcat(nameList(3), nameList(4), nameList(5));
result_dir=string(name);
disp(result_dir);
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end

logfile=strcat(result_dir, "/", name, ".log");
disp(logfile);
diary(logfile);
diary on;

opt.max_it = 200;
opt.debug = 0;
opt.viz = 0;
opt.segment = 0;
opt.metric = "ALL";

X = read_mesh(source_path);
Y = read_mesh(des_path);
X = downsample(X, 0.07);
Y = downsample(Y, 0.1);

disp('=======origin models=========');
disp(size(X));
disp(size(Y));
opt.savename = strcat(name, ".mat");
distance = Identification(X, Y, opt, result_dir);

diary off;
