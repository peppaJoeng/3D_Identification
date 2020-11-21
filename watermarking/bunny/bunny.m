clear all; close all; clc;

addpath(genpath('../../core'));
addpath(genpath('../../utils'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));
addpath(genpath('../../watermarking'));

Y_path = '../../data/Bunny/bunny.off';
modeldir = '../../data/Bunny/';
result_mat = 'bunny.mat';


ID_by_category(Y_path, modeldir,result_mat);

