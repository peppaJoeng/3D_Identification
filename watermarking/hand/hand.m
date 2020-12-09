clear all; close all; clc;

addpath(genpath('../../core'));
addpath(genpath('../../utils'));
addpath(genpath('../../thirdparty/CPD2/core'));
addpath(genpath('../../thirdparty/inexact_alm_rpca'));
addpath(genpath('../../mex'));
addpath(genpath('../../watermarking'));

Y_path = '../../data/Hand/hand.off';
modeldir = '../../data/Hand/';
result_mat = 'hand.mat';


ID_by_category_4_euclidean(Y_path, modeldir,result_mat);

