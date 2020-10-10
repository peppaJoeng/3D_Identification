function distance = Identification(X, Y, opt, result_dir)
%IDENTIFICATION 此处显示有关此函数的摘要
%   此处显示详细说明
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

opt.method = 'rigid';
opt.outliers = 0;
opt.normalize = 1;
opt.scale = 1;
opt.rot =1;
opt.corresp = 1;
opt.tol = 1e-8;

% transformation
[Transform, ~] = cpd_register(X,Y,opt);


opt.thred = 2000;
result = splitModel(X, Transform.Y, opt, Transform.sigma2, result_dir);

[~, len] = size(result);
disp(len);
distance = 0;
for i = 1 : len
    distance = distance + compute_dis(result{1,i});
end

end



