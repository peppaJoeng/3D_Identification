function distance = Identification(X, Y, opt, result_dir)
%IDENTIFICATION 此处显示有关此函数的摘要
%   此处显示详细说明
addpath(genpath('../thirdparty/CPD2/core'));
addpath(genpath('../thirdparty/inexact_alm_rpca'));
addpath(genpath('../mex'));

if nargin < 4 && opt.debug == 0
    result_dir = '';
end

opt.method = 'rigid';
opt.outliers = 0;
opt.normalize = 1;
opt.scale = 1;
opt.rot =1;
opt.corresp = 1;
opt.tol = 1e-8;

% transformation
[Transform, ~] = cpd_register(X,Y,opt);

result = splitModel(X, Transform.Y, opt, Transform.sigma2, result_dir);

[~, len] = size(result);
% disp(['Now aggregate the ', num2str(len), ' block results']);
distance = 0;
for i = 1 : len
    distance = distance + result(i);
end

end



