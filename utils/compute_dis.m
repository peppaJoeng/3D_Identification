function distance = compute_dis(rpca)
%COMPUTE_DIS 此处显示有关此函数的摘要
%   此处显示详细说明
dis = rpca(:,1);
dist = abs(dis);
distance = sum(dist);
end

