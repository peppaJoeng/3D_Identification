function X = downsample(X, percentage)
%DOWNSAMPLE 此处显示有关此函数的摘要
%   此处显示详细说明
pcx = pointCloud(X);
pcx = pcdownsample(pcx,'random',percentage);
X = pcx.Location;
end

