function X = downsample(X, percentage)
%DOWNSAMPLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
pcx = pointCloud(X);
pcx = pcdownsample(pcx,'random',percentage);
X = pcx.Location;
end

