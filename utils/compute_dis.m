function distance = compute_dis(rpca)
%COMPUTE_DIS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
dis = rpca(:,1);
dist = abs(dis);
distance = sum(dist);
end

