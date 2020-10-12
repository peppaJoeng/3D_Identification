function [pc] = read_mesh(path)
%READ_MESH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[~, ~, ext] = fileparts(path);
if ext == '.off'
    [pc, ~] = read_off(path);
elseif ext == '.ply'
    [pc, ~] = read_ply(path);
elseif ext == '.obj'
    [pc, ~] = read_obj(path);
else
    disp('System can only support .off/.ply/.obj file as input');
end

