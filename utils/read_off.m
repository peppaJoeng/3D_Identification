function [vertex, face] = read_off(filename)
% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2003 Gabriel Peyr?
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end
% ��ȡָ���ļ��е���һ�����ݣ����������з���
str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end
str = fgets(fid); 
% �����ҽ��� str��ʹ�ÿհ��ַ���Ϊ�ָ��������� a �з��ز��ֻ�ȫ���ı������strΪʣ�ಿ��
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);
% A = fscanf(fileID,formatSpec,sizeA) ���ļ����ݶ�ȡ��ά��Ϊ sizeA ������ A �У������ļ�ָ�붨λ������ȡ��ֵ֮��
% fscanf ����˳����� A��
% sizeA ����Ϊ����������� [m n] ����ʽ������ m �� n Ϊ��������
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A';
% read Face 1  1088 480 1022
[A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
if cnt~=4*nface
    warning('Problem in reading faces.');
end
A = reshape(A, 4, cnt/4);
face = A(2:4,:)+1;
fclose(fid);
end

