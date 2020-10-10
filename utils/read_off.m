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
% 读取指定文件中的下一行内容，并包含换行符。
str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end
str = fgets(fid); 
% 从左到右解析 str，使用空白字符作为分隔符，并在 a 中返回部分或全部文本。输出str为剩余部分
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);
% A = fscanf(fileID,formatSpec,sizeA) 将文件数据读取到维度为 sizeA 的数组 A 中，并将文件指针定位到最后读取的值之后。
% fscanf 按列顺序填充 A。
% sizeA 必须为正整数或采用 [m n] 的形式，其中 m 和 n 为正整数。
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

