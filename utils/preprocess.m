function result_dir = preprocess(path)
%PREPROCESS 此处显示有关此函数的摘要
%   此处显示详细说明
pathsave = pwd;
utils_path = mfilename('fullpath'); [pathstr, ~, ~] = fileparts(utils_path);
cd (pathstr);
cd ..;
[~, name, ~] = fileparts(path);
% result_dir : used to store intermediate and final results
result_dir = [pathstr,'/../result/', name];

if exist(result_dir,'dir')
    status = rmdir(result_dir, 's');
    if status ~= 1
        disp('can not remove dir');
        exit(0);
    end
end
mkdir(result_dir);


% log : used to record the process
if exist([result_dir, '/log.txt'], 'file')
    delete([result_dir, '/log.txt']);
end
cd (pathsave);
end

