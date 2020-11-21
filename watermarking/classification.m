clear all; close all; clc;

cla = ["smooth", "noise", "quantization", "elementreordering"];
intensity = {["iteration5_","iteration10_","iteration30_","iteration50_"];...
    ["0.0005", "0.0010", "0.0030", "0.0050"];...
    ["11","10","9","8","7"];...
    [];};
modledir = '../data/Bunny/';
mid_path = '../model/';
files  = dir( modledir );

S = struct([]);

fps = {};
for i = 1 : length(files)
    if( isequal( files( i ).name, '.' )||...
        isequal( files( i ).name, '..')||...
        files( i ).isdir)               % 如果不是目录则跳过
        continue;
    end
    fp = fullfile( files(i).folder, '/', files(i).name );
    fps = [fps fp];
end
len = length(fps);
disp(len);
M = containers.Map('KeyType','char','ValueType','any');

for j = 1 : length(cla)
    in = intensity{j};
    if ~isempty(in)
        for k = 1 : length(in)
            value = {};
            for i = 1 : len
                fp = fps(i);
                
                if contains(fp, cla(j)) &&  contains(fp,in(k))
                    % disp(fp);
                    value = [value fp];
                end
            end
            key = strcat(cla(j),'_', in(k));
            M(convertStringsToChars(key)) = value;
%             disp(key)
%             disp(value)
        end
    else
        value = {};
        for i = 1 : len
            fp = fps(i);
            if contains(fp, cla(j))
                value = [value fp];
            end
            key = cla(j);
            M(convertStringsToChars(key)) = value;
%             disp(key)
%             disp(value)
        end
    end
end

ks = keys(M);
keylen = length(ks);
dist = cell(keylen, 2);
for i = 1 : keylen
    key = char(ks(i));
    value = M(key);
    disp(key);
    dist{i, 1} = key;
    disp('+++++++++++++++++++++++++++++');
    d = zeros(1, length(value));
    for j = 1 : length(value)
        Y = value(j);
        disp(Y);
        d(j) = rand;
    end
    disp(d);
    dist{i, 2} = mean(d);
    disp('----------------------------');
end
result_mat = 'bunny.mat';
save(result_mat ,'dist');
