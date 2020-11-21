function ID_by_category(Y_path, modeldir,result_mat)
%ID_BY_CATEGORY 
%   Calculate the average correlation coefficient based on the attack category
disp(modeldir);
cla = ["smooth", "noise", "quantization", "elementreordering"];
intensity = {["iteration5_","iteration10_","iteration30_","iteration50_"];...
    ["0.0005", "0.0010", "0.0030", "0.0050"];...
    ["11","10","9","8","7"];...
    [];};

files  = dir( modeldir );
Y = read_mesh(Y_path);

fps = {};
for i = 1 : length(files)
    if( isequal( files( i ).name, '.' )||...
        isequal( files( i ).name, '..')||...
        files( i ).isdir)               % 如果是目录则跳过
        continue;
    end
    fp = fullfile( files(i).folder, '/', files(i).name );
    fps = [fps fp];
end
len = length(fps);
disp(len);
M = containers.Map('KeyType','char','ValueType','any');

% Categorize files by attack category 
%(if there is intensity, use intensity as the smallest unit)
for j = 1 : length(cla)
    in = intensity{j};
    if ~isempty(in)
        for k = 1 : length(in)
            value = {};
            for i = 1 : len
                fp = fps(i);
                if contains(fp, cla(j)) &&  contains(fp,in(k))
                    value = [value fp];
                end
            end
            key = strcat(cla(j),'_', in(k));
            M(convertStringsToChars(key)) = value;
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
        end
    end
end
opt.debug = 0;
opt.max_it = 200; 
opt.viz = 0;
opt.segment = 0;
opt.metric = "CORR";

ks = keys(M);
keylen = length(ks);
dist = cell(keylen, 2);
for i = 1 : keylen
    key = char(ks(i));
    value = M(key);
    disp(key);
    dist{i, 1} = key;
    disp('----------------------------');
    d = zeros(1, length(value));
    for j = 1 : length(value)
        X_path = char(value(j));
        disp(X_path);
        X = read_mesh(X_path);
        try
            d(j) = Identification(X, Y, opt);
        catch
            d(j) = 0;
        end
        % d(j) = rand;
    end
    dist{i, 2} = mean(d);
    disp(['the mean of ', key, ' : ', num2str(dist{i,2})]);
    disp('----------------------------');
end

save(result_mat ,'dist');
end

