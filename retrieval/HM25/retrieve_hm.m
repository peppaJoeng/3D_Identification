function retrieve_hm(library_path, query_path, store_name, block, block_num, metric)
    addpath(genpath('../../utils'));
    addpath(genpath('../../core'));
    addpath(genpath('../../thirdparty/CPD2/core'));
    addpath(genpath('../../thirdparty/inexact_alm_rpca'));
    addpath(genpath('../../mex'));

    
    load(library_path, 'meshes');
    libs = meshes; 
    load(query_path, 'meshes');
    queries = meshes;
    
    [~, N] = size(queries.path);
    [~, M] = size(libs.path);
    
    dpath = strcat('distance_', store_name, '_', metric, '_', num2str(block_num), '.csv');
    Ipath = strcat('I_', store_name, '_', metric, '_', num2str(block_num), '.csv');
    ppath = strcat('path_', store_name, '_', metric, '_', num2str(block_num), '.csv');

    disp(['N : ', num2str(N), ' and M :', num2str(M)]);
    n = N / block;
    d = zeros(n, M);
    start = n * (block_num - 1);
    last = n * block_num;
    disp([num2str(start + 1), ' and ', num2str(last)]);
    % retrieval
    for i = 1 : n
        % query as X
        disp(queries.path(start + i));
        Y = read_mesh(queries.path(start + i));
        % compute M results
        dist = zeros(1, M);
        for j = 1 : M
            disp([num2str(start + i), ' and ', num2str(j)]);
            disp(libs.path(j));
            X = read_mesh(libs.path(j));
            opt.debug = 0;
            opt.max_it = 200;
            opt.viz = 0;
            opt.segment = 0;
            opt.metric = metric;
            try
                if opt.metric == "LR"
                    dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
                else
                    dist(j) = Identification(X, Y, opt);
                end
            catch
                dist(j) = 10000000;
            end
        end
        d(i, :) = dist;
    end

    writematrix(d, dpath); 
    
    if opt.metric == "LR"
        [~, I] = sort(d, 2);
    else
        [~, I] = sort(d, 2, 'descend');
    end

    writematrix(I, Ipath);
    pathmat = strings(n, M + 1);

    for i = 1 : n
        pathmat(i, 1) = queries.path(start + i);
        pathmat(i, 2 : M + 1) = libs.path(I(i, 1 : M));
    end
    writematrix(pathmat, ppath); 
end

