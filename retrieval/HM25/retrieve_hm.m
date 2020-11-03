function retrieve_hm(library_path, query_path, store_name, block, block_num)
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
    dpath = ['distance_', store_name, '_', num2str(block_num), '.csv'];
    ppath = ['path_', store_name, '_', num2str(block_num), '.csv'];

    disp(['N : ', num2str(N), ' and M :', num2str(M)]);
    d = zeros(floor (N / block), M);
    start = N / block * (block_num - 1)  + 1;
    last = N / block * block_num;
    % retrieval
    for i = start : last
        % query as X
        X = read_mesh(queries.path(i));
        % compute M results
        dist = zeros(1, M);
        for j = 1 : M
            Y = read_mesh(libs.path(j));
            opt.debug = 0;
            opt.max_it = 200;
            opt.viz = 0;
            opt.split = 0;
            try
                dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
            catch
                dist(j) = 10000000;
            end
        end
        d(i, :) = dist;
    end

    writematrix(d, dpath); 

    [~, I] = sort(d,2);
    pathmat = strings(N, M);

    for i = 1 : N
        pathmat(i, 1 : M) = libs.path(I(i, 1 : M));
    end
    writematrix(pathmat, ppath); 
end

