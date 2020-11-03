function retrieve_hm5(block_num)

    addpath(genpath('../../utils'));
    addpath(genpath('../../core'));
    addpath(genpath('../../thirdparty/CPD2/core'));
    addpath(genpath('../../thirdparty/inexact_alm_rpca'));
    addpath(genpath('../../mex'));

    hm5path = '../../model/HM5.mat';
    load(hm5path);
    
    [~, meshes_num] = size(meshes.path);
    d_result_path = ['dis_result_hm5_', num2str(meshes_num), '_', num2str(block_num), '.csv'];
    p_result_path = ['path_result_hm5_', num2str(meshes_num), '_', num2str(block_num), '.csv'];

    disp(meshes_num);
    d = zeros(meshes_num / 10, meshes_num);
    s = meshes_num / 10 * (block_num - 1)  + 1;
    e = meshes_num / 10 * block_num;
    disp([num2str(s), ' ', num2str(e)]);
    for i = s : e
        Y = read_mesh(meshes.path(i));
        % clax = meshes.cla(i);
    
        dist = zeros(1, meshes_num);
        for j = 1 : meshes_num
            if j == i
                continue;
            end
            X = read_mesh(meshes.path(j));
            opt.debug = 0;
            opt.max_it = 200;
            opt.viz = 0;
            opt.split = 0;
            opt.thred = 2000;
            try
                disp([num2str(i), ' and ', num2str(j)])
                dist(j) = Identification(X, Y, opt, ['../../result/', num2str(i), num2str(j)]);
            catch
                dist(j) = 100000000;
            end
        end
        d(i, :) = dist;
    end

    writematrix(d, d_result_path);

    [~, I] = sort(d,2);
    pathmat = strings(meshes_num, meshes_num);

    for i = 1 : meshes_num
        pathmat(i ,:) = meshes.path(I(i,:));
    end

    writematrix(pathmat, p_result_path);
end
