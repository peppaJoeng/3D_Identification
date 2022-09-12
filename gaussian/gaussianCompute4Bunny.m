function gaussianCompute(mu, std)
    addpath(genpath('../core'));
    addpath(genpath('../utils'));
    addpath(genpath('../thirdparty/CPD2/core'));
    addpath(genpath('../thirdparty/inexact_alm_rpca'));
    addpath(genpath('../mex'));

    % N > M
    source_path = '../model/frequent/bunny.off';
    

    result_dir = 'bunny';
    if ~exist(result_dir,'dir')
        mkdir(result_dir);
    end

    diary([result_dir,'/bunny_', num2str(std), '.log']);
    diary on;

    opt.max_it = 200;
    opt.debug = 0;
    opt.viz = 0; 
    opt.segment = 0;
    opt.metric = "ALL";

    X = read_mesh(source_path);
    X = downsample(X, 0.2);
    
    noise1=normrnd(str2double(mu),str2double(std),size(X,1),1);    
    disp(noise1);
    noise2=normrnd(str2double(mu),str2double(std),size(X,1),1);    
    disp(noise2);
    noise3=normrnd(str2double(mu),str2double(std),size(X,1),1);    
    disp(noise3);
    
    Y=[X(:,1)+noise1,X(:,2)+noise2,X(:,3)+noise3];

    disp(size(X));
    disp(size(Y));
    path = [result_dir,'/bunny_', num2str(std), '.mat'];
    opt.savename = path;
    Identification(X, Y, opt, result_dir);

    diary off;
end