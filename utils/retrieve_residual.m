function retrieve_residual(queryIndex, dbIndex, computeIndex)
    addpath(genpath('../utils'));
    addpath(genpath('../core'));
    addpath(genpath('../thirdparty/CPD2/core'));
    addpath(genpath('../thirdparty/inexact_alm_rpca'));
    addpath(genpath('../mex'));
    
    metric="CORR";
    dpath = strcat('dis_residual_', num2str(queryIndex), '_', num2str(dbIndex), '_', num2str(computeIndex), '_',metric, '.csv');
    ppath = strcat('path_residual_', num2str(queryIndex), '_', num2str(dbIndex), '_', num2str(computeIndex), '_', metric, '.csv');
    Ipath = strcat('I_residual_', num2str(queryIndex), '_', num2str(dbIndex), '_', num2str(computeIndex), '_', metric, '.csv');
    tmppath = strcat('Tmp_residual_', num2str(queryIndex), '_', num2str(dbIndex), '_', num2str(computeIndex),  '_', metric, '.csv');
    disp(dpath)
    disp(ppath)
    disp(Ipath)
    disp(tmppath)
    
    [Q, D]=prepareData(queryIndex,dbIndex);
    prefix='path/to/benchmark_datasets';
    disp(length(Q))
   
   
    opt.debug = 0;
    opt.max_it = 200;
    opt.viz = 0;
    opt.segment = 0;
    opt.metric = metric;
    
    disp(['computeIndex : ',num2str(computeIndex)]);
    switch num2str(computeIndex)
        case '1'
            disp("Number: 1-7");
            Q=Q(1:7,:);
        case '2'
            disp("Number: 8-14");
            Q=Q(8:14,:);
        case '3'
            disp("Number: 15-21");
            Q=Q(15:21,:);
        case '4'
            disp("Number: 22-28");
            Q=Q(22:28,:);
        case '5'
            disp("Number: 29-35");
            Q=Q(29:35,:);
        case '6'
            disp("Number: 36-42");
            Q=Q(36:42,:);
        case '7'
            disp("Number: 43-49");
            Q=Q(43:49,:);
        case '8'
            disp("Number: 50-56");
            Q=Q(50:56,:);
        case '9'
            disp("Number: 57-63");
            Q=Q(57:63,:);
        case '10'
            disp("Number: 64-70");
            Q=Q(64:70,:);
        case '11'
            disp("Number: 71-74");
            Q=Q(71:75,:);
        case '12'
            disp("test");
            Q=Q(1:2,:);
            D=D(1:3,:);
    end
    
    queryNum= length(Q);
    dbNum = length(D);
    
    disp(["Begin queryindex: ", num2str(queryIndex), " dbindex : ", num2str(dbIndex)]);
    disp(["queryNum : ", num2str(queryNum)]);
    disp(["dbNum : ", num2str(dbNum)]);
    d = zeros(queryNum, dbNum);
    for i = 1:queryNum
        query=Q(i);
        queryPath = fullfile(prefix, query);
        queryPath=bin2off(queryPath);
        disp(["Now begin ", queryPath])
        Y = read_mesh(queryPath);
       
        dist = zeros(1, dbNum);
        for j = 1:dbNum
            disp([num2str(i), ' and ', num2str(j)]);
            db=D(j);
            dbPath=fullfile(prefix, db);
            dbPath = bin2off(dbPath);
            X = read_mesh(dbPath);
            try
                dist(j) = Identification(X, Y, opt, '');
            catch
                disp("error");
                dist(j) = 0;
            end    
        end
        d(i, :) = dist;
        writematrix(d, tmppath);
    end
    
    writematrix(d, dpath);
    
    [~, I] = sort(d, 2, 'descend');
    writematrix(I, Ipath);
    pathmat = strings(queryNum, dbNum + 1);

    for i = 1 : queryNum
        query=Q(i);
        queryPath=bin2off(query);
        pathmat(i, 1) = queryPath;
        pathmat(i, 2 : dbNum + 1) = D(I(i, 1 : dbNum));
    end
    writematrix(pathmat, ppath); 
end

function res = bin2off(path)
    path=strsplit(path,".");
    path=path(1);
    res = strcat(path, '.off');
    res=char(res);
end

function [Qlist,Dlist]=prepareData(queryIndex,dbIndex)
    prefix='/Users/kiki/Projects/repos/PCAN/';
    databasepath=[prefix, 'residential_evaluation_database_file.csv'];
    disp(['Read Done! ', databasepath]);
    D = importdata(databasepath);
    Dlist=[];
    for i = 1:length(D)
        d=string(D(i));
        if belong(d, dbIndex)
            Dlist=[Dlist;d];
        end
    end
    disp(length(Dlist))
    
    querypath=[prefix, 'residential_evaluation_query_file.csv'];
    disp(['Read Done! ', querypath]);
    Q = importdata(querypath);
    Qlist=[];
    for i = 1:length(Q)
        q=string(Q(i));
        if belong(q, queryIndex)
            Qlist=[Qlist;q];
        end
    end
    
end

function choosed=belong(path, index)
    choosed=0;
    folder=strsplit(path,'/');
    folder=folder(2);
    if strcmp(folder,['residential_run',num2str(index)])
        choosed=1;
    end
end
