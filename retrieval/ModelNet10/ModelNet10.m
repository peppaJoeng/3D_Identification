clear all; close all; clc;

model10dir = '../../data/ModelNet10/';
cladir  = dir( model10dir );
modelpath = '../../model/';
total = 0;
meshes.cla = [];
% cladir : classfication folder
for i = 1 : length(cladir)
    if( isequal( cladir( i ).name, '.' )||...
        isequal( cladir( i ).name, '..')||...
        ~cladir( i ).isdir)               % 如果不是目录则跳过
        continue;
    end
    meshes.claname{i} = cladir(i).name;
    % cladirpath : absolute path
    cladirpath = fullfile( cladir(i).folder, '/', cladir(i).name );
    % phasedir : train or test  
    phasedir = dir( cladirpath );
    for j = 1 : length(phasedir)
        if((~isequal( phasedir( j ).name, 'train' ) && ...
            ~isequal( phasedir( j ).name, 'test'))||...
            ~phasedir( j ).isdir)               % 如果不是目录则跳过
            continue;
        end
        phasedirpath = fullfile( phasedir(j).folder, '/', phasedir( j ).name );
        obj = dir( [phasedirpath , '/*.off'] );
        for k = 1 : length(obj)
            objpath = fullfile( obj(k).folder, '/', obj (k ).name );
            % disp(string(objpath));
            total = total + 1;
            meshes.path(total) = string(objpath);
            meshes.cla(end+1) = i;
        end
    end
end
save([modelpath, 'modelnet10.mat'], 'meshes');
