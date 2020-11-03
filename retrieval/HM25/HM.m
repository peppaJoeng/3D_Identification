function HM(HMdir, modelname)
%HM : generate mat recording the path of the file (only test directory)
%   此处显示详细说明

    cladir  = dir( HMdir );
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
            if((~isequal( phasedir( j ).name, 'test'))||...
                ~phasedir( j ).isdir)               % 如果不是目录则跳过
                continue;
            end
            % 必须为test目录下的off文件
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
    save([modelpath, [modelname,'.mat']], 'meshes');
end

