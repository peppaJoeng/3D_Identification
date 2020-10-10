function setup()

psave=pwd; p = mfilename('fullpath'); [pathstr, ~, ~] = fileparts(p);
disp('Compiling Identification Mex functions...');
disp('If this is the first time you use mex, it will ask you to choose a compiler.');
disp('Just choose the matlab default one (usually option #1).');
cd (pathstr); cd mex;

%% Compile the downsampling function and calculate the posterior probability function
try
    mex postProbablity.cpp;
    mex HGMM.cpp;
    disp('Compilation of mex functions is SUCCESSFUL.');
catch
   disp('Compilation of mex functions failed. Try to run mex -setup to adjust your compiler.');
end

disp('Compiling CPD Mex functions...');

%% Compile the mex file that CPD use
cd (pathstr); cd thirdparty; cd CPD2; cd core; 
cpd_make

%% Create directory to store resukt 
if ~exist('./result','dir')
    mkdir('./result')
end

cd (psave);
end