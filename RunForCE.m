%% Run ForCE Model one-line cross-shore / longshore model
% Mark Davidson 13/6/2024 


%% Initialise Matlab
    clear  ; close all; clc

%% Record base directory
     baseDir=which('ForCE.m');
     baseDir=baseDir(1:end-7);
     cd (baseDir)
     addpath([pwd '/linear'])

%% Get setup file(s)
    disp('Please load the model setup file found in the project directory')
    [file,path] = uigetfile('*.m','MultiSelect','on'); 
    addpath(pwd); addpath(path);
    basePath=path(1:end-10);
    cd(basePath)

%% Find number of setup files
if ischar(file)
    Nfiles=1;
else
    Nfiles=length(file);
end

disp(['No. files to analyse = ',num2str(Nfiles)])


%% Run force for each setup file
for ifile=1:Nfiles
    % Run setup file
    if Nfiles==1
        filename.setup=fullfile(path,file);
    else
        filename.setup=fullfile(path,file{ifile});
    end
        run(filename.setup)    
        disp(['Loading file: ',filename.setup])

    % Run ForCE model
    run('ForCE')

    % Clear key variables
    clear model filename constant
end

disp('All Done!')
cd (baseDir)