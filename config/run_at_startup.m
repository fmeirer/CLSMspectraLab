% RUN AT START-UP script for DiffusionLab

% add all paths (folders named private are automatically excluded
folder = fileparts(which(mfilename)); 
addpath(genpath(fullfile(folder,['..' filesep])));

javaaddpath('../bin/external/bfmatlab/bioformats_package.jar')