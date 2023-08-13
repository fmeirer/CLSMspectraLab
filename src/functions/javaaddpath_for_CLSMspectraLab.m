function javaaddpath_for_CLSMspectraLab()
% adds java path bfmatlab

[path,~,~] = fileparts(mfilename("fullpath"));
javaaddpath([path '/../../bin/external/bfmatlab/bioformats_package.jar'])

end

