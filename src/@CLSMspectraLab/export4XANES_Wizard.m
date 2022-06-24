function export4XANES_Wizard(obj,dirname,varargin)
%EXPORT4XANES_WIZARD Export function to TXM_XANES_Wizard
%   
%   Usage:
%   export4XANES_Wizard(obj,dirname) exports the images in obj for import
%   in TXM_XANES_Wizard in the folder dirname. Leave dirname empty for a
%   file selection dialog. In the folder dirname, new folders are made with
%   name of the tag of the image. Each bin is save as tiff image. The
%   images are saved without mask applied. The mask is saved separately in
%   mask.mat, if obj.maskFlag is true. Background correction is applied if
%   obj.bgCorrectionFlag = true.
%
%   export4XANES_Wizard(obj,dirname,wavelengths) with wavelengths a vector
%   with the wavelength of each bin saves these value in a energies.txt. If
%   wavelengths is a cell with vectors, each cell element is saved for each
%   image.
%   

if nargin < 2 || isempty(dirname)
    dirname = uigetdir();
    if dirname == 0
        return
    end
end

fprintf('Writing to folder %s.\n',dirname)

tags = {obj.input.tag};
[C,ia,ic] = unique(tags);
for ii = 1:numel(ia)
    idx = find(ic==ii);
    if isempty(C{ia})
        for jj = 1:numel(idx)
            tags{idx(jj)} = 'no_name';
        end
    end
    if numel(idx) < 2
        continue
    end
    for jj = 1:numel(idx)
        tags{idx(jj)} = [tags{idx(jj)} '_' num2str(jj)];
    end
end

if nargin > 2
    if isnumeric(varargin{1})
        wavelengths_bin_cell = cell(obj.nInput,1);
        wavelengths_bin_cell(:) = varargin{1};
    else
        wavelengths_bin_cell = varargin{1};
    end
end

doMask = false;
if obj.maskFlag == true
    obj.maskFlag = false;
    doMask = true;
end

for ii = 1:obj.nInput
    % only does XYC
    dirname_sub = [dirname filesep tags{ii}];
    if ~exist(dirname_sub, 'dir')
        mkdir(dirname_sub)
    end
%     I = obj.input(ii).getReshapedImage('x','y','c');
    I = obj.getImageProcessed(ii,'reshaped','x','y','c');
    fn = cell(obj.input(ii).getDim('c'),1);
    for jj = 1:obj.input(ii).getDim('c')
        fn{jj} = sprintf('%.5i.tiff',jj);
        imwrite(I(:,:,jj),fullfile(dirname_sub,fn{jj}))
    end
    
    if nargin > 2
        wavelengths_bin = wavelengths_bin_cell{ii};
        writematrix(wavelengths_bin(:),fullfile(dirname_sub,'energies.txt'))
    end
    
    if doMask
        % mask
        imwrite(obj.mask(ii).I,fullfile(dirname_sub,'mask.tif'))
    end
    
end

fprintf('Done.\n')

% write XANESmap
% 
% if doMask
%     % mask
%     filter_edgeJump = double(obj.mask(ii).I);
%     if nargin > 2
%         energies = wavelengths_bin(:);
%     else
%         energies = (1:obj.input(ii).getDim('c'))';
%     end
%     % rest is dummy variables
%     filler = ones(size(filter_edgeJump));
%     FirstPreEdgePt = 1;
%     LastPreEdgePt = 2;
%     FirstPostEdgePt = 3;
%     LastPostEdgePt = 4;
%     avgPreEdge = filler;
%     avgPostEdge = filler;
%     EdgeJump = filler;
%     fnames = {[dirname_sub filesep]}; % add filesep to prevent error in XANES Wizard
%     image_list = fn;
%     StdDev_threshold = 1;
%     stdPostEdge = 1;
%     stdPreEdge = 1;
%     experimental_error = 1;
%     save(fullfile(dirname_sub,'XANESmaps.mat'),'filter_edgeJump','energies','FirstPreEdgePt',...
%         'LastPreEdgePt','FirstPostEdgePt','LastPostEdgePt','avgPreEdge','avgPostEdge',...
%         'EdgeJump','fnames','image_list','StdDev_threshold','stdPostEdge','stdPreEdge',...
%         'experimental_error')
% end
