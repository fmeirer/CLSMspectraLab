function export4XANES_Wizard_bins(obj,dirname,varargin)
%EXPORT4XANES_WIZARD_BINS Export function to TXM_XANES_Wizard
%   
%   Usage:
%   export4XANES_Wizard(obj,dirname,varargin) exports the images in obj for import
%   in TXM_XANES_Wizard in the folder dirname. Leave dirname empty for a
%   file selection dialog. In the folder dirname, new folders are made with
%   'bin_' + name of the tag of the image . Each bin is save as tiff image. The
%   images are saved without mask applied. The mask is saved separately in
%   mask.mat, if obj.maskFlag is true. Background correction is applied if
%   obj.bgCorrectionFlag = true.
%   The exported channels are taken from the row vectors with bin indices 
%   supplied in varargin. The last exported channel is the non-nor
%   

if nargin < 4
    error('Please provide at least one bin.')
end

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
    if isempty(C(ia))
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

doMask = false;
if obj.maskFlag == true
    obj.maskFlag = false;
    doMask = true;
end

for ii = 1:obj.nInput
    if any(~ismember(obj.input(ii).dimLabel,{'x','y','c'})) % check if is XYC image
        warning('Input %i is not an image with dimensions ''x'', ''y'', ''c''. Skipping this input.',ii)
        continue
    end
    if ~isa(obj.normalization(ii),'pixelNormalization.sum2one') % check if normalization is sum
        warning('Input %i does not have a sum2one normalization. Skipping this input.',ii)
        continue
    end
    % only does XYC
    dirname_sub = [dirname filesep 'bin_' tags{ii}];
    if ~exist(dirname_sub, 'dir')
        mkdir(dirname_sub)
    end
    
    % get channels with normalization applied
    I_r = obj.getImageProcessed(ii,'reshaped','c','x','y');
    nBins = numel(varargin);
    nChannels = obj.input(ii).getDim('c'); % assume that number of c has not bee modified during processing
    I = nan(obj.input(ii).getDim('x'),obj.input(ii).getDim('y'),nBins+1);
    for jj = 1:nBins
        binIdx = varargin{jj};
        if max(binIdx) > nChannels
            error('Supplied bin indices exceed number of channels %i in input %i.',nChannels,jj)
        end
        I(:,:,jj) = sum(I_r(binIdx,:,:),1);
    end
    
    % for export convert to uint16
    if isa(I_r,'double')
        I = uint16(I.*65535);
    else
        I = uint16(I);
    end
    
    % get channel intensity without normalization applied
    normalizationFlag = obj.normalizationFlag; % store normalization flag
    
    obj.normalizationFlag = false; % get intensity without normalization
    I_r_noNormalization = obj.getImageProcessed(ii,'reshaped','c','x','y');
    obj.normalizationFlag = normalizationFlag; % restore normalization flag
    
    I_noNormalization = squeeze(sum(I_r_noNormalization(double(string([varargin{:}])),:,:),1)); % sum all channels for intensity channel
    I_noNormalization = double(I_noNormalization); % normalize as double
    I_noNormalization = I_noNormalization./max(I_noNormalization,[],'all'); % normalize intensity
    I(:,:,end) = uint16(I_noNormalization.*65535); % convert back to uint16
    
    fn = cell(nBins+1,1);
    for jj = 1:nBins+1
        fn{jj} = sprintf('%.5i.tiff',jj);
        imwrite(I(:,:,jj),fullfile(dirname_sub,fn{jj}))
    end
    
    writematrix((1:size(I,3))',fullfile(dirname_sub,'energies.txt'))
    
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
