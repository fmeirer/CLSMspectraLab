function [spectrum,nPixels,varargout] = getXANES_WizardSpectrum(obj,indices,filenames)
%getXANES_WizardSpectrum Gets spectrum from cluster results
%   
%   Usage:
%   [spectrum,nPixels] = getXANES_WizardSpectrum(obj,indices) outputs for
%   'spectrum' a cell matching the indices of the input images, which
%   contains a [nBins x nClusters] matrix. Analogously, nPixels is a cell 
%   array containing a [1 x nClusters] vector with the number of pixels per
%   cluster. The filename of the 'KmeansResults.mat' clustering result
%   output file of XANES Wizard is requested via an input dialog box. The
%   number of result files must match the number of indices.
%
%   [spectrum,nPixels] = getXANES_WizardSpectrum(...,filenames)
%   allows the user to specify the paths to the 'KmeansResults.mat'
%   file(s). Multiple files can be parsed as cell array.
%
%   [spectrum,nPixels,clusterMaps] = getXANES_WizardSpectrum(...)
%   clusterMaps is a matrix with the dimensions of the input image with the
%   cluster assignment per pixel. Pixels with a value of 0 are background.

if nargin < 2 || isempty(indices)
    indices = 1:obj.nInput;
end

if max(indices) > obj.nInput
    error('Not all specified indices do exist.')
end

if nargin < 3 || isempty(filenames)
   [file,path] = uigetfile('*.mat',...
        'Select one or more XANESWizard ''KmeansResults.mat'' output files', ...
        'MultiSelect', 'on');
   if file == 0
       spectrum = [];
       nPixels = [];
       return
   end
   filenames = fullfile(path,file);
end

if ischar(filenames) % is only one filename
    filenames = {filenames}; % save in cell
end

if numel(indices) ~= numel(filenames)
    error('The number of indices does not match the number of filenames.')
end

% preallocate
spectrum = cell(numel(indices),1);
nPixels = cell(numel(indices),1);
if nargout > 2
    clustermaps = cell(numel(indices),1);
end

% turn mask off
maskState = obj.maskFlag;
obj.maskFlag = false;

for ii = 1:numel(indices)
    index = indices(ii);
    % open files
    try
        out = open(filenames{ii});
    catch
        error('Could not open file: ''%s''',filenames{ii})
    end
    
    if ~isfield(out,'ClusterResult')
        error('Cluster results were not found in file: ''%s''',filenames{ii})
    end
    
    ClusterResult = out.ClusterResult;
    
    if obj.input(index).isVolume || obj.input(index).isTimeseries
        warning('Data set %i contains is a volume and/or timeseries. These dimensions are averaged.',indices(ii))
    end
    
    I = getImageProcessed(obj,index,'mean','x','y','c'); % in order to have the mask applied, we cannot directly request 'c'
    
    if size(ClusterResult{1,3},1) ~= size(I,1) || size(ClusterResult{1,3},2) ~= size(I,2)
        error('Dimensions input and cluster result do not match. Got X %i and %i, Y %i and %i, respectively.',size(ClusterResult{1,3},1),size(I,1),size(ClusterResult{1,3},2),size(I,2))
    end
    
    if nargout > 2 % store clustermaps
        clustermaps{ii} = ClusterResult{1,3};
    end
    
    nClusters = max(ClusterResult{1,3},[],'all');
    spectrum{ii} = nan(obj.input(index).getDim('c'),nClusters);
    nPixels{ii} = nan(1,nClusters);
    
    for jj = 1:nClusters
        thisCluster = ClusterResult{1,3} == jj;
        nPixels{ii}(jj) = sum(thisCluster,'all');
        % compute mean
        spectrum{ii}(:,jj) = sum(bsxfun(@times, cast(I,'double'), cast(thisCluster,'double')),[1 2])./nPixels{ii}(jj);
        
    end

end

% restore mask
obj.maskFlag = maskState;

% output 
if nargout > 2
    varargout{1} = clustermaps;
end

end

