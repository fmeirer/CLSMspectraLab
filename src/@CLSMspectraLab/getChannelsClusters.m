function spectrum = getChannelsClusters(obj,indices)
%getSpectrumClusters Gets spectrum from cluster results
%   
%   Usage:
%   spectrum = getChannelsClusters(obj,indices) outputs for
%   'spectrum' a cell matching the indices of the input images, which
%   contains a [nBins x nClusters] matrix.
%   The filename of the 'KmeansResults.mat' clustering result
%   output file of XANES Wizard is requested via an input dialog box. The
%   number of result files must match the number of indices.

if max(indices) > obj.nInput || min(indices) < 1
    error('Please provide valid indices')
end

% preallocate
spectrum = cell(numel(indices),1);
nPixels = cell(numel(indices),1);

% turn mask off
maskState = obj.maskFlag;
obj.maskFlag = false;

for ii = 1:numel(indices)
    index = indices(ii);
    
    ClusterResult = obj.clustering.getClusterMap(index);
    
    if obj.input(index).isVolume || obj.input(index).isTimeseries
        warning('Data set %i contains is a volume and/or timeseries. These dimensions are averaged.',indices(ii))
    end
    
    I = getImageProcessed(obj,index,'mean','x','y','c'); % in order to have the mask applied, we cannot directly request 'c'
    
    if size(ClusterResult,1) ~= size(I,1) || size(ClusterResult,2) ~= size(I,2)
        error('Dimensions input and cluster result do not match. Got X %i and %i, Y %i and %i, respectively.',size(CClusterResult,1),size(I,1),size(CClusterResult,2),size(I,2))
    end
    
    
    spectrum{ii} = nan(obj.input(index).getDim('c'),obj.clustering(ii).nClusters);
    nPixels{ii} = nan(1,obj.clustering(ii).nClusters);
    
    for jj = 1:obj.clustering(ii).nClusters
        thisCluster = ClusterResult == jj;
        nPixels{ii}(jj) = sum(thisCluster,'all');
        % compute mean
        spectrum{ii}(:,jj) = sum(bsxfun(@times, cast(I,'double'), cast(thisCluster,'double')),[1 2])./nPixels{ii}(jj);
    end

end

% restore mask
obj.maskFlag = maskState;

end

