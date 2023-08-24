function getClusterStats(obj,savepath,indices)
% GETCLUSTERSTATS saves some statistical properties of the cluster map
%   Usage:
%   obj.getClusterStats()L Saves in .xlsx file selected by user via input
%       dialog for all data.
%
%   obj.getClusterStats(savepath): Saves in path specified in 'savepath'
%
%   obj.getClusterStats(savepath,indicies): only saves the data with
%   indices in 

if nargin < 2 || isempty(savepath)
    [file,path] = uiputfile('*.xlsx',...
        'Save file name');
    if file == 0
        return
    end
    savepath = fullfile(path,file);
end

if nargin < 3 || isempty(indices)
    indices = 1:obj.nInput;
end

if max(indices) > obj.nInput || min(indices) < 1
    error('Please provide valid indices')
end

T_mean = [];
for ii = indices
    if isempty(obj.clustering(ii).I)
        continue
    end
    [pixelSize,unitName] = obj.input(1).unit.getUnitFactor('x');
    processedClusterMap = getProcessedClusterMap(obj.clustering,ii);
    for jj = 1:obj.clustering(ii).nClusters
        thisMap = processedClusterMap{jj};
        T = getStats(thisMap,pixelSize,unitName);
        writetable(T,savepath,"WriteMode","overwritesheet","sheet",sprintf('%s_Cluster%s',num2str(ii),num2str(jj)))
        T_mean = [T_mean; varfun(@mean, T, 'InputVariables', @isnumeric)];
    end
end
if ~isempty(T_mean)
    writetable(T_mean,savepath,"WriteMode","overwritesheet","sheet",sprintf('Mean'));
end

end

function T = getStats(BW,pixelSize,unitName)
    T1 = regionprops("table",BW,"Area","FilledArea","Eccentricity","Centroid");
    NN = nan(size(T1,1),1);
    for ii = 1:numel(NN)
        [~,NN(ii)] = knnsearch(T1.Centroid(ii,1),T1.Centroid(ii,2),'K',1,'Distance','euclidean');
    end
    T2 = table(NN);
    T = [T1 T2];

    if ~isempty(pixelSize)
        T_NN_unit = table(T.Area .* pixelSize^2, T.FilledArea .* pixelSize^2, T.NN .* pixelSize,'VariableNames',{['Area_' unitName '2'],['FilledArea_' unitName '2'],['NN_' unitName]});
    end

    T = [T T_NN_unit];
end