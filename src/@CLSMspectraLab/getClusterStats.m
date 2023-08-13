function getClusterStats(obj,savepath)

if nargin < 2 || isempty(savepath)
    [file,path] = uiputfile('*.xlsx',...
        'Save file name');
    if file == 0
        return
    end
    savepath = fullfile(path,file);
end

N = getPreallocSize(obj);
T = cell(N,1);

kk = 1;
for ii = 1:obj.nInput
    for jj = 1:obj.clustering(ii).nClusters
        thisMap = obj.clustering(ii).I == jj;
        T{kk} = getStats(thisMap);
        writetable(T{kk},savepath,"WriteMode","overwritesheet","sheet",sprintf('%s_Cluster%s',num2str(ii),num2str(jj)))
        kk = kk + 1;
    end
end

end

function N = getPreallocSize(obj)

N = 0;
for ii = 1:obj.nInput
    N = N + obj.clustering(ii).nClusters;
end

end

function T = getStats(BW)

    T1 = regionprops("table",BW,"Area","FilledArea","Eccentricity","Centroid");
    T1(T1.Area <= sizeThreshold,:) = [];
    NN = nan(size(T1,1),1);
    for ii = 1:numel(NN)
        [~,NN(ii)] = knnsearch(T1.Centroid(ii,1),T1.Centroid(ii,2),'K',1,'Distance','euclidean');
    end
    T2 = table(NN);
    T = [T1 T2];
%     Area = mean(T.Area);
%     FilledArea = mean(T.FilledArea);
%     Eccentricity = mean(T.Eccentricity);
%     NN = mean(NN);
end