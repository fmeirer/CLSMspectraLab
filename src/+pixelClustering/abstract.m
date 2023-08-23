classdef (Abstract) abstract < imageStack 
    %CLUSTERING Abstract class to compute clustering
    
    properties (Abstract,Constant)
        name; % name of image mask computation method
    end

    properties (Abstract)
        minPixGroupSize
        diagonalPixConnected
    end
    
    methods (Abstract)
        obj = compute(obj,imageStack); % compute the mask
    end
    
    methods
        function plotClusterMap(obj,ha,id)
            % PLOTCLUSTERMAP plots the clusters in figure handle 'ha' of
            % image with index 'id'
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            clusterMap = getClusterMap(obj,id);
            imagesc(ha,clusterMap)
            axis equal
            uni = unique(clusterMap);
            n = (max(uni) - min(uni)) + 1;
            colormap(lines(n))
            colorbar();
            title('Unprocessed cluster map')
            %Cmap = linspecer(obj(id).nClusters+1);
        end

        function plotProcessedClusterMaps(obj,id)
            % PLOTPROCESSEDCLUSTERMAP plots the clusters in figure handle 'ha' of
            % image with index 'id'
            
            if nargin < 2
                id = 1;
            end
            
            processedClusterMap = getProcessedClusterMap(obj,id);
            figure
            for ii = 1:obj(id).nClusters
                subplot(1,obj(id).nClusters,ii)
                imagesc(processedClusterMap{ii})
                colormap('gray')
                title(['Processed map cluster ' num2str(ii)])
            end
        end

        function plotProcessedClusterMapsComparison(obj,id)
            % PLOTPROCESSEDCLUSTERMAPCOMPARISON plots the clusters in 
            % figure handle 'ha' of image with index 'id'. Top row is
            % non-processed, bottom row is processed.
            
            if nargin < 2
                id = 1;
            end
            
            processedClusterMap = getProcessedClusterMap(obj,id);
            clusterMap = getClusterMap(obj,id);
            figure
            for ii = 1:obj(id).nClusters
                subplot(2,obj(id).nClusters,ii)
                imagesc(clusterMap == ii)
                colormap('gray')
                title(['Non-processed map cluster ' num2str(ii)])
                subplot(2,obj(id).nClusters,obj(id).nClusters+ii)
                imagesc(processedClusterMap{ii})
                colormap('gray')
                title(['Processed map cluster ' num2str(ii)])
            end
        end

        function clusterMap = getClusterMap(obj,id)
            if nargin < 2 || isempty(id)
                id = 1;
            end
            if isempty(obj(id).I)
                error('No cluster map available for data %s.\nPlease compute/load clustering first.', num2str(id))
            end
            clusterMap = obj(id).I;
        end

        function processedClusterMap = getProcessedClusterMap(obj,id)
            clusterMap = getClusterMap(obj,id);
            processedClusterMap = cell(obj(id).nClusters,1);

            for ii = 1:obj(id).nClusters
                processedClusterMap{ii} = process(obj,id,clusterMap == ii);
            end

            function BW = process(obj,id,BW)
                if obj(id).diagonalPixConnected
                    conn = 8;
                else
                    conn = 6;
                end
                if obj(id).minPixGroupSize > 1
                    BW = bwareaopen(BW, obj(id).minPixGroupSize, conn);
                end
            end

        end
    end
    
end

