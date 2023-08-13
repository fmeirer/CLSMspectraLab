classdef (Abstract) abstract < imageStack 
    %CLUSTERING Abstract class to compute clustering
    
    properties (Abstract,Constant)
        name; % name of image mask computation method
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
            
            if isempty(obj(id).I)
                error('Please compute/load clustering first.')
            end

            imagesc(ha,obj(id).I)
            axis equal
            n = obj(id).nClusters;
            colormap(lines(n))
            colorbar();
            %Cmap = linspecer(obj(id).nClusters+1);
            

        end
    end
    
end

