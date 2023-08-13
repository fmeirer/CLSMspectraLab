classdef importXANESwizard < pixelClustering.abstract
    %FIRSTPCMASK Computes the mask
    %   Detailed explanation goes here
    
    properties(Constant)
        name = 'importXANESwizard';
    end
    
    properties
        loadFilepath = []; % filepath to the 'KmeansResults.mat' output by XANES Wizard
        nClusters = []; % number of clusters
    end
    
    properties(Hidden)
    end
    
    methods
        function obj = importXANESwizard(varargin)
            %FIRSTPCMASK Construct an instance of this class
            %
            % Usage:
            %   obj = firstpcMask(varargin) with varargin Name-Value pairs,
            %   e.g. if we set the property 'prop' to 5: 
            %   obj = firstpcMask('prop',5).
            if nargin > 0
                nIn = numel(varargin);
                if mod(nIn,2) ~= 0
                    error('Variables can only be parsed as Name-Value pairs')
                end
                for ii = 1:nIn/2
                    try
                        obj.(varargin{2*ii-1}) = varargin{2*ii};
                    catch
                        error(['The property ''' varargin{2*ii-1} ''' is not recognised.'])
                    end
                end
            end
        end
        
        function obj = compute(obj,~)
            %COMPUTE Loads the XANES Wizard data 'KmeansResults.mat'
            %    and stores the cluster map as obj.I
            
            
            if isempty(obj.loadFilepath)
               [file,path] = uigetfile('*.mat',...
                    'Select one XANESWizard ''KmeansResults.mat'' output files', ...
                    'MultiSelect', 'off');
               if file == 0
                   error('No file loaded.')
               end
               filename = fullfile(path,file);
            else
               filename = obj.loadFilepath;
            end
            
            % open files
            try
                out = open(filename);
            catch
                error('Could not open file: ''%s''',filename{ii})
            end
            
            if ~isfield(out,'ClusterResult')
                error('Cluster results were not found in file: ''%s''',filename{ii})
            end
            
            if min(out.ClusterResult{1,3},[],'all') == 0
                clustermap = out.ClusterResult{1,3}+1; % cluster 0 is background, but treated here as a unique cluster
            else
                clustermap = out.ClusterResult{1,3};
            end
            obj = obj.addImage(clustermap,'x','y'); % 3D is not yet supported
            obj.nClusters = max(obj.I,[],'all'); % starts counting at 0
        end

        function obj = setLoadFilenames(obj,filepath)
            % SETLOADFILENAME allows to set the filepath of the XANES Wizard data 'KmeansResults.mat'
            % If this function is not run, a file selection dialog box will
            % prompt for a file path when running compute.

            obj.loadFilepath = filepath;
        end

    end

end