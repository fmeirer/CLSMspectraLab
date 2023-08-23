classdef similarityFit < pixelClustering.abstract
    %SIMILARITYFIT Computes the clustering by assigning each pixel to 
    %   Detailed explanation goes here
    
    properties(Constant)
        name = 'similarityFit';
    end
    
    properties
        loadFilepath = []; % filepath to the reference spectra (should be csv file with columns with spectra)
        nClusters = []; % number of clusters
        minPixGroupSize = 2; % minimum number of pixels in group to be considered a particle for in the processed cluster map
        diagonalPixConnected = true; % if true, connectivity of 8; otherwise, 6
        thresholdBg = 0.07;
        belowThresholdIdx = [1 4];
    end
    
    properties(Hidden)
    end
    
    methods
        function obj = similarityFit(varargin)
            %SIMILARITYFIT Construct an instance of this class
            %
            % Usage:
            %   obj = similarityFit(varargin) with varargin Name-Value pairs,
            %   e.g. if we set the property 'prop' to 5: 
            %   obj = similarityFit('prop',5).
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
        
        function obj = compute(obj,iS)
            %COMPUTE Loads the reference files and
            %    and stores the cluster map as obj.I
            
            if isempty(obj.loadFilepath)
                [file,path] = uigetfile(...
                     {'*.csv;*.dat;*.txt',...
                     'Delimited text files (*.csv;*.dat;*.txt)';...
                     '*.xls,*.xlsb,*.xlsm,*.xlsx,*.xltm,*.xltx,*.ods',...
                     'Spreadsheet files (*.xls,*.xlsb,*.xlsm,*.xlsx,*.xltm,*.xltx,*.ods)'},...
                     'Select one file with reference spectra', ...
                     'MultiSelect', 'off');
               if file == 0
                   error('No file loaded.')
               end
               filename = fullfile(path,file);
            else
               filename = obj.loadFilepath;
            end
            
            % open files
            refSpectra = readmatrix(filename); % provide address of CSV file for reference spectra with order as specified in template CSV file
            refSpectraNorm = normalize(refSpectra,1,"range"); %[bg,clay,USY,Al,Si]

            obj.nClusters = size(refSpectraNorm,2);

            if max(obj.belowThresholdIdx) > obj.nClusters
                error('Invalid spectrum index in "belowThresholdIdx"')
            end

            if any(ismember(iS.dimLabel,{'z','t'}))
                error('Dimensions ''z'' and ''t'' are not (yet) supported');
            end
            
            I = iS.getReducedImage('reshaped','c','x','y');
            I = double(I); % convert to double to compare with reference spectra
            
            Isize = iS.dimSize;
            
            Inorm = normalize(I,1,"range");
            Imean_norm = rescale(squeeze(sum(I))); % same as rescale but now over all dimensions

            residuals = nan(size(I,2),size(I,3),obj.nClusters); % default is residual of 1 per pixel, which is max possible for norm data
            
            aboveThresholdMask = Imean_norm > obj.thresholdBg;
            %aboveThresholdMask = ~logical(imsegkmeans(single(Imean_norm),2)-1);
            for ii = 1:obj.nClusters
                if any(ii == obj.belowThresholdIdx)
                    tmp = squeeze(sqrt(sum((Inorm-refSpectraNorm(:,ii)).^2)));
                    tmp(aboveThresholdMask) = Isize(1);
                    residuals(:,:,ii) = tmp;
                else
                    tmp = squeeze(sqrt(sum((Inorm-refSpectraNorm(:,ii)).^2)));
                    tmp(~aboveThresholdMask) = Isize(1);
                    residuals(:,:,ii) = tmp;
                end
            end

            [~,clustermap] = min(residuals,[],3);
            obj = obj.addImage(clustermap,'x','y'); % 3D is not yet supported

        end

        function obj = setLoadFilenames(obj,filepath)
            % SETLOADFILENAME allows to set the filepath of the reference
            % spectra.
            % If this function is not run, a file selection dialog box will
            % prompt for a file path when running compute.

            [obj.loadFilepath] = deal(filepath);
        end

    end

end