classdef firstpcMask < imageMask.abstract
    %FIRSTPCMASK Computes the mask
    %   Detailed explanation goes here
    
    properties(Constant)
        name = 'First principal component mask';
    end
    
    properties
        medianFilterSize = [7 7 7]; % 1x3 matrix / size of the median filter
        thresholds; % threshold value of eigen image. Leave empty to determine value automatically
    end
    
    properties(Hidden)
        thresholdPerSlice; % threshold per slice; threshold of the first eigenimage, indexing: [z,number of thresholds].
        eigenI; % first eigen image used for thresholding
        nThresholds = 1;
    end
    
    methods
        function obj = firstpcMask(varargin)
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
        
        function obj = compute(obj,imageStack)
            %COMPUTE computes the domains with similar spectra and stores a mask
            %
            % The first principal component of the image/volume is computed and
            % thresholded using Otsu's method. The method can be extented to
            % multithesholding, but this has not been implemented yet.
            %
            % Code is not yet compatible with timeseries.
            %
            % Usage:
            % obj = compute(obj,imageStack) computes and stores the mask of
            % the imageStack object.
            
            if imageStack.isTimeseries
                error('Code cannot (yet) process timeseries.')
            end
            
            % hardcoded variables
            nThs = obj.nThresholds; % number of thresholds (only one now, check commented sections for multithreshold behavior)
            medfiltersize = obj.medianFilterSize; % get size of the median filter
            
            % extract full data volume: X x Y x spectrum x z-index
            Vr = imageStack.getReshapedImage('x','y','c','z');

            nX = imageStack.getDim('x');
            nY = imageStack.getDim('y');
            nZ = imageStack.getDim('z');
            nC = imageStack.getDim('c');

            % PCA
            if isempty(obj.eigenI)
                I = nan(nX,nY,nZ); % preallocation
                for ii = 1:nZ
                    % get z slice: x y ch
                    Vxyc = reshape(squeeze(Vr(:,:,:,ii)),nX*nY,nC);
                    % do a PCA:
                    [~, score, ~] = pca(double(Vxyc),'Centered','on','Economy','on');
                    % get first eigenimage
                    I(:,:,ii) = reshape(score(:,1),nX,nY);
                end
                
                % save eigenimage
                obj.eigenI = I;
            end
            
            if isempty(obj.thresholds)
                th = nan(nZ,nThs);
                for ii = 1:nZ
                    % get thresholds of eigenimage using Otsu's method
                    th(ii,:) = multithresh(obj.eigenI(:,:,ii),nThs);
                end
                obj.thresholdPerSlice = th;
            end

            if imageStack.isVolume
                if isempty(obj.thresholds)
                    % use average of last 10 values (check plot if OK)
                    if nZ>9
                        for ii=1:nThs
                            obj.thresholds = mean(th(end-10:end,ii));
                        end
                    else
                        for ii=1:nThs
                            obj.thresholds = mean(th(:,ii));
                        end
                    end
                end
                % and binarize the whole stack:
                Ibin = imbinarize(obj.eigenI,obj.thresholds(1,1));

                % remove salt and pepper noise:
                Ibin_mdf = medfilt3(Ibin,medfiltersize);
                % fill holes:
                Ibin_mdf_filled = imfill(Ibin_mdf,18,'holes'); 
                % store in object
                obj = obj.addImage(Ibin_mdf_filled,'x','y','z');
            else
                if isempty(obj.thresholds)
                    obj.thresholds = th;
                end
                % binarize image:
                Ibin = imbinarize(obj.eigenI,obj.thresholds);
                % remove salt and pepper noise:
                Ibin_mdf = medfilt2(Ibin,medfiltersize(1:2));
                % fill holes:
                Ibin_mdf_filled = imfill(Ibin_mdf,8,'holes');
                % store in object
                obj = obj.addImage(Ibin_mdf_filled,'x','y');
            end
            
            obj.unit = imageStack.unit;

        end
        
        function varargout = plotThresholdHistogram(obj,ha,id)
            %PLOTTHRESHOLDHISTOGRAM plot the histogram of the eigenimage
            %with the threshold value.
            %
            % Usage:
            %
            %   plotThresholdHistogram(obj,ha,id) plots the histogram of 
            %   the eigenimage with the threshold value in axis ha. Leave
            %   ha empty to use the current axis. id is an integer number
            %   of the data set to be plotted.
            %
            %   hps = obj.plotThresholdHistogram(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj.plotThresholdHistogram(...) returns also the handle to
            %   the axes handle of the plot.
            
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            if isempty(obj(id).eigenI) || isempty(obj(id).thresholds)
                error('Please compute mask first.')
            end
            
            nThs = obj(id).nThresholds;
            th = obj(id).thresholdPerSlice;
            
            nZ = obj(id).getDim('z');
            [nR,nC] = subplotDimensions(nZ);
            for ii = 1:nZ
                thisI = obj(id).eigenI(:,:,ii);
                subplot(nR,nC,ii);
                histogram(thisI(:))
                hold on
                for jj = 1:nThs
                    xline(th(ii,jj));
                end
                title('Eigenimage threshold')
                xlabel('Score')
                ylabel('Counts')
            end
            
            % Output
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
        end
        
        function varargout = plotThresholdPerZSlice(obj,ha,id)
            %PLOTTHRESHOLDPERZSLICE plots the threshold value per z-slice
            %
            % Usage:
            %
            %   plotThresholdPerZSlice(obj,ha,id) plots the threshold value 
            %   per z-slice in axis ha. Leave ha empty to use the current axis.
            %   id is an integer number of the data set to be plotted.
            %
            %   hps = obj.plotThresholdPerZSlice(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj.plotThresholdPerZSlice(...) returns also the handle to
            %   the axes handle of the plot.
            
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            if isempty(obj(id).thresholdPerSlice)
                error('Please compute mask first with auto thresholding.')
            end
            
            th = obj(id).thresholdPerSlice;
            
            plot(ha,th);
            title('threshold(s) for each slice');
            xlabel('Slice #')
            ylabel('Threshold')
            
            % Output
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
        end
        
        function varargout = plotEigenImage(obj,ha,id)
            %PLOTEIGENIMAGE plots the first eigen image used with threshold
            %
            % Usage:
            %
            %   plotEigenImage(obj,ha,id) plots first eigen image in axis 
            %   ha. Leave ha empty to use the current axis.
            %   id is an integer number of the data set to be plotted.
            %
            %   hps = obj.plotEigenImage(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj.plotEigenImage(...) returns also the handle to
            %   the axes handle of the plot.
            
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            if isempty(obj(id).eigenI)
                error('Please compute mask first.')
            end
            
            nZ = obj(id).getDim('z');
            [nR,nC] = subplotDimensions(nZ);
            for ii = 1:nZ
                thisI = obj(id).eigenI(:,:,ii);
                subplot(nR,nC,ii);
                imagesc(thisI)
                colorbar
                axis square
            end
            
            % Output
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end

    end
end

