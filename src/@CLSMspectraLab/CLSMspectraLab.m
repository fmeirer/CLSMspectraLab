classdef CLSMspectraLab
    %CLSMSPECTRALAB is a class for analysis and plotting spectral CLSM images
    %   
    %   Order of excution processing steps:
    %       (1) background correction
    %       (2) Masking
    %       (3) Dimensionality reduction (to be added)
    %       (4) Clustering (to be added)
    
    properties
        input; % imagestack
        nInput; % number of input files: number of imageStack objects
        ROItable; % table with polygon ROIs
        bgCorrection = bgCorrection.linear; % background correction from bgCorrection domain
        bgCorrectionFlag = true; % if true, the background correction is applied
        mask = imageMask.firstpcMask; % mask from imageMask domain
        maskFlag = true; % if true, the mask is applied
        normalization = pixelNormalization.sum2one; % normalization from normalization domain for clustering
        normalizationFlag = true; % if true, the normalization is applied
        clustering = pixelClustering.importXANESwizard; % clustering of pixels into classes
    end

    properties(SetAccess = private,Hidden)
        binFactor = 1;
    end
    
    methods
        function obj = CLSMspectraLab(filepath,binFactor)
            %CLSMSPECTRALAB Construct an instance of this class
            %   obj = CLSMspectraLab(filepath): accepts 
            %   a string with the file path to the image file or a cell
            %   with the file path to multiple files. If empty a file 
            %   selector is opened. Returns an instance of the  
            %   CLSMspectraLab class.
            %
            %   obj = CLSMspectraLab(iStack): accepts an imageStack object
            %   or array of objects. Returns an instance of the  
            %   CLSMspectraLab class.
            %
            %   obj = CLSMspectraLab(...,binFactor): allows the
            %   storage of the binFactor, which can be used by other
            %   functions. Default values is 1 (no binning).
            
            if nargin < 1
                filepath = [];
            end

            if nargin > 1
                obj.binFactor = binFactor;
            end
            
            % write to object
            if isa(filepath,'imageStack')
                obj.input = filepath(:);
                obj.nInput = numel(filepath);
            else
                [obj.input,obj.nInput] = imageStack.import(filepath);
            end
            
        end
        
        function obj = computeBgCorrection(obj,id)
            %COMPUTEBGCORRECTION computes the background correction for the
            % input image
            %
            %   Usage:
            %   obj = computeBgCorrection(obj) computes the background 
            %   correction following the 1st bgCorrection object in 
            %   obj.bgCorrection for all different input images.
            %
            %   obj = computeBgCorrection(obj,id) computes only for index
            %   'id', which must be a scalar or vector of indices.

            if nargin < 2 || isempty(id)
                id = 1:obj.nInput;
            else
                if max(id) > obj.nInput
                    error('Largest value id exceeds number of images, got %i (max is %i).',max(id),obj.nInput)
                end
            end

            if numel(obj.bgCorrection) ~= numel(obj.input)
                obj = initBgCorrection(obj);
            end

            [iStack,mask_,name] = getProcessedImageStack(obj,id,'bgCorrection');
            
            % note that bgCorrection also needs mask_ for linear
            % bgCorrection
            if isempty(mask_)
                for ii = 1:numel(id)
                    obj.bgCorrection(id(ii)) = obj.bgCorrection(id(ii)).compute(iStack(ii));
                end
            else
                for ii = 1:numel(id)
                    obj.bgCorrection(id(ii)) = obj.bgCorrection(id(ii)).compute(iStack(ii),mask_(ii));
                end
            end
            fprintf('Computed background correction using %s image.\n',name)
            
%             for ii = 1:numel(obj.bgCorrection)
%                 obj.bgCorrection(ii) = obj.bgCorrection(ii).compute(obj.input(ii));
%             end
            
        end
        
        function obj = initBgCorrection(obj)
            %INITBGCORRECTION repeats the first bgCorrection for the number of input images
            %
            %   Usage:
            %   obj = initBgCorrection(obj) makes a vector of N copies of 
            %   the first bgCorrection object in obj.bgCorrection, with N the 
            %   number of input objects at obj.input.
            
            if isempty(obj.bgCorrection)
                error('No background correction method selected. Please assign bgCorrection object first.')
            end
            
            em = obj.bgCorrection.empty(0,numel(obj.input));
            em(1,:) = obj.bgCorrection(1);
            obj.bgCorrection = em;
            
        end
        
        function obj = computeMask(obj,id)
            %COMPUTEMASK computes the mask for the input image
            %
            %   Usage:
            %   obj = computeMask(obj) computes the mask following the
            %   1st inputMask object in obj.mask for all different input
            %   images.
            %
            %   obj = computeMask(obj,id) computes only for index
            %   'id', which must be a scalar or vector of indices.
            
            if nargin < 2 || isempty(id)
                id = 1:obj.nInput;
            else
                if max(id) > obj.nInput
                    error('Largest value id exceeds number of images, got %i (max is %i).',max(id),obj.nInput)
                end
            end

            if numel(obj.mask) ~= numel(obj.input)
                obj = initMask(obj);
            end

            [iStack,~,name] = getProcessedImageStack(obj,id);
            
            for ii = 1:numel(id)
                obj.mask(id(ii)) = obj.mask(id(ii)).compute(iStack(ii));
            end
            fprintf('Computed mask using %s image.\n',name)
            
%             if ~obj.bgCorrectionFlag || isempty(obj.bgCorrection(1).I) % without background correction
%                 for ii = 1:numel(obj.mask)
%                     obj.mask(ii) = obj.mask(ii).compute(obj.input(ii));
%                 end
%                 fprintf('Computed mask without background correction.\n')
%             else
%                 for ii = 1:numel(obj.mask)
%                     obj.mask(ii) = obj.mask(ii).compute(obj.bgCorrection(ii));
%                 end
%                 fprintf('Computed mask with background correction.\n')
%             end
            
        end
        
        function obj = initMask(obj)
            %INITMASK repeats the first mask for the number of input images
            %
            %   Usage:
            %   obj = initMask(obj) makes a vector of N copies of the first
            %   mask object in obj.mask, with N the number of input
            %   objects at obj.input.
            
            if isempty(obj.mask)
                error('No mask selected. Please assign mask first.')
            end
            
            em = obj.mask.empty(0,numel(obj.input));
            em(1,:) = obj.mask(1);
            obj.mask = em;
            
        end
        
        function obj = invertMask(obj,idx)
            %INVERTMASK inverts the values for an image containing 0 and 1
            %
            %   Usage:
            %   obj = invertMask(obj,idx) converts 0 -> 1 and vice versa in
            %   mask. idx is the number of the masks that have to be
            %   inverted, leave empty for all masks
            
            if nargin < 2 || isempty(idx)
                idx = 1:numel(obj.mask);
            end

            if max(idx) > numel(obj.mask)
                error('Idx is larger than the number of image masks.')
            end

            for ii = idx
                obj.mask(ii) = obj.mask(ii).invertImage;
            end
            
        end

        function obj = computeNormalization(obj,id)
            %COMPUTENORMALIZATION computes the normalization for the input image
            %
            %   Usage:
            %   obj = computeNormalization(obj) computes the normalization 
            %   following the 1st normalization object in obj.normalization
            %   for all different input images.
            %
            %   obj = computeNormalization(obj,id) computes only for index
            %   'id', which must be a scalar or vector of indices.
            
            if nargin < 2 || isempty(id)
                id = 1:obj.nInput;
            else
                if max(id) > obj.nInput
                    error('Largest value id exceeds number of images, got %i (max is %i).',max(id),obj.nInput)
                end
            end

            if numel(obj.normalization) ~= numel(obj.input)
                obj = initNormalization(obj);
            end

            [iStack,~,name] = getProcessedImageStack(obj,id,'normalization');
            
            
            for ii = 1:numel(iStack)
                obj.normalization(id(ii)) = obj.normalization(id(ii)).compute(iStack(ii));
            end
            fprintf('Computed normalization using %s image.\n',name)
            
        end

        function obj = initNormalization(obj)
            %INITNORMALIZATION repeats the first normalization for the number of input images
            %
            %   Usage:
            %   obj = initNormalization(obj) makes a vector of N copies of 
            %   the first normalization object in obj.normalization, with N the 
            %   number of input objects at obj.input.
            
            if isempty(obj.normalization)
                error('No background correction method selected. Please assign bgCorrection object first.')
            end
            
            em = obj.normalization.empty(0,numel(obj.input));
            em(1,:) = obj.normalization(1);
            obj.normalization = em;
            
        end

        function obj = computeClustering(obj,id)
            %COMPUTECLUSTERING computes the clustering for the input image
            %
            %   Usage:
            %   obj = computeClustering(obj) computes the clustering following the
            %   1st pixelClustering object in obj.clustering for all different input
            %   images.
            %
            %   obj = computeClustering(obj,id) computes only for index
            %   'id', which must be a scalar or vector of indices.
            
            if nargin < 2 || isempty(id)
                id = 1:obj.nInput;
            else
                if max(id) > obj.nInput
                    error('Largest value id exceeds number of images, got %i (max is %i).',max(id),obj.nInput)
                end
            end

            if numel(obj.clustering) ~= numel(obj.input)
                obj = initClustering(obj);
            end

            maskState = obj.maskFlag;
            obj.maskFlag = false; % do not use a mask for the clustering

            [iStack,~,name] = getProcessedImageStack(obj,id);

            obj.maskFlag = maskState;
            
            for ii = 1:numel(id)
                obj.clustering(id(ii)) = obj.clustering(id(ii)).compute(iStack(ii));
            end
            fprintf('Loaded clustering or computed clustering using %s image.\n',name)
            
        end

        function obj = initClustering(obj)
            %INITCLUSTERING repeats the first clustering for the number of input images
            %
            %   Usage:
            %   obj = initClustering(obj) makes a vector of N copies of 
            %   the first clustering object in obj.clustering, with N the 
            %   number of input objects at obj.input.
            
            if isempty(obj.clustering)
                error('No clustering method selected. Please assign clustering object first.')
            end
            
            em = obj.clustering.empty(0,numel(obj.input));
            em(1,:) = obj.clustering(1);
            obj.clustering = em;
            
        end
        
        function Ir = getImageProcessed(obj,id,modus,varargin)
            %GETIMAGEPROCESSED get the reshaped image with all processing
            % steps applied
            %
            %   Function does allows the reduction of dimensions.
            %
            %   (1) bgCorrection if bgCorrectionFlag = true
            %   (2) mask if maskFlag = true
            %
            % Usage:
            %   Ir = getImageProcessed(obj,id,modus,varargin) gives the permuted
            %   image with number 'id' in the specified order of labels specified in
            %   varargin. Non-existing dimensions are singleton if not
            %   trailing. Dimensions that do exist and are not specified in
            %   the labels are reduced. The method is defined in modus,
            %   the available options can be found with in 'help
            %   imageStack.getReducedImage'.
            %
            %   Please note that for the mask to be applied, its dimensions
            %   should be present.

            if numel(id) > 1
                error('Please provide one id, got: %i.',numel(id))
            end

            [iStack,mask_,name] = getProcessedImageStack(obj,id);

            if ~isempty(mask_)
                Ir = iStack.applyMask(mask_,modus,varargin{:});
            else
                Ir = iStack.applyMask([],modus,varargin{:});
            end

            fprintf('Plotting using %s image.\n',name)
            
%             if obj.bgCorrectionFlag && ~isempty(obj.bgCorrection(id).I)
%                 if obj.maskFlag
%                     mask_ = obj.mask(id);
%                     Ir = obj.bgCorrection(id).applyMask(mask_,modus,varargin{:});
%                 else
%                     Ir = obj.bgCorrection(id).applyMask([],modus,varargin{:});
%                 end
%             else
%                 if obj.maskFlag
%                     mask_ = obj.mask(id);
%                     Ir = obj.input(id).applyMask(mask_,modus,varargin{:});
%                 else
%                     Ir = obj.input(id).applyMask([],modus,varargin{:});
%                 end
%             end
            
        end

        function varargout = getIbinsProcessed(obj,id,modus,binLabel,binIdx,varargin)
            % GETIBINSPROCESSED gets the I from an arbitrary number of bins with all
            % processing steps applied
            %
            %   (1) bgCorrection if bgCorrectionFlag = true
            %   (2) mask if maskFlag = true
            %
            % Usage:
            %   [ICh1,ICh2,...,IChN] = getIbinProcessed(obj,id,mask,modus,binLabel,binIdx,varargin)
            %   outputs the images ICh1,...,IChN of the {Ch1,..,ChN} channel
            %   in the cell 'binIdx'. The channels are defined by the
            %   'binLabel' variable and is usually 'c', meaning the channel
            %   values. The mask must be a imageMask object and the modus. 
            %   Dimensions that do exist and are not specified in
            %   the labels are reduced. The method is defined in modus,
            %   the available options can be found with in 'help
            %   imageStack.getReducedImage'. The id is the integer number
            %   of the data set.

            if numel(id) > 1
                error('Please provide one id, got: %i.',numel(id))
            end

            [iStack,mask_,name] = getProcessedImageStack(obj,id);

            if ~isempty(mask_)
                [varargout{1:nargout}] = iStack.getIbins(mask_,modus,binLabel,binIdx,varargin{:});
            else
                [varargout{1:nargout}] = iStack.getIbins([],modus,binLabel,binIdx,varargin{:});
            end
            fprintf('Getting bins from %s image.\n',name)
            
%             if obj.bgCorrectionFlag && ~isempty(obj.bgCorrection(id).I)
%                 if obj.maskFlag
%                     mask_ = obj.mask(id);
%                     [varargout{1:nargout}] = obj.bgCorrection(id).getIbins(mask_,modus,binLabel,binIdx,varargin{:});
%                 else
%                     [varargout{1:nargout}] = obj.bgCorrection(id).getIbins([],modus,binLabel,binIdx,varargin{:});
%                 end
%             else
%                 if obj.maskFlag
%                     mask_ = obj.mask(id);
%                     [varargout{1:nargout}] = obj.input(id).getIbins(mask_,modus,binLabel,binIdx,varargin{:});
%                 else
%                     [varargout{1:nargout}] = obj.input(id).getIbins([],modus,binLabel,binIdx,varargin{:});
%                 end
%             end
        end
        
    end

    methods (Hidden)
        function [iStack,mask_,name] = getProcessedImageStack(obj,id,excludeNames)
            %GETPROCESSEDIMAGESTACK helper function to get processed iStack
            %
            %   Usage:
            %   [iStack,mask_,name] = getProcessedImageStack(obj) gets the
            %   imageStack 'iStack', and corresponding mask 'mask_' of all
            %   images. Depending of the flags and availablity of the
            %   images, either the input, background-corrected, or
            %   normalized imageStack is returned. The order of importance
            %   is: (1) normalization, (2) background, and (3) input
            %   images. The mask is returned if available.
            %   The name can be used to generate printed output for the
            %   user, e.g. fprintf('Plotting using %s image.\n',name)
            %
            %   [iStack,mask_,name] = getProcessedImageStack(obj,id) only
            %   gets the images with index 'id'. Must be a scalar or
            %   vector.
            %   
            %   [iStack,mask_,name] = getProcessedImageStack(obj,id,excludeNames)
            %   where 'exludeNames' is a char/string or cell vector with
            %   char/string containing the name of the images that should
            %   be excluded from the order of importance. Allowed tags are:
            %   'normalization' and 'bgCorrection'.
            
            if nargin < 2 || isempty(id)
                id = 1:obj.nInput;
            end

            if nargin < 3 || isempty(excludeNames)
                excludeNames = 'noMatch';
            end
            
            if ~iscell(excludeNames)
                excludeNames = {excludeNames}; % for ismember function
            end

            % get mask
            if obj.maskFlag && numel(obj.mask) >= max(id)
                emptyFlag_mask = any(cellfun(@(x) isempty(x),{obj.mask(id).I}));
            else
                emptyFlag_mask = true;
            end
            if ~emptyFlag_mask
                mask_ = obj.mask(id);
            else
                mask_ = [];
            end

            % get imageStack
            % get in order of importance: (1) normalization, (2)
            % background, (3) input
            if obj.normalizationFlag && numel(obj.normalization) >= max(id)
                emptyFlag_norm = any(cellfun(@(x) isempty(x),{obj.normalization(id).I}));
            else
                emptyFlag_norm = true;
            end
            if obj.bgCorrectionFlag && numel(obj.bgCorrection) >= max(id)
                emptyFlag_bg = any(cellfun(@(x) isempty(x),{obj.bgCorrection(id).I}));
            else
                emptyFlag_bg = true;
            end

            if ~emptyFlag_norm && any(~ismember({'normalization'},excludeNames))
                iStack = obj.normalization(id);
                name = 'normalized';
            elseif ~emptyFlag_bg && any(~ismember({'bgCorrection'},excludeNames))
                iStack = obj.bgCorrection(id);
                name = 'background corrected';
            else
                iStack = obj.input(id);
                name = 'input';
            end

        end


    end
    
    methods (Static,Hidden)
        
        function n = getDim(imageStack,varargin)
            % GETDIM Get dimension of label of imageStack array
            %
            %   Usage:
            %   n = getDim(obj,varargin) outputs vector of dimensions
            %   matching labels in varargin. Dimensions that are not in
            %   data set have n = 1. Function sums all dimensions over the
            %   imageStack
            
            if nargin ~= 2
                error('Please provide only one label.')
            end
            
            n = 0;
            for ii = 1:numel(imageStack)
                idx = getLabelPermutation(imageStack(ii),false,varargin{:});
                if isnan(idx)
                    n = n + 1; % dimensions not in data set are set to 1
                else
                    n = n + imageStack(ii).dimSize(idx);
                end
            end
        end
        
    end
    
end

