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
    end
    
    methods
        function obj = CLSMspectraLab(filepath)
            %CLSMSPECTRALAB Construct an instance of this class
            %   obj = CLSMspectraLab(filepath): accepts 
            %   a string with the file path to the image file or a cell
            %   with the file path to multiple files. If empty a file 
            %   selector is opened. Returns an instance of the  
            %   CLSMspectraLab class.
            
            if nargin < 1
                filepath = [];
            end
            
            % write to object
            [obj.input,obj.nInput] = imageStack.import(filepath);
            
        end
        
        function obj = computeBgCorrection(obj)
            %COMPUTEBGCORRECTION computes the background correction for the
            % input image
            %
            %   Usage:
            %   obj = computeBgCorrection(obj) computes the background 
            %   correction following the 1st bgCorrection object in 
            %   obj.bgCorrection for all different input images.
            
            if numel(obj.bgCorrection) ~= numel(obj.input)
                obj = initBgCorrection(obj);
            end
            
            for ii = 1:numel(obj.bgCorrection)
                obj.bgCorrection(ii) = obj.bgCorrection(ii).compute(obj.input(ii));
            end
            
        end
        
        function obj = initBgCorrection(obj)
            %INITBGCORRECTION repeats the first bgCorrection for the number of input images
            %
            %   Usage:
            %   obj = initBgCorrection(obj) makes a vector of N copies of 
            %   the first bgCorrection object in obj.bgCorrection, with N the 
            %   number of input objects at obj.input.
            
            if isempty(obj.bgCorrection)
                error('No mask selected. Please assign mask first.')
            end
            
            em = obj.bgCorrection.empty(0,numel(obj.input));
            em(1,:) = obj.bgCorrection(1);
            obj.bgCorrection = em;
            
        end
        
        function obj = computeMask(obj)
            %COMPUTEMASK computes the mask for the input image
            %
            %   Usage:
            %   obj = computeMask(obj) computes the mask following the
            %   1st inputMask object in obj.mask for all different input
            %   images.
            
            if numel(obj.mask) ~= numel(obj.input)
                obj = initMask(obj);
            end
            
            if ~obj.bgCorrectionFlag || isempty(obj.bgCorrection.I) % without background correction
                for ii = 1:numel(obj.mask)
                    obj.mask(ii) = obj.mask(ii).compute(obj.input(ii));
                end
                fprintf('Computed mask without background correction.\n')
            else
                for ii = 1:numel(obj.mask)
                    obj.mask(ii) = obj.mask(ii).compute(obj.bgCorrection(ii));
                end
                fprintf('Computed mask with background correction.\n')
            end
            
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
        
        function obj = invertMask(obj)
            %INVERTMASK inverts the values for an image containing 0 and 1
            %
            %   Usage:
            %   obj = invertMask(obj) converts 0 -> 1 and vice versa in
            %   mask.
            
            for ii = 1:numel(obj.mask)
                obj.mask(ii) = obj.mask(ii).invertImage;
            end
            
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
            
            if obj.bgCorrectionFlag && ~isempty(obj.bgCorrection(id).I)
                if obj.maskFlag
                    mask_ = obj.mask(id);
                    Ir = obj.bgCorrection(id).applyMask(mask_,modus,varargin{:});
                else
                    Ir = obj.bgCorrection(id).applyMask([],modus,varargin{:});
                end
            else
                if obj.maskFlag
                    mask_ = obj.mask(id);
                    Ir = obj.input(id).applyMask(mask_,modus,varargin{:});
                else
                    Ir = obj.input(id).applyMask([],modus,varargin{:});
                end
            end
            
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
            
            if obj.bgCorrectionFlag && ~isempty(obj.bgCorrection(id).I)
                if obj.maskFlag
                    mask_ = obj.mask(id);
                    [varargout{1:nargout}] = obj.bgCorrection(id).getIbins(mask_,modus,binLabel,binIdx,varargin{:});
                else
                    [varargout{1:nargout}] = obj.bgCorrection(id).getIbins([],modus,binLabel,binIdx,varargin{:});
                end
            else
                if obj.maskFlag
                    mask_ = obj.mask(id);
                    [varargout{1:nargout}] = obj.input(id).getIbins(mask_,modus,binLabel,binIdx,varargin{:});
                else
                    [varargout{1:nargout}] = obj.input(id).getIbins([],modus,binLabel,binIdx,varargin{:});
                end
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

