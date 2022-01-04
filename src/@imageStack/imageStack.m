classdef imageStack
    %IMAGESTACK Class for saving a single image (stack)
    %   
    
    properties(SetAccess = private)
        I; % Image (stack) that is part of a single series (time/volume)
        nDim; % number of dimensions I
        dimSize; % size of each dimension in I
        dimLabel; % labels of the dimensions, such as x,y,z,t,c
        globalMaxMin; % global maximum and minimum of stack
    end
    
    properties(SetAccess = protected)
        unit = unitProps; % instance of the unitProps class for unit handling
    end
    
    properties(SetAccess = private,Hidden)
        allowedLabels = {'x','y','z','t','c'}; % labels allowed for dimLabel
        isInverted = false; % if true, mask is inverted w.r.t. original settings
    end
    
    methods
        function obj = imageStack()
            %IMAGESTACK Construct an instance of this class
        end
        
        function obj = addImage(obj,I,varargin)
            %ADDIMAGE Adds an image to the object
            %
            %   Usage:
            %   obj = addImage(I,varargin) produces an instance of the 
            %   imageStack class for an image I. Varargin is the labels for
            %   the dimensions. Allowed labels are 'x','y','z','t','c',
            %   respectively voxel x, voxel y, voxel z, time, and channel.
            %   Each label can only be assigned once.
            
            obj.I = I;
            obj.nDim = ndims(I);
            obj.dimSize = size(I);
            obj.globalMaxMin = obj.getGlobalMaxMin();
            
            % check correct parsing labels via varargin
            if (nargin - 2) ~= obj.nDim
                error('Number of labels must be equal to number of dimensions image.\n Got %i labels and %i dimensions.',nargin-1,obj.nDim)
            end
            if any(~cellfun(@ischar,varargin)) % is not a char
                error('Labels must be a character or string.')
            end
            if numel(unique(varargin)) ~= obj.nDim
                error('Each label can only be used once.')
            end            
            if any(~ismember(varargin,obj.allowedLabels))
                error('Label(s) not allowed. Please use one of these labels: %s',strjoin(obj.allowedLabels,', '))
            end
            
            % save label
            obj.dimLabel = varargin;
            
%             % set units; initialize empty
%             u = symunit;
%             obj.unit = obj.unit.setUnit('x',[].*u.m,'y',[].*u.m,...
%                 'z',[].*u.m,'t',[].*u.s,'c',[].*u.m);
        end
        
        function obj = setUnit(obj,unitX,unitY,unitZ,unitT,unitC)
            %SETUNIT Set the units of the image I
            %
            %   Usage:
            %   obj = setUnits(obj,unitX,unitY,unitZ,unitT,unitC) sets the
            %   units in obj for respectively x,y,z,t,c. Leave empty if
            %   unknown or irrelevant. Use SI units (meter, meter, meter,
            %   seconds, meter).
            u = symunit;
            obj.unit = obj.unit.setUnit('x',unitX*u.m,'y',unitY*u.m,...
                'z',unitZ*u.m,'t',unitT*u.s,'c',unitC*u.m);
        end
        
        function bool = isVolume(obj)
            % ISVOLUME Helper function to d(etermine whether the z dimension
            % is present and the image stack can be treated as volume.
            %
            %   Usage:
            %   bool = isVolume(obj) outputs a bool wether image in object
            %   is volume.
            
            if obj.getDim('z') > 1
                bool = true;
            else
                bool = false;
            end
            
        end
        
        function bool = isTimeseries(obj)
            % ISTIMESERIES Helper function to determine whether the t dimension
            % is present and the image stack can be treated as timeseries.
            %
            %   Usage:
            %   bool = isTimeseries(obj) outputs a bool wether image in object
            %   is timeseries.
            
            if obj.getDim('t') > 1
                bool = true;
            else
                bool = false;
            end
            
        end
        
        function bool = isChannel(obj)
            % ISCHANNEL Helper function to determine whether the c dimension
            % is present and the image stack can be treated as channel set.
            %
            %   Usage:
            %   bool = isChannel(obj) outputs a bool wether image in object
            %   is channel data set.
            
            if obj.getDim('c') > 1
                bool = true;
            else
                bool = false;
            end
            
        end
        
        function voxelSize = getVoxelSize(obj)
            % GETVOXELSIZE extracts the voxel size from the unitProps
            % functionality
            %
            % Usage:
            % voxelSize = getVoxelSize(obj) outputs a row vector with the
            % voxel size of X, Y, Z. Missing values are given as NaN. 2D
            % data (identified by a isVolume = false) will produce a 1x2
            % matrix, and 3D data will produce a 1x3 matrix.
            
            if isempty(obj.unit)
                if ~obj.isVolume % 2D
                    voxelSize = nan(1,2);
                else
                    voxelSize = nan(1,3);
                end
            else
                X = obj.unit.getUnitFactor('x');
                if isempty(X); X = nan; end
                Y = obj.unit.getUnitFactor('y');
                if isempty(Y); Y = nan; end
                Z = obj.unit.getUnitFactor('z');
                if isempty(Z); Z = nan; end

                if ~obj.isVolume % 2D
                    voxelSize = [X,Y];
                else
                    voxelSize = [X,Y,Z];
                end
            end
        end
        
        function globalMaxMin = getGlobalMaxMin(obj)
            % GLOBALMAXMIN Exacts the intensity minima and maxima of whole
            % data set.
            % extract full data volume:
            % reshape into X x Y x spectrum x z-index:
            Ir = getReshapedImage(obj,'x','c','z','t','y'); % y at end to prevent tailing singleton dimension
            for ii=1:obj.getDim('z')
               sumI = squeeze(sum(Ir(:,:,ii,:,:),2)); 
               globalMaxMin(1,2) = max(sumI,[],'all');
               globalMaxMin(1,1) = min(sumI,[],'all');
            end
        end
        
        function permuteOrderNan = getLabelPermutation(obj,mustExistFlag,varargin)
            % GETLABELPERMUTATION Helper function: extracts permutation
            % order for labels
            %
            % SEE ALSO: getImage, getReshapedImage
            %
            %   Usage:
            %   permuteOrderNan = getLabelPermutation(obj,mustExistFlag,varargin)  
            %   gives the permutation order of labels specified in
            %   varargin. [index dimLabel varargin{1}, ... , index dimLabel
            %   varargin{n}] If nan, then the requested label is singleton.
            %   If mustExistFlag the function throws an error if
            %   the label does not exist. If mustExistFlag is true, the
            %   non-existing dimensions are labeled with nan in the
            %   permuteOrderNan variable.
            
            if isempty(varargin)
                error('Please provide a label')
            end
            
            if numel(unique(varargin)) ~= numel(varargin)
                error('Each label can only be used once.')
            end
            
            % determine order of permutations
            permuteOrder = nan(nargin-2,1);
            doesNotExist = zeros(nargin-2,1);
            for ii = 1:(nargin - 2)
                idx = find(strcmp(obj.dimLabel,varargin{ii}));
                if numel(idx) < 1; doesNotExist(ii) = true; idx = nan; end
                permuteOrder(ii) = idx;
            end
            permuteOrderNan = permuteOrder;
            
            % raise error if label does not exist
            if mustExistFlag && any(doesNotExist)
                error('Label(s) %s do(es) not exist in image.',strjoin(varargin(find(doesNotExist)),', '));
            end
            
        end
        
        function Ir = getReshapedImage(obj,varargin)
            % GETRESHAPEDIMAGE reshapes image based on labels.
            %   Function does not allow the reduction of dimensions.
            %
            %   Usage:
            %   Ir = getReshaedImage(obj,varargin) gives the permuted
            %   image in the specified order of labels specified in
            %   varargin. Non-existing dimensions are singleton if not
            %   trailing.
            
            if any(~ismember(varargin,obj.allowedLabels))
                error('Label(s) not allowed. Please use one of these labels: %s',strjoin(obj.allowedLabels,', '))
            end
            
            permuteOrderNaN = getLabelPermutation(obj,false,varargin{:});
            if numel(permuteOrderNaN) ~= numel(varargin)
                error('Dimensions of the image are reduced.\nPlease use the ''getImage'' function to specify how the trailing dimension should be treated.')
            end
            % permute image
            nNan = numel(permuteOrderNaN(isnan(permuteOrderNaN)));
            nNotNan = numel(permuteOrderNaN)-nNan;
            permuteOrderNaN(isnan(permuteOrderNaN)) = nNotNan+1:nNotNan+nNan; % fill nan with dimensions outside range image
            Ir = permute(obj.I,permuteOrderNaN); % dimensions outside range image are made singleton
            
        end
            
        function Ir = getReducedImage(obj,modus,varargin)
            % GETREDUCEDIMAGE reshapes image based on labels.
            %   Function does allows the reduction of dimensions.
            %
            %   Usage:
            %   Ir = getReshaedImage(obj,modus,varargin) gives the permuted
            %   image in the specified order of labels specified in
            %   varargin. Non-existing dimensions are singleton if not
            %   trailing. Dimensions that do exist and are not specified in
            %   the labels are reduced. The method is defined in modus,
            %   options are:
            %       - 'sum': dimensions are summed together
            %       - 'mean': the mean of the dimensions is taken
            %       - 'reshaped': the getReshapedImage function is called
            
            if any(~ismember(varargin,obj.allowedLabels))
                error('Label(s) not allowed. Please use one of these labels: %s',strjoin(obj.allowedLabels,', '))
            end
            
            permuteOrderNaN = getLabelPermutation(obj,false,varargin{:});
            
            % --- permute image
            % find dimensions in data set that have to be removed; these 
            % are either not present data set (nan in permuteOrderNaN), or
            % present in data set, but not present requested labels
            % (obj.nDim).
            trailIdx = find(~ismember(1:obj.nDim,permuteOrderNaN)); 
            nNan = numel(permuteOrderNaN(isnan(permuteOrderNaN)));
            % fill nan with dimensions outside range image
            permuteOrderNaN(isnan(permuteOrderNaN)) = obj.nDim+1:obj.nDim+nNan; 
            % dimensions outside range image are made singleton during
            % permutation
            Ir = permute(obj.I,[permuteOrderNaN' trailIdx]); 
            
            if ~isempty(trailIdx)
                idxDimThrow = numel(varargin)+1:numel(varargin)+numel(trailIdx); % indices of to be removed dimensions
                switch modus
                    case 'sum'
                        Ir = sum(Ir,idxDimThrow);
                    case 'mean'
                        Ir = mean(Ir,idxDimThrow);
                    case 'reshaped'
                        Ir = obj.getReshapedImage(varargin{:}); % call function
                    otherwise
                        error('The modus ''%s'' is not recognised by the function ''getImage''.',modus)
                end
            end
            
        end
        
        function [varargout] = getIbins(obj,mask,modus,binLabel,binIdx,varargin)
            % GETIBINS gets the I from an arbitrary number of bins
            %
            % Usage:
            %   [ICh1,ICh2,...,IChN] = getIbins(obj,mask,modus,binLabel,binIdx,varargin)
            %   outputs the images ICh1,...,IChN of the {Ch1,..,ChN} channel
            %   in the cell 'binIdx'. The channels are defined by the
            %   'binLabel' variable and is usually 'c', meaning the channel
            %   values. The mask must be a imageMask object and the modus. 
            %   Dimensions that do exist and are not specified in
            %   the labels are reduced. The method is defined in modus,
            %   the available options can be found with in 'help
            %   imageStack.getReducedImage'.
            
            if nargin < 2
                mask = [];
            end
            
            % -- check labels
            
            if any(~ismember(varargin,obj.allowedLabels))
                error('Label(s) not allowed. Please use one of these labels: %s',strjoin(obj.allowedLabels,', '))
            end
            
            if iscell(binLabel)
                if numel(binLabel) ~= 1
                    error('Please provide one bin label.')
                end
                binLabel = binLabel{1}; % unpack cell
            end
            
            if any(~ismember(binLabel,obj.allowedLabels))
                error('Bin label is not allowed. Please use one of these labels: %s',strjoin(obj.allowedLabels,', '))
            end
            
            
            if ~ismember(binLabel,varargin)
                warning('Bin label does not exist in the requested dimensions. Adding binLabel as last dimension.')
                varargin{end+1} = binLabel; % add dimension
            end
            
            % --- check binIdx
            if ~iscell(binIdx)
                binIdx = {binIdx}; % store in cell
            end
            
            % --- apply mask

            if ~isempty(mask) % request masked image including binLabel
                Ir = applyMask(obj,mask,modus,varargin{:});
            else
                Ir = obj.getReducedImage(modus,varargin{:});
            end
            
            % --- take requested bins
            
            binLabelIdx = ismember(varargin,binLabel); % index of the binLabel in the requested dimensions (varargin)
            
            varargout = cell(1,numel(binIdx)); % number of bins requested
            for ii = 1:numel(binIdx)
                thisChannelIdx = binIdx{ii};
                if ~ismember(thisChannelIdx,1:getDim(obj,binLabel))
                    error('Channel index %i with a value of %i is not a valid index.',ii,thisChannelIdx);
                end
                % make array with desired indices Ir array
                takeIrIdx = cell(1,ndims(Ir));
                for jj = 1:ndims(Ir)
                    takeIrIdx{jj} = 1:size(Ir,jj);
                end
                takeIrIdx{binLabelIdx} = thisChannelIdx; % set the bins in binLabel dimension
                varargout{ii} = Ir(takeIrIdx{:});
            end
            
        end
        
        function n = getDim(obj,varargin)
            % GETDIM Get dimension of label
            %
            %   Usage:
            %   n = getDim(obj,varargin) outputs vector of dimensions
            %   matching labels in varargin. Dimensions that are not in
            %   data set have n = 1.
            
            if nargin > 2
                error('Please provide only one label.')
            end
            
            idx = getLabelPermutation(obj,false,varargin{:});
            if isnan(idx)
                n = 1; % dimensions not in data set are set to 1
            else
                n = obj.dimSize(idx);
            end
            
        end
        
        function V = getMeanChannel(obj)
            % GETMEANCHANNEL gets the value per channel for the full data
            % set
            %
            %   Usage:
            %   V = getMeanChannel(obj) returns a vector with the mean
            %   values per channel
            
            I_ = obj.getReshapedImage('c','z','t','x','y');
            nChannels = obj.getDim('c');
            V = nan(nChannels,1);
            for ii = 1:nChannels
                V(ii) = mean(I_(ii,:,:,:,:),'all');
            end
            
        end
        

        function I = applyMask(obj,mask,modus,varargin)
            %APPLYMASK Applies the mask in obj to the image stack
            %
            %   Usage:
            %   I = applyMask(obj,mask,modus,varargin) outputs an array
            %   with dimensions specified in varargin. This is the masked
            %   image of obj with mask.
            %   The options for the variable 'modus' can be found in the
            %   function imageStack.getReducedImage
            
            
            I = obj.getReducedImage(modus,varargin{:}); % was: obj.getReshapedImage(varargin{:});
            if isempty(mask) % check if mask object is available
                return
            end
            
            Imask = mask.getReshapedImage(varargin{:});
            
            % An alternate method to multiplication channel by channel.
            % Mask the image using bsxfun() function
            % https://nl.mathworks.com/matlabcentral/answers/38547-masking-out-image-area-using-binary-mask
            if ~isempty(Imask) % check if mask is loaded and computed
                if any(~ismember(mask.dimLabel,varargin)) % mask dimension is not requested
                    % in the fix, the missing mask dimenions should be
                    % added and averaged/summed in a second step in this
                    % function
                    error('Mask dimension must be requested, otherwise masking fails. Can be fixed in a future release.')
                end
                I = bsxfun(@times, I, cast(Imask, 'like', I));
            end
            
        end
        
        function iS = applyMaskImageStack(obj,mask,modus,varargin)
            %APPLYMASKIMAGESTACK Applies the mask in obj to the image stack
            %
            %   Usage:
            %   I = applyMask(obj,mask,modus,varargin) outputs an 
            %   imageStack with dimensions specified in varargin. This is 
            %   the masked image of obj with mask.
            %   The options for the variable 'modus' can be found in the
            %   function imageStack.getReducedImage
            
            
            I_ = obj.getReducedImage(modus,varargin{:}); % was: obj.getReshapedImage(varargin{:});
            Imask = mask.getReshapedImage(varargin{:});
            
            % An alternate method to multiplication channel by channel.
            % Mask the image using bsxfun() function
            % https://nl.mathworks.com/matlabcentral/answers/38547-masking-out-image-area-using-binary-mask
            I_ = bsxfun(@times, I_, cast(Imask, 'like', I_));
            
            % write to image stack
            iS = imageStack();
            iS = iS.addImage(I_,varargin{:});
            iS.unit = obj.unit;
            
        end
        
        function varargout = plotImage(obj,ha,mask)
            % PLOTIMAGE plots the image in XY or XYZ.
            %
            %   Usage:
            %   plotImage(obj) plots the image.
            %
            %   plotImage(obj,ha) plots the image or volume in the axis 
            %   specified by the axis handle ha. All non-spatial dimensions are
            %   summed.
            %
            %   plotImage(obj,ha,mask) uses the imageMask object in 
            %   vargargin.
            %
            %   hps = obj(ii).plotImage(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj(ii).plotImage(...) returns also the handle to
            %   the axes handle of the plot.


            if nargin < 2 || isempty(ha)
                ha = gca;
            end 

            if nargin < 3
                mask = [];
            end

            if ~isempty(mask)
                Ir = applyMask(obj,mask,'sum','x','y','z');
            else
                Ir = obj.getReducedImage('sum','x','y','z');
            end


            if ~obj.isVolume
                hps = imagesc(Ir); 
                axis image;
            else
                voxelSize = obj.getVoxelSize;
                voxelSize = voxelSize./voxelSize(1); % normalize to X
                voxelSize(isnan(voxelSize)) = 1; % set al unknown values to 1
    %             input('Displaying volume. Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
                volumeViewer(Ir,'ScaleFactors',voxelSize);
                hps = [];
            end

            % Output
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end

        end
        
        function obj = invertImage(obj)
            %INVERTIMAGE inverts the values for an image containing 0 and 1
            %
            %   Usage:
            %   obj = invertImage(obj) converts 0 -> 1 and vice versa in
            %   image.
            
            if any(~ismember(unique(obj.I(:)),[0 1]))
                error('Image contains values apart from 0 and 1. Cannot invert.')
            end
            
            % invert
            obj.I = ~obj.I;
            
            % administration, change isInverted property
            if obj.isInverted
                obj.isInverted = false;
            else
                obj.isInverted = true;
            end
            
        end
        
    end
    
    methods (Static)
        function [imageStackArray,varargout] = import(filepath)
            %IMPORT Loads image and construct an instance of this class
            %   imageStackArray = import(filepath) accepts 
            %   a string with the file path to the image file or a cell
            %   with the file path to multiple files. If empty a file 
            %   selector is opened. Returns a column vector of imageStack
            %   objects.
            %
            %   [imageStackArray,nImageStack] = import(filepath) returns
            %   the size of the image stack array.
            %
            %   WARNING: Not tested with XYZTC data
            
            if nargin < 1
                filepath = [];
            end
            
            if ischar(filepath)
                % single file; string
                dataAll = bfopen(filepath);
            elseif iscell(filepath)
                % multiple files; cell string
                  nFiles = numel(filepath);      
                  dataAll={};
                  for ii=1:nFiles
                      data = bfopen(filepath{ii});           
                      dataAll = [dataAll; data];
                  end
            else
                % file selection dialog
                  [files, path, ~] = uigetfile(bfGetFileExtensions, 'Choose a file to open','Multiselect','on');
                  nFiles = size(files,2);      
                  dataAll={};
                  if iscell(files) % user selected multiple files
                      for ii=1:nFiles
                          data = bfopen(fullfile(path, files{ii}));             
                          dataAll = [dataAll; data];
                      end
                  else % user selected single file
                      dataAll = bfopen(fullfile(path, files));
                  end
            end
            
            % now we have rows with the files and columns with the file
            % properties and data. nFiles tells us the number of files
            % loaded
            nSeries =  size(dataAll,1);
            
            % put data in imageStack object
            imageStackArray = imageStack.empty(0,nSeries);
            for ii = 1:nSeries
                omeMeta = dataAll{ii, 4};
                imageX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
                imageY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
                imageZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
                imageT = omeMeta.getPixelsSizeT(0).getValue(); % number of T slices
                imageC = omeMeta.getPixelsSizeC(0).getValue(); % number of C slices
                imageData = cat(3,dataAll{ii,1}{:,1}); % X x Y x [X,Z]
                % reshape to X x Y x C x Z x T, I don't know what the right
                % order is in all different dimensionalities (e.g. XYCT
                % data set).
                imageData = squeeze(reshape(imageData,imageY,imageX,imageC,imageZ,imageT)); 
                if imageT > 1; warning('Timeseries detected. Import timeseries in not tested.'); end
                
                dimLabels = {'x','y','c','z','t'};
                if imageC == 1; dimLabels(strcmp(dimLabels,'c')) = []; end % remove channel label if singleton dimension
                if imageZ == 1; dimLabels(strcmp(dimLabels,'z')) = []; end % remove z label if singleton dimension
                if imageT == 1; dimLabels(strcmp(dimLabels,'t')) = []; end % remove t label if singleton dimension
                imageStackArray(ii) = imageStack();
                imageStackArray(ii) = imageStackArray(ii).addImage(imageData,dimLabels{:});
                [voxelSizeX,voxelSizeY,voxelSizeZ] = getVoxelSizeMetaData(omeMeta);
                imageStackArray(ii) = imageStackArray(ii).setUnit(voxelSizeY,voxelSizeX,voxelSizeZ,[],[]);
            end
            
            function [voxelSizeX,voxelSizeY,voxelSizeZ] = getVoxelSizeMetaData(omeMeta)

                try
                    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.METER); % in m
                    voxelSizeX = double(voxelSizeX);
                catch
                    disp('WARNING: Voxel size X could not be read from file!');
                    voxelSizeX = [];    
                end
                try
                    voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.METER); % in m
                    voxelSizeY = double(voxelSizeY);
                catch
                    disp('WARNING: Voxel size Y could not be read from file!');
                    voxelSizeY = [];
                end
                try
                    voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.METER); % in m
                    voxelSizeZ = double(voxelSizeZ);
                catch
                    disp('WARNING: Voxel size Z could not be read from file!');
                    voxelSizeZ = [];
                end

            end
            
            if nargout > 1
                varargout{1} = nSeries;
            end
        end
        
        function imageStackArray = importTifStack(filepath,stackLabel)
            %IMPORTTIFSTACK Loads image and construct an instance of this class
            %   imageStackArray = importTifStack(filepath,stackLabel): accepts 
            %   a string with the file path to the image file or a cell
            %   with the file path to multiple files.
            %   stackLabel is a character with the label of the third
            %   dimension, which can be either 'z','t', or 'c'.
            %   Returns a column vector of imageStack objects.
            
            if ischar(filepath)
                filepath = {filepath};
            end
            
            imageStackArray = imageStack.empty(numel(filepath),0);
            for ii = 1:numel(filepath)
            warning('off','all') % Suppress all the tiff warnings
            tstack  = Tiff(filepath{ii});
            [I,J] = size(tstack.read());
            K = length(imfinfo(filepath{ii}));
            image = zeros(I,J,K);
            image(:,:,1)  = tstack.read();
            for n = 2:K
                tstack.nextDirectory()
                image(:,:,n) = tstack.read();
            end
            warning('on','all')
                imageStackArray(ii) = imageStack();
                imageStackArray(ii) = imageStackArray(ii).addImage(image,'x','y',stackLabel);
                imageStackArray(ii) = imageStackArray(ii).setUnit([],[],[],[],[]); % unknown
            end
            
        end
        
    end
end

