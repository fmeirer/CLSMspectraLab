classdef emptyMask < imageMask.abstract
    %CUSTOM Custom loaded mask
    
    properties(Constant)
        name = 'Custom mask';
    end
    
    methods
        function obj = emptyMask(varargin)
            %CUSTOMMASK Construct an instance of this class
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
            %COMPUTE stores a mask without modification
            %
            % Usage:
            % obj = compute(obj,imageStack) stores the mask of imageStack,
            % which is an imageStack object.
            
            % get labels without C
            dimLabel = imageStack.dimLabel;
            dimSize = imageStack.dimSize;
            idx_c = strcmp(dimLabel,'c');
            dimSize(idx_c) = [];
            dimLabel(idx_c) = [];
            
            % create image with ones, all is included in mask
            I = ones(dimSize);
            obj = obj.addImage(I,dimLabel{:});
            obj.unit = imageStack.unit;
        end

    end
end

