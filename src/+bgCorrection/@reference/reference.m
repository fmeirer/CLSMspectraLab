classdef reference < bgCorrection.abstract
    %LINEAR linear background correction spectra
    
    properties(Constant)
        name = 'reference background correction';
    end
    
    properties
        referenceimage% = imageStack()
        uniformBackgroundFlag = true; % if true, compute one background for the full image, if false do per pixel
    end
    
    properties(Hidden)
    end
    
    methods
        function obj = reference(varargin)
            % reference Construct an instance of this class
            %
            %   obj = freference(varargin) with varargin Name-Value pairs,
            %   e.g. if we set the property 'prop' to 5:
            %   obj = linear('prop',5).
            
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
        
        function obj = compute(obj,iStack)
            %COMPUTE computes the domains with similar spectra and stores a mask
            %
            % Usage:
            % obj = compute(obj,iStack) computes and stores the mask of
            % the imageStack object.
            
            %[file,path] = uigetfile([],'Select Flatfield image')
            
            if ~any(strcmp('c',iStack.dimLabel)) % does not contain 'c'
                error('Provided ''imageStack'' does not contain channels ''c''. ')
            end
            
            
            if isempty(obj.referenceimage) % does not contain 'c'
                error('No reference image was loaded, please load one before doing the reference background correction')
            end
            
            % get existing dimLabels and place 'c' first
            dimLabel = iStack.dimLabel;
            dimLabel(strcmp('c',dimLabel)) = [];
            dimLabel = ['c' dimLabel];
            I = iStack.getReshapedImage(dimLabel{:});
            dimSize = size(I);
            % reshape into [channels X pixels]
            Ir = reshape(I,iStack.getDim('c'),prod(dimSize(2:end)));
            
            % get existing dimLabels and place 'c' first for reference image
            iStackreference=obj.referenceimage;
            dimLabelreference = iStackreference.dimLabel;
            dimLabelreference(strcmp('c',dimLabelreference)) = [];
            dimLabelreference = ['c' dimLabelreference];
            Ireference = iStackreference.getReshapedImage(dimLabelreference{:});
            dimSizereference = size(Ireference);
            % reshape into [channels X pixels]
            Irreference = reshape(Ireference,iStackreference.getDim('c'),prod(dimSize(2:end)));
            
            % compute a and b of y = ax + b
            if obj.uniformBackgroundFlag
                % get the spectrum of the whole image
                Ispectrum = iStack.getReducedImage('mean','c');
                Ispectrumreference = iStackreference.getReducedImage('mean','c');
                
                
                
                
                % substract reference image spectrum 
                if isinteger(Ir)
                    if isa(class(Ir),'uint64')
                        warning('Data loss may occur during the conversion of uint64 to int64.')
                    end
                    Ispectrumreference=cast(Ispectrumreference,'int64');
                    Ir=cast(Ir,'int64');
                    %try
                    Ir = Ir-  repmat(Ispectrumreference,[1,size(Ir,2)]);
                   % catch
                    %end
                else
                    % Isubstract =  cast(((1:iStack.getDim('c'))'.*double(obj.a) + double(obj.b)),'like',Ir); % substract and convert to data type Ir
                    Ir = Ir - repmat(Ispectrumreference,[1,size(Ir,2)]);
                end
                
            else
                
                %if ~obj.uniformBackgroundFlag
                Ir = (Ir - Irreference);
            end
            
            I = reshape(Ir,dimSize);
            
            obj = obj.addImage(I,dimLabel{:});
            obj.unit = iStack.unit;
        end
    end
    
end


