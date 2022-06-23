classdef reference < bgCorrection.abstract
    %LINEAR linear background correction spectra
    
    properties(Constant)
        name = 'reference';
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
            %   obj = reference(varargin) with varargin Name-Value pairs,
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
            
            if any(~ismember(iStack.dimLabel,{'c','x','y'})) % does not contain 'c'
                error('Provided ''imageStack'' does not contain channels ''c'', ''x'', and ''y''.')
            end
            
            if isempty(obj.Iref) % does not contain 'c'
                error('No reference image was loaded, please load one before doing the reference background correction')
            end
            
            if iStack.getDim('c') ~= obj.Iref.getDim('c')
                error('The input imageStack does not have the same number of channels as the reference')
            end
            
            if ~obj.uniformBackgroundFlag
                if iStack.getDim('x') ~= obj.Iref.getDim('x') || iStack.getDim('y') ~= obj.Iref.getDim('y')
                    error('Dimensions input imageStack and reference do not match dimensions in x and y. Use uniformBackgroundFlag = false instead.')
                end
            end
            
            % get existing dimLabels and place 'c', 'x', 'y' first
            dimLabel_input = iStack.dimLabel;
            dimLabel_input(ismember(dimLabel_input,{'c','x','y'})) = [];
            dimLabel_input = ['c','x','y' dimLabel_input];
            Ir = iStack.getReshapedImage(dimLabel_input{:});
%             dimSize = size(I);
%             % reshape into [channels X pixels]
%             Ir = reshape(I,iStack.getDim('c'),prod(dimSize(2:end)));
            
            if obj.uniformBackgroundFlag
                % get the spectrum of the whole image
                Iref_spec = obj.Iref.getReducedImage('mean','c');
                
                % substract reference image spectrum 
                if isinteger(Ir)
                    if isa(class(Ir),'uint64')
                        warning('Data loss may occur during the conversion of uint64 to int64.')
                    end
                    Ir = cast(Ir,'int64') - cast(Iref_spec,'int64');
                else
                    % Isubstract =  cast(((1:iStack.getDim('c'))'.*double(obj.a) + double(obj.b)),'like',Ir); % substract and convert to data type Ir
                    Ir = Ir - repmat(Iref_spec,[1,size(Ir,2)]);
                end
                
            else
                
                % get existing dimLabels and place 'c' first for reference image
                dimLabelreference = obj.Iref.dimLabel;
                dimLabelreference(strcmp('c',dimLabelreference)) = [];
                dimLabelreference = ['c' dimLabelreference];
                Iref_im = obj.Iref.getReshapedImage(dimLabelreference{:});
                % reshape into [channels X Y rest]
                Irreference = reshape(Iref_im,obj.Iref.getDim('c'),prod(dimSize(2:numel(dimSize))));
                
                %if ~obj.uniformBackgroundFlag
                Ir = (Ir - Irreference);
            end
            
            obj = obj.addImage(Ir,dimLabel_input{:});
            obj.unit = iStack.unit;
        end
    end
    
end


