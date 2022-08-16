classdef reference < bgCorrection.abstract
    %LINEAR linear background correction spectra
    
    properties(Constant)
        name = 'reference';
    end
    
    properties
        Iref% = imageStack()
        uniformBackgroundFlag = true; % if true, compute one background for the full image, if false do per pixel
        binFactor = 1; % bin factor of the reference spectrum is diectly applied after import.
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
        
        function obj = compute(obj,iStack,varargin)
            %COMPUTE computes the domains with similar spectra and stores a mask
            %
            % Usage:
            % obj = compute(obj,iStack) computes and stores the mask of
            % the imageStack object.

            if isempty(obj.Iref) % does not contain 'c'
                error('No reference image was loaded, please load one before doing the reference background correction')
            end

            % perform binning
            if obj.binFactor ~= 1
                Iref_ = obj.Iref.binImageFirst2Dims(obj.binFactor,'x','y','z','t','c');
            else
                Iref_ = obj.Iref;
            end
            
            if any(~ismember(iStack.dimLabel,{'c','x','y'})) % does not contain 'c'
                error('Provided ''imageStack'' does not contain channels ''c'', ''x'', and ''y''.')
            end
            
            if iStack.getDim('c') ~= Iref_.getDim('c')
                error('The input imageStack does not have the same number of channels as the reference')
            end
            
            if ~obj.uniformBackgroundFlag
                if iStack.getDim('x') ~= Iref_.getDim('x') || iStack.getDim('y') ~= Iref_.getDim('y')
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
                Iref_spec = Iref_.getReducedImage('mean','c');
                
                % substract reference image spectrum 
                if isinteger(Ir)
                    Ir = Ir - cast(Iref_spec,'like',Ir); % do not allow negative values, assuming input is unsigned integer;

                    % code that allows for negative values. At this point,
                    % I do not see the advantage of this approach.
%                     if isa(class(Ir),'uint64')
%                         warning('Data loss may occur during the conversion of uint64 to int64.')
%                     end
%                     Ir = cast(Ir,'int64') - cast(Iref_spec,'int64');
                else
                    % Isubstract =  cast(((1:iStack.getDim('c'))'.*double(obj.a) + double(obj.b)),'like',Ir); % substract and convert to data type Ir
                    Ir = Ir - repmat(Iref_spec,[1,size(Ir,2)]);
                end
                
            else
                
                % get existing dimLabels and place 'c' first for reference image
                dimLabelreference = Iref_.dimLabel;
                dimLabelreference(strcmp('c',dimLabelreference)) = [];
                dimLabelreference = ['c' dimLabelreference];
                Iref_im = Iref_.getReshapedImage(dimLabelreference{:});
                % reshape into [channels X Y rest]
%                 Irreference = reshape(Iref_im,obj.Iref.getDim('c'),[]);
                
                %if ~obj.uniformBackgroundFlag
                Ir = (Ir - Iref_im);
            end
            
            obj = obj.addImage(Ir,dimLabel_input{:});
            obj.unit = iStack.unit;
        end
    end
    
end


