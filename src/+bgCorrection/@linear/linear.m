classdef linear < bgCorrection.abstract
    %LINEAR linear background correction spectra
    
    properties(Constant)
        name = 'Linear background correction';
    end
    
    properties
        c1 = 1; % first bin
        c2 = 2; % second bin
        uniformBackgroundFlag = true; % if true, compute one background for the full image, if false do per pixel
    end
    
    properties(Hidden)
        a = []; % slope dy/dx; labels are dimLabel{2:end}
        b = []; % offset y; labels are dimLabel{2:end}
    end
    
    methods
        function obj = linear(varargin)
            %LINEAR Construct an instance of this class
            %   
            %   obj = linear(varargin) with varargin Name-Value pairs,
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
            % obj = compute(obj,iStack,mask) computes and stores the
            % background correction, with mask applied
            % the imageStack object.
            
            if ~any(strcmp('c',iStack.dimLabel)) % does not contain 'c'
                error('Provided ''imageStack'' does not contain channels ''c''. ')
            end
            
            if max(obj.c1,obj.c2) > iStack.getDim('c')
                error('The coefficients are out of range and cannot be larger than % i.',iStack.getDim('c'))
            end
            
            if nargin < 3
                mask = [];
            else
                mask = varargin{1};
            end
            
            % get existing dimLabels and place 'c' first
            dimLabel = iStack.dimLabel;
            dimLabel(strcmp('c',dimLabel)) = [];
            dimLabel = ['c' dimLabel];
            I = iStack.getReshapedImage(dimLabel{:});
            dimSize = size(I);
            % reshape into [channels X pixels]
            Ir = reshape(I,iStack.getDim('c'),prod(dimSize(2:end)));
            
            % compute a and b of y = ax + b
            if obj.uniformBackgroundFlag
                % get the spectrum of the whole image
                if isempty(mask)
                    Ispectrum = iStack.getReducedImage('mean','c');
                else
                    Ispectrum = iStack.applyMask(mask,'mean','x','y','z','t','c'); % do with mask
                    Ispectrum = squeeze(mean(Ispectrum,[1 2 3 4],'omitnan'));
                end
                obj.a = (Ispectrum(obj.c2)-Ispectrum(obj.c1))/(obj.c2-obj.c1);
                obj.b = Ispectrum(obj.c1)-obj.a*obj.c1;
            else
                % compute a and b
                obj.a = (Ir(obj.c2,:)-Ir(obj.c1,:))/(obj.c2-obj.c1); % row vector
                obj.b = Ir(obj.c1,:)-obj.a*obj.c1; % row vector
            end
            
            % substract linear fit: [channels X pixels]
            if isinteger(Ir)
                if isa(class(Ir),'uint64')
                    warning('Data loss may occur during the conversion of uint64 to int64.')
                end
                % cast to int64 to allow negative values
                Isubstract =  cast(((1:iStack.getDim('c'))'.*double(obj.a) + double(obj.b)),'int64'); % substract and convert to data type Ir
                Ir = cast(Ir,'int64') - Isubstract;
            else
                Isubstract =  cast(((1:iStack.getDim('c'))'.*double(obj.a) + double(obj.b)),'like',Ir); % substract and convert to data type Ir
                Ir = Ir - Isubstract;
            end
                
            if ~obj.uniformBackgroundFlag
                % reshape to stack array
                obj.a = reshape(obj.a,dimSize(2:end)); 
                obj.b = reshape(obj.b,dimSize(2:end));
            end
            
            I = reshape(Ir,dimSize);
            
            obj = obj.addImage(I,dimLabel{:});
            obj.unit = iStack.unit;
            
        end
        
        function varargout = plotA(obj,ha,id)
            %PLOTA plots the a coeffient
            %
            % Usage:
            %
            %   plotA(obj,ha,id) plots a coefficient in 2D/3D (no
            %   timeseries).
            %   ha. Leave ha empty to use the current axis.
            %   id is an integer number of the data set to be plotted.
            %
            %   hps = obj.plotA(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj.plotA(...) returns also the handle to
            %   the axes handle of the plot.
            
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            if isscalar(obj(id).a)
                fprintf('Coefficient ''a'' is %f.\n',obj(id).a)
                return
            end
            
            % load in imageStack
            iStack = imageStack();
            iStack = iStack.addImage(obj(id).a,obj.dimLabel{2:end}); % first label is c, ignore
            iStack.unit = obj.unit;
            
            % plot
            [varargout{1:nargout}] = iStack.plotImage(ha);
            
        end
        
        function varargout = plotB(obj,ha,id)
            %PLOTB plots the b coeffient
            %
            % Usage:
            %
            %   plotA(obj,ha,id) plots b coefficient in 2D/3D (no
            %   timeseries).
            %   ha. Leave ha empty to use the current axis.
            %   id is an integer number of the data set to be plotted.
            %
            %   hps = obj.plotA(...) returns the handles to the
            %   plot.
            %
            %   [hps, ha] = obj.plotA(...) returns also the handle to
            %   the axes handle of the plot.
            
            if nargin < 2 || isempty(ha)
                ha = gca;
            end
            
            if nargin < 3
                id = 1;
            end
            
            if isscalar(obj(id).a)
                fprintf('Coefficient ''b'' is %f.\n',obj(id).b)
                return
            end
            
            % load in imageStack
            iStack = imageStack();
            iStack = iStack.addImage(obj(id).b,obj.dimLabel{2:end}); % first label is c, ignore
            iStack.unit = obj.unit;
            
            % plot
            [varargout{1:nargout}] = iStack.plotImage(ha);
            
        end
        
        
            
    end
end

