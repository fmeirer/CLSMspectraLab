classdef sum2one < pixelNormalization.abstract
    %NORMALIZATION sum-of-all-channels normalization.
    %
    % Absolute value of negative pixel values are taken to prevent
    % sum of pixel cose to zero.
    
    properties(Constant)
        name = 'sum2one';
    end
    
    properties
        dataType = 'integer'; % use: 'integer' for an unsigned 16 bit (normalizes to bit range) or 'float' for a double
    end
    
    properties(Hidden)

    end
    
    methods
        function obj = sum2one(varargin)
            %SUM2ONE Construct an instance of this class
            %   
            %   obj = sum2one(varargin) with varargin Name-Value pairs,
            %   e.g. if we set the property 'prop' to 5: 
            %   obj = sum2one('prop',5).
            
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
            %COMPUTE computes the domains with similar spectra and stores a
            %normalized imagestack
            %
            % Usage:
            % obj = compute(obj,iStack) computes and stores the normalized 
            % imageStack object.
            
            if ~any(strcmp('c',iStack.dimLabel)) % does not contain 'c'
                error('Provided ''imageStack'' does not contain channels ''c''. ')
            end
            
            % reshape to get 'c' channel as first dimension
            dimLabel_input = iStack.dimLabel;
            dimLabel_input(ismember(dimLabel_input,{'c'})) = [];
            dimLabel_input = ['c' dimLabel_input];
            I = iStack.getReshapedImage(dimLabel_input{:});

            % get c x Npixels matrix
            sz = size(I);
            Ir = reshape(I,iStack.getDim('c'),prod(iStack.dimSize)/iStack.getDim('c'));
            normFactor = sum(abs(Ir),1); % outputs double
            % take abs to deal with negative values, leading to normalization factor close to zero
            normFactor(normFactor == 0) = 1; % skip zeros

            switch obj.dataType
                case 'integer'
                    Ir = uint16(double(Ir)./normFactor.*65535);
                case 'float'
                    Ir = double(Ir)./normFactor;
                otherwise
                    error('dataType not ''%s'' is not recognised. Use ''integer'' or ''float''.',obj.dataType)
            end
            % get original shape back
            I = reshape(Ir,sz);
            
            obj = obj.addImage(I,dimLabel_input{:});
            obj.unit = iStack.unit;
            
        end
            
    end
end

%             if sum(Ir_sum == 0) > 0
%                 N = Ir_sum == 0; % gives Inf
%                 Ir_sum(N) = 1; % devide by one or negative value = ignore
%                 warning('Found %i pixels in which the normalization factor is 0. Ignoring Pixels.\nThis warning can be ignored if pixel values >> 1.\nOtherwise perform normalization without background correction.',sum(N))
%             end

