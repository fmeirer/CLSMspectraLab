classdef (Abstract) abstract < imageStack 
    %NORMALIZATION Abstract class to compute a normalization
    
    properties (Abstract,Constant)
        name; % name of normalization computation method
    end
    
    methods (Abstract)
        obj = compute(obj,imageStack); % compute the normalization
    end
    
    methods

    end
    
end

